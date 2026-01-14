#![allow(unused_imports, unused_variables)]
use anyhow::{anyhow, Result as AnyResult};
use crossbeam::channel;
use minimap2::{ Aligner, Built};
use minimap2_sys::*;
use rust_htslib::bam::{ 
    self,
    record::Aux, record::CigarStringView, 
    record::Cigar, record::CigarString,
    Read, Reader, Record, HeaderView, 
    Header, header::HeaderRecord,
    Writer, ext::BamRecordExtensions
};
use std::collections::HashMap;
use std::sync::Arc;
use std::sync::mpsc;
use std::thread;
use std::ffi::{CStr, CString};
use std::io::{BufRead, BufReader};
use std::path::Path;
use rayon::prelude::*;


// fn build_bam_header_from_fai(reference_fasta: &str) -> AnyResult<Header> {
//     let fai_path = format!("{reference_fasta}.fai");
//     if !Path::new(&fai_path).exists() {
//         return Err(anyhow!(
//             "FASTA index not found: {fai_path}. Please run: samtools faidx {reference_fasta}"
//         ));
//     }

//     let f = common_reader(&fai_path);

//     let mut header = Header::new();

//     for (lineno, line) in f.lines().enumerate() {
//         let line = line?;
//         if line.trim().is_empty() {
//             continue;
//         }
//         let cols: Vec<&str> = line.split('\t').collect();
//         if cols.len() < 2 {
//             return Err(anyhow!("Bad .fai line {}: {}", lineno + 1, line));
//         }
//         let name = cols[0];
//         let len: i64 = cols[1]
//             .parse()
//             .map_err(|_| anyhow!("Bad length in .fai line {}: {}", lineno + 1, line))?;

//         let mut sq = HeaderRecord::new(b"SQ");
//         sq.push_tag(b"SN", name);
//         sq.push_tag(b"LN", &len.to_string());
//         sq.push_tag(b"UR", reference_fasta);

//         header.push_record(&sq);
//     }

  
//     let mut pg = HeaderRecord::new(b"PG");
//     pg.push_tag(b"ID", "bamaligner");
//     pg.push_tag(b"PN", "bamaligner");
//     header.push_record(&pg);

//     Ok(header)
// }

// fn build_tid_map_from_fai(reference_fasta: &str) -> AnyResult<HashMap<String, i32>> {
//     let fai_path = format!("{reference_fasta}.fai");
//     let r = common_reader(&fai_path);
//     let r = BufReader::new(r);

//     let mut map = HashMap::new();
//     let mut tid: i32 = 0;
//     for line in r.lines() {
//         let line = line?;
//         if line.trim().is_empty() { continue; }
//         let name = line.split('\t').next().unwrap();
//         map.insert(name.to_string(), tid);
//         tid += 1;
//     }
//     Ok(map)
// }

fn parse_cigar_string(cigar: &str) -> AnyResult<CigarString> {
    let mut ops: Vec<Cigar> = Vec::new();
    let mut num: u32 = 0;

    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            num = num * 10 + (ch as u8 - b'0') as u32;
            continue;
        }
        if num == 0 {
            return Err(anyhow!("Invalid CIGAR string: {cigar}"));
        }
        let op = match ch {
            'M' => Cigar::Match(num),
            'I' => Cigar::Ins(num),
            'D' => Cigar::Del(num),
            'N' => Cigar::RefSkip(num),
            'S' => Cigar::SoftClip(num),
            'H' => Cigar::HardClip(num),
            'P' => Cigar::Pad(num),
            '=' => Cigar::Equal(num),
            'X' => Cigar::Diff(num),
            _ => return Err(anyhow!("Unsupported CIGAR op {ch} in {cigar}")),
        };
        ops.push(op);
        num = 0;
    }
    Ok(CigarString(ops))
}

fn cigar_to_cigarstr(cigar: &Vec<(u32, u8)>) -> CigarString {
    let op_vec: Vec<Cigar> = cigar
        .to_owned()
        .iter()
        .map(|(len, op)| match op {
            0 => Cigar::Match(*len),
            1 => Cigar::Ins(*len),
            2 => Cigar::Del(*len),
            3 => Cigar::RefSkip(*len),
            4 => Cigar::SoftClip(*len),
            5 => Cigar::HardClip(*len),
            6 => Cigar::Pad(*len),
            7 => Cigar::Equal(*len),
            8 => Cigar::Diff(*len),
            _ => panic!("Unexpected cigar operation"),
        })
        .collect();
    CigarString(op_vec)
}

fn revcomp_in_place(seq: &mut [u8]) {
    seq.reverse();
    for b in seq.iter_mut() {
        *b = match *b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            b'N' | b'n' => b'N',
            _ => b'N',
        };
    }
}

fn reverse_in_place<T>(v: &mut [T]) {
    v.reverse();
}

fn reverse_qual(qual: &[u8]) -> Vec<u8> {
    let mut q = qual.to_vec();
    q.reverse();
    q
}

fn count_base(seq: &[u8], base: u8) -> usize {
    let b = base.to_ascii_uppercase();
    seq.iter().filter(|&&x| x.to_ascii_uppercase() == b).count()
}

fn parse_mm_part(part: &str) -> AnyResult<(u8, String, Vec<i64>)> {
    let mut it = part.split(',');
    let header = it.next().ok_or_else(|| anyhow!("Bad MM part: {part}"))?.to_string();
    let base = header
        .as_bytes()
        .get(0)
        .copied()
        .ok_or_else(|| anyhow!("Bad MM header: {header}"))?;

    let mut deltas: Vec<i64> = Vec::new();
    for x in it {
        if x.is_empty() {
            continue;
        }
        deltas.push(x.parse::<i64>()?);
    }
    Ok((base, header, deltas))
}

fn mm_deltas_to_abs(deltas: &[i64]) -> Vec<i64> {
    let mut abs = Vec::with_capacity(deltas.len());
    let mut cur: i64 = 0;
    for &d in deltas {
        cur += d;      // skip d bases of that type
        abs.push(cur); // current base index is modified
        cur += 1;      // move past this modified base
    }
    abs
}

fn mm_abs_to_deltas(abs: &[i64]) -> Vec<i64> {
    let mut deltas = Vec::with_capacity(abs.len());
    let mut last: i64 = 0;
    for &idx in abs {
        deltas.push(idx - last);
        last = idx + 1;
    }
    deltas
}

// fn reverse_mm_part(part: &str, seq_fwd_oriented: &[u8]) -> AnyResult<String> {
//     let (base, header, deltas) = parse_mm_part(part)?;
//     let total = count_base(seq_fwd_oriented, base) as i64;

//     if total <= 0 || deltas.is_empty() {
//         return Ok(format!("{header}"));
//     }

//     let mut abs = mm_deltas_to_abs(&deltas);

//     for x in abs.iter_mut() {
//         *x = (total - 1) - *x;
//     }
//     abs.sort_unstable();

//     let new_deltas = mm_abs_to_deltas(&abs);
//     let mut out = String::new();
//     out.push_str(&header);
//     for d in new_deltas {
//         out.push(',');
//         out.push_str(&d.to_string());
//     }
//     Ok(out)
// }

fn reverse_mm_part(part: &str, seq_fwd_oriented: &[u8]) -> AnyResult<String> {
    let (base, header, deltas) = parse_mm_part(part)?;
    let total = count_base(seq_fwd_oriented, base) as i64;

    let comp_base = match base {
        b'A' => b'T', b'C' => b'G', b'G' => b'C', b'T' => b'A',
        b'a' => b't', b'c' => b'g', b'g' => b'c', b't' => b'a',
        other => other,
    };

    let mut new_header = header.clone();
    if !new_header.is_empty() {
        unsafe { new_header.as_bytes_mut()[0] = comp_base; }
    }

    if total <= 0 || deltas.is_empty() {
        return Ok(new_header);
    }

    let mut abs = mm_deltas_to_abs(&deltas);

    for x in abs.iter_mut() {
        *x = (total - 1) - *x;
    }
    abs.sort_unstable();

    let new_deltas = mm_abs_to_deltas(&abs);
    let mut out = String::new();
    out.push_str(&new_header);
    for d in new_deltas {
        out.push(',');
        out.push_str(&d.to_string());
    }
    Ok(out)
}


fn reverse_mm_tag(mm: &str, seq_fwd_oriented: &[u8]) -> AnyResult<String> {
    let mm = mm.trim_end_matches('\0'); 
    let parts: Vec<&str> = mm.split(';').filter(|s| !s.is_empty()).collect();
    if parts.is_empty() {
        return Ok(mm.to_string());
    }
    let mut new_parts = Vec::with_capacity(parts.len());
    for p in parts {
        new_parts.push(reverse_mm_part(p, seq_fwd_oriented)?);
    }
    new_parts.reverse();
    Ok(format!("{};", new_parts.join(";")))
}


fn reverse_ml(ml: &[u8]) -> Vec<u8> {
    let mut v = ml.to_vec();
    v.reverse();
    v
}


fn copy_all_aux_except(rec_src: &Record, rec_dst: &mut Record, skip: &[&[u8; 2]]) -> AnyResult<()> {
    for item in rec_src.aux_iter() {
        let (tag, aux) = item?;
        if skip.iter().any(|t| t.as_slice() == tag) {
            continue;
        }
        rec_dst.push_aux(tag, aux)?;
    }
    Ok(())
}


#[derive(Debug, PartialEq, Eq)]
pub struct SeqMetaData {
    pub name: String,
    pub length: u32,
    pub is_alt: bool,
}

pub struct MMIndex {
    pub inner: Arc<MmIdx>,
}

impl MMIndex {
    pub fn n_seq(&self) -> u32 {
        unsafe { (**self.inner).n_seq }
    }

    pub fn seqs(&self) -> Vec<SeqMetaData> {
        let mut seqs: Vec<SeqMetaData> = Vec::with_capacity(self.n_seq() as usize);
        for i in 0..self.n_seq() {
            let _seq = unsafe { *(**self.inner).seq.offset(i as isize) };
            let c_str = unsafe { CStr::from_ptr(_seq.name) };
            let rust_str = c_str.to_str().unwrap().to_string();
            seqs.push(SeqMetaData {
                name: rust_str,
                length: _seq.len,
                is_alt: _seq.is_alt != 0,
            });
        }
        seqs
    }

    pub fn get_header(&self) -> Header {
        let mut header = Header::new();
        for seq in self.seqs() {
            header.push_record(
                HeaderRecord::new(b"SQ")
                    .push_tag(b"SN", &seq.name)
                    .push_tag(b"LN", &seq.length),
            );
        }
        header.push_record(
            HeaderRecord::new(b"PG")
                .push_tag(b"ID", "bamaligner")
                .push_tag(b"PN", "bamaligner"),
        );
        header
    }
}

impl From<&Aligner<Built>> for MMIndex {
    fn from(aligner: &Aligner<Built>) -> Self {
        MMIndex {
            inner: std::sync::Arc::clone(aligner.idx.as_ref().unwrap()),
        }
    }
}


fn process_single_read(
    original_record: &Record,
    aligner: &Aligner<Built>,
    tid_map: &HashMap<String, i32>,
) -> AnyResult<Vec<Record>> {
    let secondary = aligner.mapopt.flag & MM_F_NO_PRINT_2ND as i64 == 0;
    let soft_clip = aligner.mapopt.flag & MM_F_SOFTCLIP as i64 != 0;
    let qname = original_record.qname();
    let seq_buf = original_record.seq().as_bytes();
    let qual_buf = original_record.qual();

    let mm_opt: Option<String> = original_record.aux(b"MM").ok().and_then(|a| match a {
        Aux::String(s) => Some(s.to_string()),
        _ => None,
    });

    let ml_opt: Option<Vec<u8>> = original_record.aux(b"ML").ok().and_then(|a| match a {
        Aux::ArrayU8(arr) => Some(arr.iter().collect::<Vec<u8>>()),
        _ => None,
    });

    let alns = aligner
        .map(&seq_buf, false, false, None, None, None)
        .expect("Failed to map sequence");

    if alns.is_empty() {
        return Ok(Vec::new());
    }

    let mut output_records = Vec::with_capacity(alns.len());

    for (_i, aln) in alns.iter().enumerate() {
        let is_rev = aln.strand == minimap2::Strand::Reverse;
        let is_primary = aln.is_primary;
        let is_supplementary = aln.is_supplementary;
        
        if !secondary && !is_primary && !is_supplementary { continue; }
        
        let aln_info = aln
            .alignment
            .as_ref()
            .ok_or_else(|| anyhow!("Missing alignment info"))?;
        
        let query_start = aln.query_start;
        let query_end = aln.query_end;
        let raw_cigar = parse_cigar_string(aln_info.cigar_str.as_deref().unwrap_or("*"))?;
          

        
        let use_soft_clip = soft_clip || is_primary;
  
        let (final_seq, final_qual, final_cigar) = if use_soft_clip {
            let mut s = seq_buf.to_vec();
            let mut q = qual_buf.to_vec();
            if is_rev {
                revcomp_in_place(&mut s);
                q.reverse();
            }

            let ops: Vec<Cigar> = raw_cigar.iter().map(|op| match op {
                Cigar::HardClip(n) => Cigar::SoftClip(*n),
                _ => op.clone(),
            }).collect();
            (s, q, CigarString(ops))
        } else {
            let q_start = aln.query_start as usize;
            let q_end = aln.query_end as usize;
            let mut s = seq_buf[q_start..q_end].to_vec();
            let mut q = qual_buf[q_start..q_end].to_vec();
            if is_rev {
                revcomp_in_place(&mut s);
                q.reverse();
            }
            let ops: Vec<Cigar> = raw_cigar.iter()
                .filter(|c| !matches!(c, Cigar::SoftClip(_) | Cigar::HardClip(_)))
                .cloned().collect();
            (s, q, CigarString(ops))
        };
        
        let tid: i32 = if let Some(ref tname_arc) = aln.target_name {
            let tname_str: &str = tname_arc.as_str();
            *tid_map
                .get(tname_str)
                .ok_or_else(|| anyhow!("Unknown reference name in alignment: {tname_str}"))?
        } else {
            aln.target_id
        };

        let mut rec = Record::new();
        rec.set(qname, Some(&final_cigar), &final_seq, &final_qual);
        rec.set_tid(tid);
        rec.set_pos(aln.target_start as i64);
        rec.set_mtid(-1);
        rec.set_mpos(-1);
        rec.set_insert_size(0);
        rec.set_mapq(aln.mapq.try_into().unwrap_or(0));

        let mut flag: u16 = 0;
        if is_primary {
            flag |= 0x0;
        } else if is_supplementary {
            flag |= 0x800;
        } else {
            flag |= 0x100;
        }
        if is_rev {
            flag |= 0x10;
        }
        rec.set_flags(flag);

        copy_all_aux_except(&original_record, &mut rec, &[b"MM", b"ML"])?;

        rec.push_aux(b"NM", Aux::I32(aln_info.nm as i32))?;
        if let Some(ascore) = aln_info.alignment_score {
            rec.push_aux(b"AS", Aux::I32(ascore as i32))?;
        }

        let tp = if is_primary { b'P' } else { b'S' };
        rec.push_aux(b"tp", Aux::Char(tp))?;

        if is_primary || soft_clip {

            if let Some(mm) = mm_opt.as_ref() {
                rec.push_aux(b"MM", Aux::String(mm))?;
            }
            if let Some(ml) = ml_opt.as_ref() {
                rec.push_aux(b"ML", Aux::ArrayU8(ml.into()))?;
             
            }
        }
        output_records.push(rec);
    }
    Ok(output_records)
}

pub fn align(reference: &String, 
            input_bams: &Vec<String>, 
            output_bam: &String,
            k: Option<i16>, w: Option<i16>, 
            batch_size: u64,
            soft_clip: bool, 
            secondary: bool,
            rank_frac: Option<f32>,
            bw: Option<(i32, Option<i32>)>,
            min_cnt: Option<i32>,
            min_chain_score: Option<i32>,
            matching_score: Option<i32>, 
            mismatch_penalty: Option<i32>,
            gap_open: Option<(i32, Option<i32>)>, 
            gap_extension: Option<(i32, Option<i32>)>,
            z_drop: Option<(i32, Option<i32>)>,
            min_dp_max: Option<i32>,
            best_n: Option<i32>, 
            pri_ratio: Option<f32>, 
            mini_batch_size: i64,
            preset: &str,
            threads: usize,
        ) -> AnyResult<()>  {

    let preset = match preset {
        "lr:hq" => minimap2::Preset::LrHq,
        "map-hifi" => minimap2::Preset::MapHifi,
        "map-ont" => minimap2::Preset::MapOnt,
        "map-pb" => minimap2::Preset::MapPb,
        "sr" => minimap2::Preset::Sr,
        "splice" => minimap2::Preset::Splice,
        "splice:hq" => minimap2::Preset::SpliceHq,
        "asm5" => minimap2::Preset::Asm5,
        "asm10" => minimap2::Preset::Asm10,
        "asm20" => minimap2::Preset::Asm20,
        "ava-pb" => minimap2::Preset::AvaPb,
        "ava-ont" => minimap2::Preset::AvaOnt,
        _ => {
            log::warn!("Unknown preset: {}. Using default 'map-ont'", preset);
            minimap2::Preset::MapOnt
        }
    };
    
    let mut aligner = Aligner::builder()
                        .preset(preset)
                        .with_index_threads(threads)
                        .with_cigar()
                        .with_index(reference, None)
                        .unwrap();
    // set parameters
    if let Some(k) = k {
        aligner.idxopt.k = k;
    }
    if let Some(w) = w {
        aligner.idxopt.w = w;
    }

    aligner.idxopt.batch_size = batch_size;
    aligner.idxopt.mini_batch_size = mini_batch_size;

    if let Some(rf) = rank_frac {
        aligner.mapopt.rank_frac = rf;
    }

    if let Some((bw_val, bw2_opt)) = bw {
        aligner.mapopt.bw = bw_val;
        
        aligner.mapopt.bw_long = bw2_opt.unwrap_or(bw_val);
    }

    if let Some(n) = min_cnt {
        aligner.mapopt.min_cnt = n;
    }

    if let Some(m) = min_chain_score {
        aligner.mapopt.min_chain_score = m;
    }

    if let Some(a) = matching_score {
        aligner.mapopt.a = a;
    }

    if let Some(b) = mismatch_penalty {
        aligner.mapopt.b = b;
    }

    if let Some((e_val, e2_opt)) = gap_extension {
        aligner.mapopt.e = e_val;
        
        aligner.mapopt.e2 = e2_opt.unwrap_or(e_val);
    }

      
    if let Some((e_val, e2_opt)) = gap_extension {
        aligner.mapopt.e = e_val;
       
        aligner.mapopt.e2 = e2_opt.unwrap_or(e_val);
    }

    if let Some((z_drop, z_drop_inv)) = z_drop {
        aligner.mapopt.zdrop = z_drop;
        aligner.mapopt.zdrop_inv = z_drop_inv.unwrap_or(z_drop);
    }

    if let Some(min_dp_max) = min_dp_max {
        aligner.mapopt.min_dp_max = min_dp_max;
    }
    if secondary {
        aligner.mapopt.flag |= MM_F_SECONDARY_SEQ as i64;
    } else {
        aligner.mapopt.flag |= MM_F_NO_PRINT_2ND as i64;
    }

    if soft_clip {
        aligner.mapopt.flag |= MM_F_SOFTCLIP as i64;
    }
    
    if let Some(best_n) = best_n {
        aligner.mapopt.best_n = best_n;
    }

    if let Some(pri_ratio) = pri_ratio {
        aligner.mapopt.pri_ratio = pri_ratio;
    }
  
    let idx = aligner.idx.clone().unwrap();


    let mut header = Header::new();
    let mmindex = MMIndex::from(&aligner);
    let header = mmindex.get_header();
    let mut tid_map: HashMap<String, i32> = HashMap::new();
    for (i, seq) in mmindex.seqs().iter().enumerate() {
        tid_map.insert(seq.name.clone(), i as i32);
    }


    let mut writer = Writer::from_path(output_bam, &header, bam::Format::Bam)?;
    let _ = writer.set_threads((threads / 4).max(1));

    let (tx, rx) = channel::bounded::<AnyResult<Vec<Record>>>(threads * 2);

    thread::scope(|s| {
        s.spawn(move || {
            
            // bam.records().par_bridge().for_each_with(tx, |tx, r| {
            //     if let Ok(rec) = r {
            //         let res = process_single_read(&rec, &aligner, &tid_map);
            //         let _ = tx.send(res);
            //     }
            // });
            for bam_path in input_bams {

                log::info!("Starting alignment for: {}", bam_path);
                let mut reader = if bam_path == "-" {
                    Reader::from_stdin().expect("Failed to read from stdin")
                } else {
                    Reader::from_path(bam_path).expect("Failed to read BAM file")
                };
                let _ = reader.set_threads((threads / 4).max(1));

     
                reader.records().par_bridge().for_each_with(tx.clone(), |tx_inner, r| {
                    if let Ok(rec) = r {
                        let res = process_single_read(&rec, &aligner, &tid_map);
                        let _ = tx_inner.send(res);
                    }
                });
            }
        });

    
        let mut count = 0;
        for result in rx {
            let records = result?;
            for rec in records {
                writer.write(&rec)?;
            }
            count += 1;
            if count % mini_batch_size == 0 {
                log::info!("Processed {} reads...", count);
            }
        }
        log::info!("Finished. Total processed records: {}", count);

        Ok::<(), anyhow::Error>(())
    })?;


    // let (tx, rx) = mpsc::sync_channel(threads * 2);
    // thread::spawn(move || {
    //     bam.records().par_bridge().for_each(|r| {
    //         if let Ok(rec) = r {
    //             let res = process_single_read(&rec, &aligner, &tid_map);
    //             tx.send(res).ok();
    //         }
    //     });
    // });

    // let mut count = 0;
    // while let Ok(res) = rx.recv() {
    //     for out_rec in res? {
    //         writer.write(&out_rec)?;
    //     }
    //     count += 1;
    //     if count % mini_batch_size == 0 { log::info!("Mapped {} reads", count); }
    // }

    // let batch_size = mini_batch_size as usize;
    // let mut batch: Vec<Record> = Vec::with_capacity(batch_size);

    // for r in bam.records() {
    //     let original_record = r?;
    //     batch.push(original_record);

    //     if batch.len() >= batch_size {
    //         let results: Vec<AnyResult<Vec<Record>>> = batch
    //             .par_iter()
    //             .map(|rec| process_single_read(rec, &aligner, &tid_map))
    //             .collect();

    //         for res in results {
    //             let output_records = res?; 
    //             for out_rec in output_records {
    //                 writer.write(&out_rec)?;
    //             }
    //         }
    //         log::info!("mapped {} sequences", batch_size);
    //         batch.clear();
    //     }
    // }

    // if !batch.is_empty() {
    //     let results: Vec<AnyResult<Vec<Record>>> = batch
    //         .par_iter()
    //         .map(|rec| process_single_read(rec, &aligner, &tid_map))
    //         .collect();

    //     for res in results {
    //         let output_records = res?;
    //         for out_rec in output_records {
    //             writer.write(&out_rec)?;
    //         }
    //     }
        
    //     log::info!("mapped {} sequences", batch.len());
    // }

    Ok(())

}