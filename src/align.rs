#![allow(unused_imports, unused_variables, unused_unsafe)]
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
    Writer
};
use std::collections::HashMap;
use std::sync::Arc;
use std::sync::mpsc;
use std::thread;
use std::ffi::{CStr, CString};
use std::io::{BufRead, BufReader, Read as stdRead};
use std::path::Path;
use rayon::prelude::*;
use needletail::{ parse_fastx_file, parse_fastx_stdin, parser::SequenceRecord};



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
        let is_secondary_aln = !is_primary && !is_supplementary;
        if !secondary && is_secondary_aln { continue; }
        
        let aln_info = aln
            .alignment
            .as_ref()
            .ok_or_else(|| anyhow!("Missing alignment info"))?;
        
        let raw_cigar = parse_cigar_string(aln_info.cigar_str.as_deref().unwrap_or("*"))?;
        
        let (final_seq, final_qual, final_cigar) = if soft_clip {
            let mut s = seq_buf.to_vec();
            let mut q = qual_buf.to_vec();
            if is_rev {
                revcomp_in_place(&mut s);
                q.reverse();
            }
            (s, q, raw_cigar)
        } else if is_secondary_aln {
            (vec![], vec![], raw_cigar)
        } else if is_primary {
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
        } else if is_supplementary {
            let q_start = aln.query_start as usize;
            let q_end = aln.query_end as usize;
            let mut s = seq_buf[q_start..q_end].to_vec();
            let mut q = qual_buf[q_start..q_end].to_vec();
            if is_rev {
                revcomp_in_place(&mut s);
                q.reverse();
            }
            let ops: Vec<Cigar> = raw_cigar.iter().map(|op| match op {
                Cigar::SoftClip(n) => Cigar::HardClip(*n),
                _ => op.clone(),
            }).collect();
            (s, q, CigarString(ops))
        } else {
            continue;
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
            mask_level: Option<f32>,
            max_gap: Option<i32>,
            max_gap_ref: Option<i32>,
            max_frag_len: Option<i32>,
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
            max_qlen: Option<i32>,
            mini_batch_size: i64,
            seed: Option<i32>,
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
    
    log::info!("Indexing reference: {} ...", reference);
    let mut aligner = Aligner::builder()
                        .preset(preset);
    // set parameters
    if let Some(k) = k {
        aligner.idxopt.k = k;
    }
    if let Some(w) = w {
        aligner.idxopt.w = w;
    }

    aligner.idxopt.batch_size = batch_size;
    aligner.idxopt.mini_batch_size = mini_batch_size;

    if let Some(f) = mask_level {
        aligner.mapopt.mask_level = f;
    }

    if let Some(max_gap) = max_gap {
        aligner.mapopt.max_gap = max_gap;
    }

    if let Some(max_gap_ref) = max_gap_ref {
        aligner.mapopt.max_gap_ref = max_gap_ref;
    }
    if let Some(max_frag_len) = max_frag_len {
        aligner.mapopt.max_frag_len = max_frag_len;
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
        // aligner.mapopt.set_secondary_seq();
        // aligner.mapopt.unset_no_print_2nd();
    } else {
        aligner.mapopt.flag |= MM_F_NO_PRINT_2ND as i64;
        // aligner.mapopt.unset_secondary_seq();
        // aligner.mapopt.set_no_print_2nd();
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

    if let Some(max_qlen) = max_qlen {
        aligner.mapopt.max_qlen = max_qlen;
    }
    
    if let Some(seed) = seed {
        aligner.mapopt.seed = seed;
    }
    
    let aligner = aligner.with_cigar()
        .with_index_threads(threads)
        .with_index(reference, None)
        .unwrap();

    let idx = aligner.idx.clone().unwrap();

    log::info!("Finished indexing. Number of sequences in reference: {}", unsafe { (**idx).n_seq });
    let header = Header::new();
    let mmindex = MMIndex::from(&aligner);
    let header = mmindex.get_header();

    let tid_map: HashMap<String, i32> = mmindex
        .seqs()
        .into_par_iter()
        .enumerate()
        .map(|(i, seq)| (seq.name, i as i32))
        .collect();
    if input_bams.len() > 1 {
        log::info!("Starting alignment of input files...");
    }


    let mut writer = Writer::from_path(output_bam, &header, bam::Format::Bam)?;
    let _ = writer.set_threads(threads.min(32));

    let cap = (threads * 2).clamp(32, 256);
    let (tx, rx) = channel::bounded::<AnyResult<Vec<Record>>>(cap);

    let aligner_ref = &aligner;
    let tid_map_ref = &tid_map;

    thread::scope(|s| {
        let tx_producer = tx.clone();
        s.spawn(move || {
            for path in input_bams {

                log::info!("Starting alignment for: {}", path);

                let is_stdin = path == "-" || path == "/dev/stdin";
                let is_pipe = path.starts_with("/dev/fd/");

                let is_bam = if is_stdin || is_pipe {
                    true
                } else {
                    let p = Path::new(path);
                    match p.extension() {
                        Some(ext) => {
                            let ext_str = ext.to_string_lossy().to_lowercase();
                            ext_str == "bam" || ext_str == "cram" || ext_str == "sam"
                        }
                        None => false,
                    }
                };
                
                if is_bam {
                    log::info!("Processing HTS file for: {}", path);
                    let mut reader = if is_stdin {
                        Reader::from_stdin().expect("Failed to read from stdin")
                    } else {
                        Reader::from_path(path).expect("Failed to read BAM file")
                    };
                    let _ = reader.set_threads(threads.min(16));
                    reader.records().par_bridge().for_each_with(tx_producer.clone(), |tx_inner, r| {
                        if let Ok(rec) = r {
                            if let Ok(res) = process_single_read(&rec, aligner_ref, tid_map_ref) {
                                if !res.is_empty() {
                                    let _ = tx_inner.send(Ok(res));
                                }
                            }
                        }
                    });
                    continue; 
                }

                if let Ok(mut n_reader) = parse_fastx_file(path) {
                    log::info!("Detected Fastx format for: {}", path);
                    let iter = std::iter::from_fn(move || {
                        n_reader.next().map(|r| {
                            let rec = r.expect("Invalid fastx record");
                            let mut bam_rec = Record::new();
                            let qual = rec.qual();
                            if let Some(q) = qual.as_deref() {
                                bam_rec.set(rec.id(), None, &rec.seq(), q);
                            } else {
                                let dummy_qual = vec![255u8; rec.seq().len()];
                                bam_rec.set(rec.id(), None, &rec.seq(), &dummy_qual);
                            }
                            bam_rec
                        })
                    });
                    iter.par_bridge().for_each_with(tx_producer.clone(), |tx_inner, bam_rec| {
                        if let Ok(res) = process_single_read(&bam_rec, aligner_ref, tid_map_ref) {
                            if !res.is_empty() {
                                let _ = tx_inner.send(Ok(res));
                            }
                        }
                    });
                    continue; 
                } else {
                    log::info!("Processing HTS file for: {}", path);
                    let mut reader = if is_stdin {
                        Reader::from_stdin().expect("Failed to read from stdin")
                    } else {
                        Reader::from_path(path).expect("Failed to read BAM file")
                    };
                    let _ = reader.set_threads(threads.min(16));
                    reader.records().par_bridge().for_each_with(tx_producer.clone(), |tx_inner, r| {
                        if let Ok(rec) = r {
                          
                            if let Ok(res) = process_single_read(&rec, aligner_ref, tid_map_ref) {
                                if !res.is_empty() {
                                    let _ = tx_inner.send(Ok(res));
                                }
                            }
                        }
                    });
                }
                
            }
        }); 

        drop(tx);

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

    Ok(())

}