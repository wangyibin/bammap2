
use clap::{arg, Arg, ArgAction, 
            builder::{
                styling::{AnsiColor, Effects},
                Styles,
            },
            Command, 
            value_parser};

const VERSION: &str = "0.1.4";

const STYLES: Styles = Styles::styled()
    .header(AnsiColor::Green.on_default().effects(Effects::BOLD))
    .usage(AnsiColor::Green.on_default().effects(Effects::BOLD))
    .literal(AnsiColor::Cyan.on_default().effects(Effects::BOLD))
    .placeholder(AnsiColor::Yellow.on_default());

fn parse_si_size_u64(s: &str) -> Result<u64, String> {
    let s = s.trim();
    if s.is_empty() {
        return Err("Empty input".to_string());
    }

    let last = s.chars().last().unwrap();
    let (num_str, multiplier) = match last {
        'G' | 'g' => (&s[..s.len() - 1], 1_000_000_000.0),
        'M' | 'm' => (&s[..s.len() - 1], 1_000_000.0),
        'K' | 'k' => (&s[..s.len() - 1], 1_000.0),
        _ => (s, 1.0),
    };

    let num = num_str.parse::<f64>()
        .map_err(|_| format!("Invalid number format: {}", s))?;

    Ok((num * multiplier) as u64)
}

fn parse_si_size_i64(s: &str) -> Result<i64, String> {
    let s = s.trim();
    if s.is_empty() {
        return Err("Empty input".to_string());
    }

    let last = s.chars().last().unwrap();
    let (num_str, multiplier) = match last {
        'G' | 'g' => (&s[..s.len() - 1], 1_000_000_000.0),
        'M' | 'm' => (&s[..s.len() - 1], 1_000_000.0),
        'K' | 'k' => (&s[..s.len() - 1], 1_000.0),
        _ => (s, 1.0),
    };

    let num = num_str.parse::<f64>()
        .map_err(|_| format!("Invalid number format: {}", s))?;

    Ok((num * multiplier) as i64)
}

fn parse_si_size_i32(s: &str) -> Result<i32, String> {
    let s = s.trim();
    if s.is_empty() {
        return Err("Empty input".to_string());
    }

    let last = s.chars().last().unwrap();
    let (num_str, multiplier) = match last {
        'G' | 'g' => (&s[..s.len() - 1], 1_000_000_000.0),
        'M' | 'm' => (&s[..s.len() - 1], 1_000_000.0),
        'K' | 'k' => (&s[..s.len() - 1], 1_000.0),
        _ => (s, 1.0),
    };

    let num = num_str.parse::<f64>()
        .map_err(|_| format!("Invalid number format: {}", s))?;

    Ok((num * multiplier) as i32)
}


fn parse_int_pair(s: &str) -> Result<(i32, Option<i32>), String> {
    let s = s.trim();
    if s.is_empty() {
        return Err("Empty input".to_string());
    }

    if let Some((first, second)) = s.split_once(',') {
        let v1 = first.trim().parse::<i32>()
            .map_err(|_| format!("Invalid integer format: {}", first))?;
        let v2 = second.trim().parse::<i32>()
            .map_err(|_| format!("Invalid integer format: {}", second))?;
        Ok((v1, Some(v2)))
    } else {
        let v1 = s.parse::<i32>()
            .map_err(|_| format!("Invalid integer format: {}", s))?;
        Ok((v1, None))
    }
}


pub fn cli() -> Command {
    Command::new("bammap2")
        .version(VERSION)
        .styles(STYLES)
        .about("A lightweight wrapper to run minimap2 (v2.30) alignments directly on BAM file")
        .arg(
            arg!(<reference> "target.fa or target.idx" )
        )
        .arg(
            Arg::new("query")
                .help("query sequences with BAM or fastq(.gz)/fasta(.gz), multiple files allowed. Use '-' for stdin. Stdin and pipe are only supported for BAM input and single file.")
                .required(true)
                .num_args(1..)
        )
        .next_help_heading("Indexing")
        .arg(
            Arg::new("hpc")
                .short('H')
                .hide(true)
                .help("Use homopolymer-compressed k-mer (always true for map-pb/ont)")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("kmer")
                .short('k')
                .help("k-mer size (no larger than 28) [15]")
                .value_parser(value_parser!(i16))
                .value_name("INT")
        )
        .arg(
            Arg::new("window")
                .short('w')
                .help("minimizer window size [5]")
                .value_parser(value_parser!(i16))
                .value_name("INT")
        )
        .arg(
            Arg::new("batch_size")
            .short('I')
            .help("split index for every ~NUM input bases [16G]")
            .value_parser(parse_si_size_u64)
            .value_name("NUM")
            .default_value("16G"),
        
        )
        .next_help_heading("Mapping")
        .arg(
            Arg::new("mask_level")
                .short('f')
                .help("filter out top FLOAT fraction of repetitive minimizers [0.0002]")
                .value_parser(value_parser!(f32))
                .value_name("FLOAT")
            
        )
        .arg(
            Arg::new("max_gap")
                .short('g')
                .help("stop chain enlongation if there are no minimizers in INT-bp [5000]")
                .value_parser(parse_si_size_i32)
                .value_name("INT")
        )
        .arg(
            Arg::new("max_gap_ref")
                .short('G')
                .help("max intron length (effective with -xsplice; changing -r) [200k]")
                .value_parser(parse_si_size_i32)
                .value_name("INT")
        )
        .arg(
            Arg::new("max_frag_len")
                .short('F')
                .help("max fragment length (effective with -xsr or in the fragment mode) [800]")
                .value_parser(parse_si_size_i32)
                .value_name("INT")
        )
        .arg(
            Arg::new("bw")
                .short('r')
                .help("chaining/alignment bandwidth and long-join bandwidth [500,20000]")
                .value_parser(parse_int_pair)
                .value_name("INT,[INT]"),
        )
        .arg(
            Arg::new("min_cnt")
                .short('n')
                .help("minimal number of minimizers on a chain [3]")
                .value_parser(value_parser!(i32))
                .value_name("INT")
        )
        .arg(
            Arg::new("min_chain_score")
                .short('m')
                .help("minimal chaining score (matching bases minus log gap penalty) [40]")
                .value_parser(value_parser!(i32))
                .value_name("INT")
        )
        .arg(
            Arg::new("pri_ratio")
                .short('p')
                .help("Min secondary-to-primary score ratio [0.8]")
                .value_parser(value_parser!(f32))
                .value_name("FLOAT"),
        )
        .arg(
            Arg::new("best_n")
                .short('N')
                .help("Retain at most N secondary alignments [5]")
                .value_parser(value_parser!(i32))
                .value_name("INT"),
        )
        
        .next_help_heading("Alignments")
        .arg(
            Arg::new("matching_score")
                .short('A')
                .help("matching score [2]")
                .value_parser(value_parser!(i32))
                .value_name("INT"),
        )
        .arg(
            Arg::new("mismatch_penalty")
                .short('B')
                .help("mismatch penalty (larger value for lower divergence) [4]")
                .value_parser(value_parser!(i32))
                .value_name("INT"),
        )
        .arg(
            Arg::new("gap_open")
                .short('O')
                .help("gap open penalty [4,24]")
                .value_parser(parse_int_pair)
                .value_name("INT,[INT]"),
        )
        .arg(
            Arg::new("gap_extension")
                .short('E')
                .help("gap extension penalty; a k-long gap costs min{O1+k*E1,O2+k*E2} [2,1]")
                .value_parser(parse_int_pair)
                .value_name("INT,[INT]"),
        )
        .arg(
            Arg::new("z_drop")
                .short('z')
                .help("Z-drop score and inversion Z-drop score [400,200]")
                .value_parser(parse_int_pair)
                .value_name("INT,[INT]"),
        )
        .arg(
            Arg::new("min_dp_max")
                .short('s')
                .help("minimal peak DP alignment score [80]")
                .value_parser(value_parser!(i32))
                .value_name("INT")
        )

        .next_help_heading("Input/Output")
        
        .arg(
            Arg::new("output_bam")
                .short('a')
                .help("output in BAM format")
                .action(ArgAction::SetTrue)
                .hide(true)
        )
        .arg(
            Arg::new("output")
                .short('o')
                .help("output file path, currently only support output in BAM formst [stdout]")
                .value_parser(value_parser!(String))
                .default_value("-")
                .value_name("FILE"),
        )

        .arg(
            Arg::new("soft_clip")
                .short('Y')
                .help("use soft clipping for supplementary alignments")
                .action(ArgAction::SetTrue),
                
        )
        .arg(
            Arg::new("seed")
                .long("seed")
                .help("Integer seed for randomizing equally best hits. Minimap2 hashes INT and read name when choosing between equally best hits. [11]" )
                .value_parser(value_parser!(i32))
                .value_name("INT")
                .default_value("11")
        )
        .arg(
            Arg::new("max_qlen")
                .long("max-qlen")
                .help("skip reads longer than INT [0 (disabled)]")
                .value_parser(value_parser!(i64))
                .value_name("INT")
        )
        .arg(
            Arg::new("secondary")
                .long("secondary")
                .help("Whether to output secondary alignments" )
                .value_parser(["yes", "no"])
                .default_value("yes")
        )
        .arg(
            Arg::new("threads")
                .short('t')
                .help("number of threads")
                .value_parser(value_parser!(usize))
                .value_name("INT")
                .default_value("8"),
        )
        .arg(
            Arg::new("mini_batch_size")
                .short('K')
                .help("minibatch size logging when mapping")
                .value_parser(parse_si_size_i64)
                .value_name("STR")
                .default_value("10k"),

        )
        .next_help_heading("Presets")
        .arg(
            Arg::new("preset")
                .short('x')
                .value_name("STR")
                .help(
r##"- lr:hq - accurate long reads (error rate <1%) against a reference genome
- splice/splice:hq - spliced alignment for long reads/accurate long reads
- asm5/asm10/asm20 - asm-to-ref mapping, for ~0.1/1/5% sequence divergence
- sr - short reads against a reference
- map-pb/map-hifi/map-ont - CLR/HiFi/Nanopore vs reference mapping
- ava-pb/ava-ont - PacBio CLR/Nanopore read overlap
"##
            )
                .value_parser([
                    "lr:hq", "map-hifi", "map-pb", "map-ont", "asm5", "asm10", "asm20", "sr", "splice", "splice:hq"
                ])
                .default_value("lr:hq"),
        )
        .arg_required_else_help(true)

}

