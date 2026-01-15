#![allow(unused_imports, unused_variables)]
use env_logger::Builder;
use log::LevelFilter;
use jemallocator::Jemalloc;
use chrono::Local;
use bammap2::cli::cli;
use bammap2::align::align;
use std::io::Write;
use std::sync::Once;
static INIT: Once = Once::new();

#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;

pub fn setup_logging() {
    INIT.call_once(|| {
        env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info"))
            .is_test(true) 
            .try_init()
            .ok();
    });
}


fn main() {
    let start_time = std::time::Instant::now();

    setup_logging();

    let matches = cli().get_matches();
    let reference = matches.get_one::<String>("reference").unwrap();
    let queries: Vec<String> = matches
                                .get_many::<String>("query")
                                .unwrap()
                                .cloned()
                                .collect();
    let output = matches.get_one::<String>("output").unwrap();
    let is_hpc = matches.get_flag("hpc");
    let kmer = matches.get_one::<i16>("kmer").copied();
    let window = matches.get_one::<i16>("window").copied();
    let mask_level = matches.get_one::<f32>("mask_level").copied();
    let max_gap = matches.get_one::<i32>("max_gap").copied();
    let max_gap_ref = matches.get_one::<i32>("max_gap_ref").copied();
    let max_frag_len = matches.get_one::<i32>("max_frag_len").copied();
    let bw = matches.get_one::<(i32, Option<i32>)>("bw").copied();
    let min_cnt = matches.get_one::<i32>("min_cnt").copied();
    let min_chain_score = matches.get_one::<i32>("min_chain_score").copied();

    let matching_score = matches.get_one::<i32>("matching_score").copied();
    let mismatch_penalty = matches.get_one::<i32>("mismatch_penalty").copied();
    let gap_open = matches.get_one::<(i32, Option<i32>)>("gap_open").copied();
    let gap_extension = matches.get_one::<(i32, Option<i32>)>("gap_extension").copied();
    let z_drop = matches.get_one::<(i32, Option<i32>)>("z_drop").copied();
    let min_dp_max = matches.get_one::<i32>("min_dp_max").copied();
    
    let max_qlen = matches.get_one::<i32>("max_qlen").copied();
    let batch_size = matches.get_one::<u64>("batch_size").unwrap();
    let soft_clip = matches.get_flag("soft_clip");
    let secondary = matches.get_one::<String>("secondary").unwrap();
    let best_n = matches.get_one::<i32>("best_n").copied();
    let pri_ratio = matches.get_one::<f32>("pri_ratio").copied();
    let mini_batch_size = matches.get_one::<i64>("mini_batch_size").unwrap();
    let seed = matches.get_one::<i32>("seed").copied();
    let preset = matches.get_one::<String>("preset").unwrap();

    let threads = matches.get_one::<usize>("threads").unwrap();

    rayon::ThreadPoolBuilder::new().num_threads(*threads).build_global().unwrap();

    let secondary = match secondary.as_str() {
        "yes" => true,
        "no" => false,
        _ => false,
    };

    log::info!("Starting bammap2...");
   
    align(reference, &queries, output,
            kmer, window, *batch_size,
            soft_clip, secondary,
            mask_level, 
            max_gap, 
            max_gap_ref,
            max_frag_len,
            bw, min_cnt,
            min_chain_score,
            matching_score, mismatch_penalty,
            gap_open, gap_extension, z_drop,
            min_dp_max,
            best_n, pri_ratio,
            max_qlen,
            *mini_batch_size,
            seed,
            preset.as_str(),
            *threads
        ).unwrap();
    

    // output command line 
    let args: Vec<String> = std::env::args().collect();
    log::info!("CMD: {}", args.join(" "));

    let duration = start_time.elapsed();
    unsafe {
        let mut usage = std::mem::zeroed::<libc::rusage>();
        if libc::getrusage(libc::RUSAGE_SELF, &mut usage) == 0 {
            let cpu_time = (usage.ru_utime.tv_sec as f64 + usage.ru_utime.tv_usec as f64 / 1_000_000.0)
                + (usage.ru_stime.tv_sec as f64 + usage.ru_stime.tv_usec as f64 / 1_000_000.0);
            let peak_rss = usage.ru_maxrss as f64 / 1024.0 / 1024.0; 
            log::info!(
                "Real time: {:.3} sec; CPU: {:.3} sec; Peak RSS: {:.3} GB",
                duration.as_secs_f64(),
                cpu_time,
                peak_rss
            );
        }
    }
} 
