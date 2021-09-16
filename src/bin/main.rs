use clap::{App, Arg};
use ecm::ecm_multithreaded;
use ecm::ecm_singlethreaded;
use log::info;
use rug::Integer;
use simple_logger;
use std::str::FromStr;
use std::sync::Arc;
extern crate hwloc;
use hwloc::{ObjectType, Topology};

fn main() {
    let matches = App::new("ECM Factorization")
        .version("1.0")
        .author("")
        .about("Factors integers using the Elliptic Curve Method")
        .arg(
            Arg::new("number")
                .about("Number to factor using ECM")
                .takes_value(true)
                .short('n')
                .long("number")
                .required(true),
        )
        .arg(
            Arg::new("num_curves")
                .about("Number of curves to try out")
                .takes_value(true)
                .short('c')
                .long("num_curves")
                .required(false),
        )
        .arg(
            Arg::new("verbose")
                .about("Detailed execution")
                .takes_value(false)
                .short('v')
                .long("verbose")
                .required(false),
        )
        .arg(
            Arg::new("debug")
                .about("Debug information")
                .takes_value(false)
                .short('d')
                .long("debug")
                .required(false),
        )
        .arg(
            Arg::new("b1_bound")
                .about("Stage 1 bound")
                .takes_value(true)
                .long("b1")
                .required(false),
        )
        .arg(
            Arg::new("b2_bound")
                .about("Stage 2 bound")
                .takes_value(true)
                .long("b2")
                .required(false),
        )
        .arg(
            Arg::new("sigma")
                .about("Curve's parameter")
                .takes_value(true)
                .short('s')
                .long("sigma")
                .required(false),
        )
        .arg(
            Arg::new("single_threaded")
                .about("Run on a single thread\nNote: the program is multi-threaded by default, using as many threads as there are cores available")
                .takes_value(false)
                .long("single_threaded")
                .required(false),
        )
        .get_matches();

    if matches.is_present("debug") {
        simple_logger::SimpleLogger::new()
            .with_level(log::LevelFilter::Debug)
            .init()
            .unwrap();
    } else if matches.is_present("verbose") {
        simple_logger::SimpleLogger::new()
            .with_level(log::LevelFilter::Info)
            .init()
            .unwrap();
    };

    if let Some(n) = matches.value_of("number") {
        match Integer::from_str(n) {
            Err(_) => println!("Wrong input"),
            Ok(n) => {
                let b1: u64 = match matches.value_of("b1_bound") {
                    Some(s) => s.parse::<u64>().unwrap(),
                    None => 10000,
                };
                let b2: u64 = match matches.value_of("b2_bound") {
                    Some(s) => s.parse::<u64>().unwrap(),
                    None => 100 * b1,
                };
                let curves = Arc::new(
                    matches
                        .value_of("num_curves")
                        .and_then(|s| Integer::from_str(s).ok()),
                );
                let sigma = matches
                    .value_of("sigma")
                    .and_then(|s| Integer::from_str(s).ok());
                if matches.is_present("single_threaded") || !sigma.is_none() {
                    match ecm_singlethreaded(&n, &curves, b1, b2, &Arc::new(sigma)) {
                        Some(f) => print!("Found factor {}.\n", f),
                        None => print!("No factor found.\n"),
                    }
                } else {
                    let topology = Topology::new();

                    // Get all objects with type "Core"
                    let cores = topology.objects_with_type(&ObjectType::Core);
                    let nthreads = match cores {
                        Ok(c) => c.len(),
                        Err(_) => 1, // fallback to one thread
                    };
                    info!("Found {} cores, spawning {} threads", nthreads, nthreads);

                    match ecm_multithreaded(&n, &curves, b1, b2, &Arc::new(sigma), nthreads) {
                        Some(f) => print!("Found factor {}.\n", f),
                        None => print!("No factor found.\n"),
                    }
                }
            }
        }
    }
}
