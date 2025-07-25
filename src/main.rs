use clap::{value_parser, Arg, Command};
use std::fs::File;
use std::io;
use std::io::BufReader;
use std::io::Write;
use std::io::{BufRead};
use std::collections::HashMap;
use std::io::BufWriter;

fn main() {
    let matches = Command::new("GTFTPMCombinator")
        .version("0.1")
        .author("Frankie B")
        .about("Converts Combining GTF & TPM mstrix to star end aligned format")
        .arg(
            clap::Arg::new("gtf")
                .short('g')
                .long("gtf")
                .num_args(1)
                .help("Input GTF file")
                .required(true),
        )
        .arg(
            clap::Arg::new("tpm")
                .short('m')
                .long("tpm")
                .num_args(1)
                .help("Input TPM Matrix(tsv)")
                .required(true),
        )
        .arg(
            clap::Arg::new("output")
                .short('o')
                .long("output")
                .num_args(1)
                .help("Output file(tsv)")
                .required(false)
                .default_value("./output.tsv"),
        )
        .arg(
            Arg::new("threads")
                .short('t')
                .long("threads")
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("Number of threads")
                .default_value("1"),
        )
        .arg(
            Arg::new("ignore")
                .long("ignore")
                .short('i')
                .num_args(1)
                .value_parser(value_parser!(u8).range(0..=4))
                .default_value("4")
                .help("Filter level (0-4): 0 keep all; 1 exclude missing; 2 also exclude non-'chr*'; 3 also exclude special chromosomes (chrX/Y/M); 4 only keep chromosomes 1-22"),
        )
        .get_matches();

    let threads: usize = *matches
        .get_one::<usize>("threads")
        .expect("threads has default");

    let gtf_file = matches
        .get_one::<String>("gtf")
        .expect("required")
        .as_str();

    let tpm_matrix = matches
        .get_one::<String>("tpm")
        .expect("required")
        .as_str();

    let output_file = matches
        .get_one::<String>("output")
        .expect("output has default")
        .as_str();

    let ignore_level = *matches.get_one::<u8>("ignore").unwrap();
    eprintln!("Starting");
    eprintln!("Threads: {}", threads);
    eprintln!("GTF File: {}", gtf_file);
    eprintln!("TPM Matrix: {}", tpm_matrix);
    eprintln!("Output: {}", output_file);
    eprintln!("Ignore level: {}", ignore_level);
    run(gtf_file, tpm_matrix, output_file, threads, ignore_level);
}

pub fn run(
    gtf_file: &str,
    tpm_matrix: &str,
    output_file: &str,
    _threads: usize,
    ignore_level: u8,
) {
    // Build a lookup table: gene_id -> (chr, start, end)
    let gene_map = parse_gtf(gtf_file).expect("Failed to parse GTF");

    let input = File::open(tpm_matrix).expect("Error opening TPM matrix");
    let reader = BufReader::new(input);
    let output = File::create(output_file).expect("Error creating output file");
    let mut writer = BufWriter::new(output);

    // Read header
    let mut lines_iter = reader.lines();
    let header_line = lines_iter.next().unwrap().unwrap();
    let mut header_cols: Vec<String> = header_line.split('\t').map(|s| s.to_string()).collect();

    // Collect rows for filtering and sorting
    let mut rows: Vec<(String, String, String, Vec<String>)> = Vec::new();
    for line in lines_iter {
        let line = line.expect("Error reading TPM matrix");
        if line.trim().is_empty() {
            continue;
        }
        let cols: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
        let gene_id = &cols[0];
        let (raw_chr, raw_start, raw_end) = gene_map
            .get(gene_id)
            .cloned()
            .unwrap_or_else(|| (String::new(), String::new(), String::new()));

        // Define categories
        let is_empty = raw_chr.is_empty() || raw_start.is_empty() || raw_end.is_empty();
        let is_independent = !raw_chr.starts_with("chr");
        let suffix = raw_chr.strip_prefix("chr").unwrap_or(raw_chr.as_str());
        let is_special = matches!(suffix, "X" | "Y" | "M");
        let num_val = suffix.parse::<u64>().unwrap_or(0);
        let is_numeric = suffix.parse::<u64>().is_ok();

        // Apply ignore_level filter
        match ignore_level {
            0 => {},
            1 if is_empty => continue,
            2 if is_empty || is_independent => continue,
            3 if is_empty || is_independent || is_special => continue,
            4 if is_empty || is_independent || is_special || !is_numeric || num_val > 22 => continue,
            _ => {},
        }

        // Determine display chromosome: strip "chr" prefix for level 4
        let chr_display = if ignore_level == 4 {
            raw_chr.strip_prefix("chr").unwrap_or(&raw_chr).to_string()
        } else {
            raw_chr.clone()
        };

        rows.push((chr_display, raw_start.clone(), raw_end.clone(), cols.clone()));
    }

    // Always sort by chromosome order then start
    rows.sort_by(|a, b| {
        let a_suffix = a.0.strip_prefix("chr").unwrap_or(&a.0);
        let b_suffix = b.0.strip_prefix("chr").unwrap_or(&b.0);
        let a_key = match a_suffix.parse::<u64>() {
            Ok(n) => n,
            Err(_) => match a_suffix {
                "X" => 23,
                "Y" => 24,
                "M" => 25,
                _ => u64::MAX,
            },
        };
        let b_key = match b_suffix.parse::<u64>() {
            Ok(n) => n,
            Err(_) => match b_suffix {
                "X" => 23,
                "Y" => 24,
                "M" => 25,
                _ => u64::MAX,
            },
        };
        if a_key != b_key {
            a_key.cmp(&b_key)
        } else {
            let a_start = a.1.parse::<u64>().unwrap_or(0);
            let b_start = b.1.parse::<u64>().unwrap_or(0);
            a_start.cmp(&b_start)
        }
    });

    // Write output: header then rows
    let mut full_header = vec![
        "Chr".to_string(),
        "start".to_string(),
        "end".to_string(),
    ];
    full_header.append(&mut header_cols);
    writeln!(writer, "#{}", full_header.join("\t")).expect("Error writing header");

    for (chr_val, start, end, mut cols) in rows {
        let mut out_cols = vec![chr_val, start, end];
        out_cols.append(&mut cols);
        writeln!(writer, "{}", out_cols.join("\t")).expect("Error writing output");
    }
}

fn parse_gtf(gtf_file: &str) -> io::Result<HashMap<String, (String, String, String)>> {
    let mut genes = HashMap::new();
    let file = File::open(gtf_file)?;
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 || fields[2] != "gene" {
            continue;
        }

        let start = fields[3].to_string();
        let end = fields[4].to_string();

        // Extract gene_id from the first attribute (gene_id "<id>";
        let gene_id = fields[8]
            .split(';')
            .next()
            .and_then(|attr| attr.split('"').nth(1))
            .unwrap_or("")
            .to_string();

        if !gene_id.is_empty() {
            let chr = fields[0].to_string();
            genes.insert(gene_id, (chr, start, end));
        }
    }

    Ok(genes)
}
