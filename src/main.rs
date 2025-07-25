use clap::Arg;
use std::fs::File;
use std::io;
use std::io::BufReader;
use std::io::Write;
use std::io::{BufRead};
use std::collections::HashMap;
use std::io::BufWriter;

fn main() {
    let matches = clap::App::new("GTFTPMCombinator")
        .version("0.1")
        .author("Frankie B")
        .about("Converts Combining GTF & TPM mstrix to star end aligned format")
        .arg(
            clap::Arg::new("gtf")
                .short('g')
                .long("gtf")
                .takes_value(true)
                .about("Input GTF file")
                .required(true),
        )
        .arg(
            clap::Arg::new("tpm")
                .short('m')
                .long("tpm")
                .takes_value(true)
                .about("Input TPM Matrix(tsv)")
                .required(true),
        )
        .arg(
            clap::Arg::new("output")
                .short('o')
                .long("output")
                .takes_value(true)
                .about("Output file(tsv)")
                .required(false)
                .default_value("./output.tsv"),
        )
        .arg(
            Arg::new("threads")
                .short('t')
                .long("threads")
                .takes_value(true)
                .required(false)
                .about("Number of threads")
                .default_value("1"),
        )
        .get_matches();

    let threads = matches
        .value_of("threads")
        .unwrap()
        .parse::<usize>()
        .unwrap();
    let gtf_file = matches.value_of("gtf").unwrap();
    let tpm_matrix = matches.value_of("tpm").unwrap();
    let output_file = matches.value_of("output").unwrap();
    eprintln!("Starting");
    eprintln!("Threads: {}", threads);
    eprintln!("Graph: {}", gtf_file);
    eprintln!("Alignment: {}", tpm_matrix);
    eprintln!("Output: {}", output_file);
        run(gtf_file, tpm_matrix, output_file, threads);
}


pub fn run(gtf_file: &str, tpm_matrix: &str, output_file: &str, _threads: usize) {
    // Build a lookup table: gene_id -> (start, end)
    let gene_map = parse_gtf(gtf_file).expect("Failed to parse GTF");

    let input = File::open(tpm_matrix).expect("Error opening TPM matrix");
    let reader = BufReader::new(input);
    let output = File::create(output_file).expect("Error creating output file");
    let mut writer = BufWriter::new(output);

    // Handle header row separately so we can insert "start" and "end" column names
    let mut lines_iter = reader.lines();
    if let Some(Ok(header_line)) = lines_iter.next() {
        let mut header_cols: Vec<String> =
            header_line.split('\t').map(|s| s.to_string()).collect();

        // Only add the extra columns if this indeed looks like a header row
        if header_cols.get(0).map(|s| s.as_str()) == Some("gene_id") {
            header_cols.insert(1, "start".to_string());
            header_cols.insert(2, "end".to_string());
        }

        writeln!(writer, "{}", header_cols.join("\t")).expect("Error writing header");
    }

    // Process the remaining (data) rows
    for line in lines_iter {
        let line = line.expect("Error reading TPM matrix");
        if line.trim().is_empty() {
            continue;
        }

        let mut cols: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
        let gene_id = &cols[0];

        if let Some((start, end)) = gene_map.get(gene_id) {
            cols.insert(1, start.clone());
            cols.insert(2, end.clone());
        } else {
            cols.insert(1, String::from(""));
            cols.insert(2, String::from(""));
        }

        writeln!(writer, "{}", cols.join("\t")).expect("Error writing output");
    }
}

fn parse_gtf(gtf_file: &str) -> io::Result<HashMap<String, (String, String)>> {
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
            genes.insert(gene_id, (start, end));
        }
    }

    Ok(genes)
}
