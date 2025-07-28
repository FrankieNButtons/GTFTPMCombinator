# GTFTPMCombinator

`GTFTPMCombinator` is a command line utility written in Rust. It merges a GTF gene annotation file with a TPM expression matrix to produce a chromosome sorted TSV table. The output is useful when downstream analysis requires gene position information alongside expression levels.

## Features

- Extracts `gene_id`, chromosome, start and end position from a GTF file.
- Reads a TPM expression matrix where the first column is `gene_id`.
- Prepends `Chr`, `start` and `end` columns before the TPM data.
- Applies filtering rules based on the chosen ignore level and sorts the final rows by chromosome and start coordinate.

## Building

This project uses [Cargo](https://www.rust-lang.org/). After cloning the repository, run the following in the project root:

```bash
cargo build --release
```

The executable will be placed in `target/release/GTFTPMCombinator`.

> If dependencies fail to download you may need to run `cargo fetch` or ensure access to crates.io.

## Usage

```bash
GTFTPMCombinator \
  --gtf <annotation.gtf> \
  --tpm <matrix.tsv> \
  [--output <out.tsv>] \
  [--threads <n>] \
  [--ignore <level>]
```

The program prints diagnostic information to stderr and writes the result to the selected output file.

### Arguments

- `-g, --gtf <file>`: **required** path to the input GTF file.
- `-m, --tpm <file>`: **required** TPM matrix in TSV format; the first column must be `gene_id`.
- `-o, --output <file>`: output file path. Default is `./output.tsv`.
- `-t, --threads <n>`: number of threads (currently reserved). Default is `1`.
- `-i, --ignore <level>`: filtering level from `0` to `4`. Default is `4`.

### Ignore Levels

The `ignore` option controls how gene entries are filtered:

| Level | Description |
|-------|-------------|
| 0 | Keep all entries. |
| 1 | Exclude genes missing chromosome or position information. |
| 2 | Additionally exclude entries whose chromosome name does not start with `chr`. |
| 3 | Additionally exclude `chrX`, `chrY` and `chrM`. |
| 4 | Keep only chromosomes `1–22` and remove the `chr` prefix. |

## Output

The result is a tab-separated file. The first line begins with `#` and lists the column headers:

```
#Chr    start   end     gene_id ...
```

Subsequent lines contain gene expression values with location information, sorted by chromosome and start coordinate.

## Example

```bash
GTFTPMCombinator -g genes.gtf -m expression.tsv -o combined.tsv -i 4
```

The above command reads `genes.gtf` and `expression.tsv`, filters and sorts them, and writes the result to `combined.tsv`.

## License

This project is available under the MIT license. See the `LICENSE` file if present.
