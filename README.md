# Freyja-SRA

[![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.04.3-brightgreen.svg)](http://nextflow.io) ![](https://img.shields.io/docker/image-size/dylanpilz/freyja-sra/latest)


Automated SRA downloading, processing and [Freyja](https://github.com/andersen-lab/Freyja) analysis pipeline.

## Installation
### Local Install via Git

```bash
git clone https://github.com/dylanpilz/Freyja-SRA.git
cd Freyja-SRA
```

## Usage
```bash
nextflow run main.nf -entry [sra|rerun_demix] -profile [docker|singularity] --accession_list [accession_list.csv] --output_dir [output_dir] --num_samples [num_samples]
```
### Parameters
* `-entry` - The pipeline entry point. 

    * `sra` will download, process and run Freyja on the provided SRA accessions.
        * `--accession_list` - A CSV file containing a list of SRA accessions to download and process. The CSV file should have a header row and the first column should be named `accession`.

    * `rerun_demix` will run freyja demix step on previously generated variants output files in the provided variants directory. This is useful if you want to run Freyja on existing data with a different barcode set.
        * `--variants_dir` must contain files in the format `[base_name].variants.tsv [base_name].depths.tsv` for each sample.

    * `--output_dir` - The final output directory. Creates `variants`, `demix`, and `covariants` subdirectories containing respective output files. (default: `./outputs`)
    
    * `--num_samples` - The number of samples to process. (default: 200)

## Configuration

Addtional configuration options can be found in `nextflow.config`