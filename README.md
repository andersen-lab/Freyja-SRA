# Freyja-SRA

[![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.04.3-brightgreen.svg)](http://nextflow.io) ![](https://img.shields.io/docker/image-size/dylanpilz/freyja-sra/latest) [![Process samples](https://github.com/andersen-lab/Freyja-SRA/actions/workflows/process_samples.yml/badge.svg)](https://github.com/andersen-lab/Freyja-SRA/actions/workflows/process_samples.yml)


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

## Data Availability

Freyja-SRA is currently in the process of downloading and processing all publicly available SRA data, fetched with the following search terms:
```
'(Wastewater[All Fields] OR wastewater metagenome[All Fields]) AND ("Severe acute respiratory syndrome coronavirus 2"[Organism] OR SARS-CoV-2[All Fields])

```
In addition, to the above search terms, we exclude accessions that don't meet the following metadata requirements:
* Missing collection date
* Missing catchment size (`ww_population`)
* Missing location (`geo_loc_name`)

To check the status of each accession, please refer to the `sample_status` column in `data/all_metadata.csv`. All currently processed freyja outputs are publicly available via Google Cloud Storage at `gs://outbreak-ww-data`
