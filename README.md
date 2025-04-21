# Freyja-SC2

[![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.04.3-brightgreen.svg)](http://nextflow.io) [![Process samples](https://github.com/andersen-lab/freyja-sc2/actions/workflows/process_samples.yml/badge.svg)](https://github.com/andersen-lab/freyja-sc2/actions/workflows/process_samples.yml)


Automated SRA downloading, processing and [Freyja](https://github.com/andersen-lab/Freyja) analysis pipeline for SARS-CoV-2 wastewater sequencing data.

## Installation
```bash
git clone https://github.com/dylanpilz/freyja-sc2.git
cd freyja-sc2
```

## Usage
```bash
nextflow run main.nf -entry [fetch|rerun_demix] -profile [docker|singularity] --accession_list [accession_list.csv] --output_dir [output_dir] --num_samples [num_samples]
```
### Parameters
* `-entry` - The pipeline entry point. 

    * `fetch` will download, process and run Freyja on the provided SRA accessions.
        * `--sra_metadata` - A CSV file containing a list of SRA accessions to download and process. The CSV file should have a header row, containing at least the following metadata fields:
            * `accession` - The SRA accession number (e.g. SRR1234567).
            * `sample_status` - The status of the sample (e.g. "to_run", "completed").
            * `amplicon_PCR_primer_scheme` - The amplicon primer scheme used for the sample (e.g. "ARTIC V4", ARTIC V5.3.2, unknown etc.)
        * `--num_samples` - The number of samples to process per run. (default: 100)

    * `rerun_demix` will run freyja demix step on previously generated variants output files in the provided variants directory. This is useful if you want to run Freyja on existing data with a different barcode set.
        * `--variants_dir` must contain files in the format `[base_name].variants.tsv [base_name].depths.tsv` for each sample.

    * `--output_dir` - The final output directory. Creates `variants` and `demix`  subdirectories containing respective output files. (default: `./outputs`)
    
## Configuration

Addtional configuration options can be found in `nextflow.config`

## Data Availability

freyja-sc2 is currently in the process of downloading and processing all publicly available SARS-CoV-2 wastewater data, fetched with the following search terms:
```
'(Wastewater[All Fields] OR wastewater metagenome[All Fields]) AND ("Severe acute respiratory syndrome coronavirus 2"[Organism] OR SARS-CoV-2[All Fields])

```
In addition, to the above search terms, we exclude accessions that don't meet the following metadata requirements:
* Missing collection date
* Missing catchment size (`ww_population`)
* Missing location (`geo_loc_name`)

To check the status of each accession, please refer to the `sample_status` column in `data/all_metadata.csv`. All currently processed freyja outputs are publicly available via Google Cloud Storage at `gs://outbreak-ww-data`
