# for accession in data/samples_to_rerun.csv, copy {accession}.variants.tsv and {accession}.depths.tsv demix_rerun/

while read -r accession; do
    gsutil stat gs://outbreak-ww-data/variants/${accession}.variants.tsv
    if [ $? -ne 0 ]; then
        echo "gs://outbreak-ww-data/variants/${accession}.variants.tsv does not exist"
        continue
    fi
    gcloud storage cp gs://outbreak-ww-data/variants/${accession}.variants.tsv demix_rerun/
    gcloud storage cp gs://outbreak-ww-data/variants/${accession}.depths.tsv demix_rerun/
done < "data/samples_to_rerun.csv"