# for accession in data/samples_to_rerun.txt, copy {accession}.variants.tsv and {accession}.depths.tsv demix_rerun/

while read -r accession; do
    gcloud storage cp gs://outbreak-ww-data/variants/${accession}.variants.tsv" "demix_rerun/"
    gcloud storage cp gs://outbreak-ww-data/variants/${accession}.depths.tsv" "demix_rerun/"
done < "data/samples_to_rerun.txt"