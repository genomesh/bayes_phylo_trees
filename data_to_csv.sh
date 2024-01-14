cd ~/dev/stori/data/new_housekeeping
for gene_id in {2..3408}
do
    cd "gene_$gene_id"
    if [ $(ls | wc -l) -gt 1 ]
    then
        echo "gene_$gene_id"
        Rscript ~/dev/stori/data_to_csv.R $gene_id
    fi
    cd ..
done