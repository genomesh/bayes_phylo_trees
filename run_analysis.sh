cp *.fasta gene.fasta;
Rscript ~/dev/stori/convert_fasta_to_nex.R gene.fasta;
~/dev/stori/iqtree/bin/iqtree -s gene.fasta -m GTR+F+I+G4;
sed -i -e 's/|/_/g' gene.nex;
mb ~/dev/stori/mb_script.nex;