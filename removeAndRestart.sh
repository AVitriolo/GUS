dirs="results
logs
data
resources/reference_data/gencode
resources/reference_data/mappability/*.filtered.bedgraph
resources/RNA/*_counts_*_*.tsv
resources/composite
resources/checks"

tokeep="resources/TxIDs/cancerGeneList_processed.tsv
resources/TxIDs/cancerGeneList.tsv
resources/TxIDs/selected_TxIDs"

for i in $tokeep; do
	mv $i .
done

rm -rf $dirs
rm resources/TxIDs/*

