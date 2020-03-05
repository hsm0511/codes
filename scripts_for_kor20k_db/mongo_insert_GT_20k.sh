gts_f=`ls | grep _T.tsv`

for result in $gts_f
    do
	# replace dot character to underbar in probe names(header)
	sed '1!b;s/\./_/g' $result > tmp_dot_replaced
	mv tmp_dot_replaced $result
	echo $result
	# insert tsv file to mongo DB collection
        mongoimport -d Repository -c cathy_test --type tsv --file $result --headerline
    done
