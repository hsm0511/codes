task_path=$1
gts_f=`ls $task_path | grep _T.tsv`

for result in $gts_f
    do
	# replace dot character to underbar in probe names(header)
	sed '1!b;s/\./_/g' $result > tmp_dot_replaced
	mv tmp_dot_replaced $result
	echo $result
	# insert tsv file to mongo DB collection
        mongoimport --host 172.19.87.50 --port 27017 -d Repository -c tmp_raw_gt --type tsv --file $result --headerline
    done
