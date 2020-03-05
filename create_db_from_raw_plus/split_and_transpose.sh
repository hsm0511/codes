task_path=$1
gts_f=`ls $task_path | grep output`

for result in $gts_f
    do
        result=$task_path$result
        echo $result
        # split files with line 100k with header for each file
        tail -n +2 $result | split -l 100000 - GTsplit_
        for file in GTsplit_*
        do
            head -n 1 $result > tmp_file
                   echo "spliting"
            cat "$file" >> tmp_file
            mv -f tmp_file "$file"
        done
        # transpose file, make tsv 
        for file in GTsplit_a*
        do
            /mnt/clinix/Tools/Dev/python/bin/python /clinix/Analysis/Projects/PersonalGenomics_RnD/add_DB_from_RawPlus/Scripts/transpose.py "$file" "$result" "$task_path"
            echo "converting"
        done
        rm GTsplit_*
        rm $result
    done
