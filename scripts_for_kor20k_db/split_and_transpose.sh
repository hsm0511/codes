gts_f=`ls | grep output`

for result in $gts_f
    do
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
            python transpose.py "$file" "$result"
            echo "converting"
        done
        rm GTsplit_*
    done
