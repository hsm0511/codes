merge_RawPlus_by_callRate_limit.py : collect genotype files of sample with high quality
get_raw_plus_and_merge.py : merge samples' genotype files by multiprocessing
merge_genotype.py : merge file by pandas module
split_and_transpose.py : split files not to exceed maximum BSON size of Mongo document
transpose.py : transpose the content in file by pandas module
extend_tsv.py : extend newly created tsv to existing tsv
mongo_insert_GT.sh : insert genotypes to DB by mongoimport tsv
split_colls_time_avoiding.py : split genotypes to 350 collections