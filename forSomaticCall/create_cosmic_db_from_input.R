args = commandArgs(trailingOnly=TRUE)

input_f<-args[1]
output_f <-args[2]


write_somatic_info_from_biomart <- function(mart_obj,somatic_attributes,chr_start_end_list,output_f)
{
    cse_filter <- c('chromosome_name','start','end')
    bm <- getBM(attributes=somatic_attributes,filters=cse_filter,values=chr_start_end_list,mart=mart_obj)
    print(head(bm))
    write.table(bm,file=output_f,row.names=FALSE,col.names=FALSE,sep='\t', quote=FALSE,append=TRUE)
}


somatic_input <- read.table(input_f, header=T, sep='\t')
library('biomaRt')
ensembl=useMart('ENSEMBL_MART_ENSEMBL',dataset='hsapiens_gene_ensembl')
soma_att <- c("somatic_variation_name","somatic_source_description","somatic_clinical_significance","somatic_chromosome_start","somatic_chromosome_end","somatic_synonymous_status")

for(i in 1:nrow(somatic_input)){
    if(i==1){
        soma_att_list<-as.list(soma_att)
        write.table(soma_att_list,file=output_f,row.names=FALSE,col.names=FALSE,sep='\t', quote=FALSE)
    }
    chr<-somatic_input[i,]$CHROM
    pos<-somatic_input[i,]$POS
    chr_num <- gsub('chr','',chr)
    l <- list(chr_num,pos,pos)
    write_somatic_info_from_biomart(ensembl,soma_att,l,output_f)
}

