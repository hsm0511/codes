import vcf
import sys
from multiprocessing import *
from pymongo import MongoClient
import requests
from bs4 import BeautifulSoup as bs
import subprocess
import time

def create_cosmic_DB_by_biomart(input_f,cosmic_db_name):
    subprocess.call(['C:/Program Files/R/R-3.6.3/bin/x64/Rscript','--vanilla','create_cosmic_db_from_input.R',input_f,cosmic_db_name])

def check_quote_exists(soup):
    soup_q=soup.find('p',class_='quote')
    if str(soup_q)=='None':
        return False
    else:
        if 'has been flagged as a SNP' in str(soup_q.text):
            return True
        else:
            return 'Unknown'

def check_somatic_pathogenic(soup):
    soup_div=soup.find('div',class_='section-content')
    somatic_check=''
    pathogenic_check=''
    if str(soup_div)=='None':
        quote_check=check_quote_exists(soup)
        if quote_check==True:
            return 'germline','-'
        else:
            return 'Unknown','-'
    else:
        soup_dt_dd=soup.find_all(['dt','dd'])
        for i in range(len(soup_dt_dd)):
            content=soup_dt_dd[i]
            if 'Ever confirmed' in content.text:
                somatic_check=soup_dt_dd[i+1].text
            elif 'FATHMM prediction' in content.text:
                pathogenic_check=soup_dt_dd[i+1].text
                
        if somatic_check=='Yes':
            pathogenic_check=pathogenic_check.strip()
            pathogenic_check=pathogenic_check.replace(' ','_')
            pathogenic_check=pathogenic_check.replace('\n',' ')
            return 'somatic',pathogenic_check
        else:
            print('weird answer : '+somatic_check)
            sys.exit(0)
            return somatic_check,pathogenic_check.strip()

def get_mutation_info_from_cosmic(cosmic_num):
    cosmic_id=cosmic_num[4:]
    cosmic_info={'mutation':'','pathogenicity':''}
    if cosmic_num=='':
        return cosmic_info
    c= requests.Session()
    cosmic_url='https://cancer.sanger.ac.uk/cosmic/mutation/overview?id='+cosmic_id
    # get html source or get error message
    try:
        cosmic_source=c.get(cosmic_url,headers = {'User-agent': 'your bot 0.1'})
        soup=bs(cosmic_source.content,'html.parser')
        cosmic_info['mutation'],cosmic_info['pathogenicity']=check_somatic_pathogenic(soup)
        
        return cosmic_info
    except:
        print(cosmic_num+' : Cannot get URL source.')
        cosmic_info['error']='get url error'
        return cosmic_info

def check_somatic_by_freq(pop_freq,var_freq,pop_ref,var_ref):
    somatic_score=0
    var_freq=float(var_freq)
    if pop_freq=='None':
        if var_freq<0.1:
            somatic_score+=20
        return 'Unknown',somatic_score
    pop_freq=float(pop_freq)
    if pop_ref==var_ref:
        if pop_freq<0.1 and var_freq<0.1:
            somatic_score+=20
            return 'somatic',somatic_score
        else:
            return 'germline',somatic_score
    else:
        if pop_freq>0.9 and var_freq<0.1:
            somatic_score+=20
            return 'somatic', somatic_score
        else:
            return 'germline', somatic_score

def get_population_dict(pop_file_name):
    pop_dict={}
    pf=open(pop_file_name,'r')
    for line in pf:
        if line.rstrip()!='':
            l=line.rstrip().split('\t')
            if 'PopFreqMax_gnomAD.ExAC.1000g' in line:
                chr_idx=l.index('CHROM')
                pos_idx=l.index('POS')
                ref_idx=l.index('REF')
                freq_idx=l.index('PopFreqMax_gnomAD.ExAC.1000g')
            else:
                chr=l[chr_idx].split('chr')[-1]
                pos=l[pos_idx]
                freq=l[freq_idx]
                ref=l[ref_idx]
                pop_dict[chr+'\t'+pos]=[freq,ref]
    pf.close()
    return pop_dict

def get_cosmic_dict(cosmic_file_name):
    cosmic_dict={}
    cf=open(cosmic_file_name,'r')
    for line in cf:
        if line.rstrip()!='':
            l=list(map(lambda x:x.strip(),line.strip('\n').split('\t')))
            if 'cosmic.v84.from.vcf' in line:
                print(l)
                chr_idx=l.index('CHROM')
                pos_idx=l.index('POS')
                cosmic_idx=l.index('cosmic.v84.from.vcf')
            else:
                chr=l[chr_idx].split('chr')[-1]
                pos=l[pos_idx]
                cosmic_id=list(map(lambda x:x.strip(),l[cosmic_idx].split(',')))
                cosmic_dict[chr+'\t'+pos]=cosmic_id
    cf.close()
    return cosmic_dict

def get_pop_freq_from_1000genome(chr_pos,alt,col_1000g):
    chr,pos=chr_pos.split('\t')
    find_1000g=col_1000g.find_one({'chr':chr,'pos':pos},{'AF':1,'alt':1})
    if str(find_1000g)=='None':
        return 'None'
    else:
        alt_l=find_1000g['alt']
        if alt in alt_l:
            alt_idx=alt_l.index(alt)
            allele_freq=find_1000g['AF'][alt_idx]
            return allele_freq
        else:
            return 'None'

def get_variant_from_vcf(vcf_name):
    variant_info={}
    vf=open(vcf_name,'r')
    for line in vf:
        if line.rstrip()!='':
            l=line.rstrip().split('\t')
            if 'HGVSc' in line and 'VAF' in line and 'Depth' in line:
                chr_idx=l.index('CHROM')
                pos_idx=l.index('POS')
                ref_idx=l.index('REF')
                alt_idx=l.index('ALT')
                vaf_idx=l.index('VAF')
            else:
                chr=l[chr_idx].split('chr')[-1]
                pos=l[pos_idx]
                ref=l[ref_idx]
                alt=l[alt_idx]
                vaf=l[vaf_idx]
                variant_info[chr+'\t'+pos]=[vaf,ref,alt]
    vf.close()
    return variant_info

def get_cosmic_id_mut_from_db(cosmic_db_file):
    cosmic_id_somatic=[]
    f=open(cosmic_db_file,'r')
    for line in f:
        if line.rstrip()!='':
            l=line.strip('\n').split('\t')
            print(line)
            cosmic_id=l[0]
            mutation_info=l[1]
            if cosmic_id=='' or cosmic_id=='NA':
                continue
            elif 'somatic' in mutation_info or 'Somatic' in mutation_info:
                cosmic_id_somatic.append(cosmic_id)
    f.close()
    return list(set(cosmic_id_somatic))

def main():
    
    mongo = MongoClient('211.192.85.209', 27017)
    db_2=mongo['public_db']
    col_1000g=db_2['1000genome']
    
    input_vcf_txt=sys.argv[1]
    pop_file=sys.argv[2]
    cosmic_file=sys.argv[3]
    cosmic_db_file =sys.argv[4]

    # create cosmic ID Database by biomart
    create_cosmic_DB_by_biomart(input_vcf_txt,cosmic_db_file)
    cosmic_id_somatic=get_cosmic_id_mut_from_db(cosmic_db_file)
    print(cosmic_id_somatic)
    
    # get information from input file
    pop_freq_dict=get_population_dict(pop_file)
    cosmic_id_dict=get_cosmic_dict(cosmic_file)
    vaf_dict=get_variant_from_vcf(input_vcf_txt)

    output_f=open('variant_patient_info.txt','w')
    header='CHROM	POS\tVariant_type\tPathogenicity\n'
    output_f.write(header)

    for chr_pos in vaf_dict:
        somatic_score_by_cosmic=[]
        somatic_score_by_freq=0
        mut_info=''
        patho_info=''
        variant_freq,v_ref,v_alt=vaf_dict[chr_pos]
        
        # get somatic score by frequency
        if chr_pos not in pop_freq_dict:
            pop_freq='None'
            pop_ref='None'
            pop_freq_from_1000g=get_pop_freq_from_1000genome(chr_pos,v_alt,col_1000g)
            if pop_freq_from_1000g!='None':
                pop_freq=pop_freq_from_1000g
        else:
            pop_freq,pop_ref=pop_freq_dict[chr_pos]
        type_from_pop,somatic_score_by_freq=check_somatic_by_freq(pop_freq,variant_freq,pop_ref,v_ref)

        # get somatic score by cosmic id
        patho_l=['']
        if chr_pos not in cosmic_id_dict or cosmic_id_dict[chr_pos]==['']:
            pass
        else:
            cosmic_num_l=cosmic_id_dict[chr_pos]
            type_info_from_cosmic={'mutation':[],'pathogenicity':[]}
            for cos_id in cosmic_num_l:
                mutation_type,patho_type=['','']
                info=get_mutation_info_from_cosmic(cos_id)
                if 'error' in info:
                    print(chr_pos)
                    print(info['error'])
                else:
                    mutation_type=info['mutation']
                    if cos_id in cosmic_id_somatic or mutation_type=='somatic':
                        somatic_score_by_cosmic.append(10)
                    patho_type=info['pathogenicity']
                    patho_type=patho_type.replace('_','')
                type_info_from_cosmic['mutation'].append(mutation_type)
                type_info_from_cosmic['pathogenicity'].append(patho_type)
            mut_l=type_info_from_cosmic['mutation']
            patho_l=type_info_from_cosmic['pathogenicity']
        
        # calculate final somatic score
        if 10 in somatic_score_by_cosmic:
            somatic_score_by_cosmic=10
        else:
            somatic_score_by_cosmic=0
        somatic_score_final=somatic_score_by_freq+somatic_score_by_cosmic
        if somatic_score_final>=20:
            mut_info='somatic'
        else:
            mut_info='germline'
        patho_l=list(set(patho_l))
        patho_info=','.join(patho_l)
        output_f.write(chr_pos+'\t'+mut_info+'\t'+patho_info+'\n')
    output_f.close()


main()






