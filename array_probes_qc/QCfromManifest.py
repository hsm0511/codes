import re
import subprocess
from multiprocessing import Pool
import pymongo
import vcf

complement_pair={'A':'T','T':'A','C':'G','G':'C',
                'a':'t','t':'a','c':'g','g':'c','':''}

def get_info_from_manifest(manifest_csv):
    probes_manifest_info={}
    header_found=False
    for line in manifest_csv:
        l=line.rstrip().split(',')
        if line.rstrip()!='':
            if 'AlleleA_ProbeSeq' in line:
                header_found=True
                name_idx=l.index('Name')
                chr_idx=l.index('Chr')
                pos_idx=l.index('MapInfo')
                probeseq_a_idx=l.index('AlleleA_ProbeSeq')
                sourceseq_idx=l.index('SourceSeq')
                snp_idx=l.index('SNP')
                max_info_idx=max([name_idx,chr_idx,pos_idx,probeseq_a_idx,sourceseq_idx,snp_idx])
            elif header_found==False:
                pass
            elif line.rstrip()=='[Controls]':
                break
            else:
                if len(l)<max_info_idx+1:
                    continue
                probe_name=l[name_idx]
                snp_l=l[snp_idx][1:-1].split('/')
                sourceseq=l[sourceseq_idx]
                probe_seq_a=l[probeseq_a_idx]
                chr_num=l[chr_idx]
                pos=l[pos_idx]
                if chr_num not in probes_manifest_info:
                    probes_manifest_info[chr_num]={}
                probes_manifest_info[chr_num][probe_name]={'probe_seq_a':probe_seq_a,'snp':snp_l,'source_seq':sourceseq,'pos':pos}
    return probes_manifest_info


def calcluate_ngs_pon(vcf_list_and_probe_info):
    vcf_file_list=vcf_list_and_probe_info['vcf']
    chr_pos_probe=vcf_list_and_probe_info['probe']
    chr_pos_ngs_pon={}
    for vcf_file in vcf_file_list:
        vcf_reader = vcf.Reader(open(vcf_file, 'r'))
        samples = vcf_reader.samples
        for rec in vcf_reader:
            chr_pos_vcf=str(rec.CHROM)+'\t'+str(rec.POS)
            #cnt += 1
            if chr_pos_vcf in chr_pos_probe:
                if chr_pos_vcf not in chr_pos_ngs_pon:
                    chr_pos_ngs_pon[chr_pos_vcf]={'Ref':str(rec.REF),'genotypes':{}}
                genotypes = []
                for s in samples:
                    if str(rec.FILTER)=='[]':
                        genotypes.append(rec.genotype(s).gt_bases)
                    else:
                        genotypes.append('./.')
                for gt in genotypes:
                    if gt not in chr_pos_ngs_pon[chr_pos_vcf]['genotypes']:
                        chr_pos_ngs_pon[chr_pos_vcf]['genotypes'][gt]=0
                    chr_pos_ngs_pon[chr_pos_vcf]['genotypes'][gt]+=1
                    if '75012985' in chr_pos_vcf:
                        print(chr_pos_ngs_pon[chr_pos_vcf])
    return chr_pos_ngs_pon


class InfiniumDesign(object):
    def __init__(self,manifest_csv):
        self.manifest=manifest_csv
    
    @staticmethod
    def reverse_complement(sequence):   # make sequence reverse complement
        seq_rv=sequence[::-1]   # reverse
        seq_rv_l=list(seq_rv)
        seq_rv_cp_l=map(lambda x:complement_pair[x],seq_rv_l)   # complement
        seq_rv_cp=''.join(seq_rv_cp_l)
        return seq_rv_cp

    @staticmethod
    def sourceSeq_match_and_inf_design(probe_seq_a,source_seq): # decide infinium design and match probeSeqA with sourceseq
        inf_design='cannot decide'
        flag=False  # flag for probeSeqA and sourseq match
        
        # reverse complement of probeSeq
        probe_seq_rv_cp_a=InfiniumDesign.reverse_complement(probe_seq_a)
        
        source_seq_half_l=[]
        snp_start_idx=source_seq.index('[')
        snp_last_idx=source_seq.index(']')

        # source seq half without snp
        source_seq_half_l.append(source_seq[:snp_start_idx])
        source_seq_half_l.append(source_seq[snp_last_idx+1:])
        snp=source_seq[snp_start_idx:snp_last_idx+1]
        snp_l=source_seq[snp_start_idx+1:snp_last_idx].split('/')
        
        # source seq half with snp
        source_seq_with_snp_l=[]
        for base in snp_l:
            source_seq_with_snp_l.append(source_seq_half_l[0]+base)
            source_seq_with_snp_l.append(base+source_seq_half_l[1])
        
        # design and source seq half
        check_source_seqs={'infinium_2_match':source_seq_half_l,'infinium_1_match':source_seq_with_snp_l}

        # decide infinium design by probeSeqA and sourceSeq match
        for design in check_source_seqs:
            for i in range(len(check_source_seqs[design])):
                check_seq=check_source_seqs[design][i]
                if design=='infinium_2_match' and i==0 or design=='infinium_1_match' and i in [0,2]:    # match probeSeqA with sourceSeq left half
                    a=re.compile(probe_seq_a+'$',re.I)
                    if probe_seq_a!='' and a.search(check_seq):
                        flag=True
                        inf_design=design[:design.index('_match')]
                        break
                else:   # match probeSeqA with soruceSeq right half
                    rv_a=re.compile(probe_seq_rv_cp_a,re.I)
                    if probe_seq_rv_cp_a!='' and rv_a.match(check_seq):
                        flag=True
                        inf_design=design[:design.index('_match')]
                        break
            if flag==True:  # infinium design decide complete and probeSeqA matches sourceSeq
                break
        return flag, inf_design
    

    def get_probes_infinium_design(self,probes_manifest_info={}):
        probe_design={}
        if probes_manifest_info=={}:
            probes_manifest_info=get_info_from_manifest(self.manifest)
        for chr_num in probes_manifest_info:
            for probe in probes_manifest_info[chr_num]:
                probe_seq_a=probes_manifest_info[chr_num][probe]['probe_seq_a']
                source_seq=probes_manifest_info[chr_num][probe]['source_seq']
                match_flag,infinium_design=InfiniumDesign.sourceSeq_match_and_inf_design(probe_seq_a,source_seq)
                probe_design[probe]=infinium_design
        return probe_design

class VariantNearTarget(object):
    
    def __init__(self,manifest_csv):
        self.manifest=manifest_csv

    @staticmethod
    def direction_from_target(probe_seq_a,source_seq):  # return direction of probeSeq from target SNP and flag whether target SNP is included in probeSeq
        direction=''
        target_included=False
        check_source_seqs=[source_seq]
        probe_seq_rv_cp_a=InfiniumDesign.reverse_complement(probe_seq_a)
        snp_start_idx=source_seq.index('[')
        snp_last_idx=source_seq.index(']')
        snp=source_seq[snp_start_idx:snp_last_idx+1]
        snp_l=source_seq[snp_start_idx+1:snp_last_idx].split('/')
        source_seq_with_snp_l=[]
        for base in snp_l:
            source_seq_modified=source_seq.replace(snp,base)
            source_seq_with_snp_l.append(source_seq_modified)
        # get list for source seq without target SNP and source seq with target SNP
        check_source_seqs=[source_seq]
        check_source_seqs.extend(source_seq_with_snp_l)
        for i in range(len(check_source_seqs)):
            check_seq=check_source_seqs[i]
            # minus direction if probeSeq A matches sourceSeq
            if re.search(probe_seq_a,check_seq,flags=re.I):
                direction='minus'
                break
            # plus direction if reverse complemented probeSeq A matches sourceSeq
            elif re.search(probe_seq_rv_cp_a,check_seq,flags=re.I):
                direction='plus'
                break
            else:
                continue
        if i==0:    # probeSeq or probeSeq_rv_cp matches sourceSeq without target SNP
            pass
        else:   # probeSeq or probeSeq_rv_cp matches sourceSeq with target SNP
            target_included=True
        return direction, target_included

    def get_probe_variant_info(self):

        # search variant from 1000genome collection from public DB in 50 server
        connection = pymongo.MongoClient('172.19.87.50', 27017)
        db_1000g = connection.public_db
        coll_1000g = db_1000g['1000genome']
        
        probe_variant={}
        
        # get manifest info
        probes_manifest_info=get_info_from_manifest(self.manifest)
        for chr in probes_manifest_info:
            for probe_name in probes_manifest_info[chr]:
                pos=int(probes_manifest_info[chr][probe_name]['pos'])
                probe_seq_a=probes_manifest_info[chr][probe_name]['probe_seq_a']
                pb_seq_a_len=len(probe_seq_a)
                source_seq=probes_manifest_info[chr][probe_name]['source_seq']

                # get direction and target inclusion info for the probe
                direction, target_include_flag=VariantNearTarget.direction_from_target(probe_seq_a,source_seq)

                # set range of position to search variant within probe sequence
                if target_include_flag==False:
                    find_pos_range=[pos-pb_seq_a_len,pos+pb_seq_a_len]
                else:
                    find_pos_range=[pos-(pb_seq_a_len-1),pos+(pb_seq_a_len-1)]
                
                # search variant within probeSeq in DB 1000genome collection
                info_1000g=coll_1000g.find({"chr":str(chr),"pos":{'$gt':find_pos_range[0]-1,'$lt':find_pos_range[1]+1}}, {"ref":1,"alt":1,"VT":1,"EAS":1,'pos':1})
                variant_info={}
                if info_1000g=='None':  # no variant in probeSeq
                    probe_variant[probe_name]='No Pos range in 1000g'
                    continue
                else:   # found variants in probeSeq
                    for info in info_1000g:
                        info_pos=info['pos']-pos    # distance from target SNP
                        # skip variant same as target SNP
                        if info_pos ==0:
                            continue
                        # plus direction variant found
                        elif info_pos > 0:
                            info_pos='+'+str(info_pos)
                            if direction=='minus':  # skip if direction of probeSeq is minus
                                continue
                        # minus direction variant found
                        else:
                            info_pos=str(info_pos)
                            if direction=='plus':   # skip if direction of probeSeq is plus
                                continue
                        # get EAS allele freq. of variant found
                        info_eas=info['EAS']
                        # insert variant distance and allele freq info in dict
                        for freq in info_eas:
                            if info_pos not in variant_info:
                                variant_info[str(info_pos)]=[]
                            variant_info[str(info_pos)].append(str(freq))

                # no variant found in probeSeq except target SNP
                if len(variant_info.keys())==0:
                    probe_variant[probe_name]='No variant found'
                    continue

                # sort distances of variants in order of absolute value
                int_sorted_vps=list(map(lambda x : int(x),variant_info.keys()))
                int_sorted_vps.sort()
                if direction=='plus':
                    sorted_vps=list(map(lambda x : '+'+str(x),int_sorted_vps))
                else:
                    sorted_vps=list(map(lambda x : str(x),int_sorted_vps))
                    sorted_vps=sorted_vps[::-1]
                
                # write variant info with allele freq info
                freq_l=[]
                for vp in sorted_vps:
                    freq_l.append(vp+' ('+','.join(variant_info[vp])+')')
                probe_variant[probe_name]=' / '.join(freq_l)

        return probe_variant

class Indel(object):
    @staticmethod
    def get_probes_indel_info(probes_manifest_info):
        probe_indel={}
        for chr in probes_manifest_info:
            for probe_name in probes_manifest_info[chr]:
                snp=probes_manifest_info[chr][probe_name]['snp']
                if 'I' in snp or 'D' in snp:    # write 'indel' if I or D in SNP
                    indel_info='indel'
                probe_indel[probe_name]=indel_info
        return probe_indel

class GenotypeCountNGS(object):
    def __init__(self,manifest_csv,vcf_file_path,probe_ref_info):
        self.manifest=manifest_csv
        self.vcf_paths=vcf_file_path
        self.probe_ref=probe_ref_info
    
    @staticmethod
    def get_chr_pos_probe(manifest_csv):
        probe_manifest_info=get_info_from_manifest(manifest_csv)
        chr_pos_probe={}
        for chr in probe_manifest_info:
            for probe in probe_manifest_info[chr]:
                pos=str(probe_manifest_info[chr][probe]['pos'])
                cp='chr'+chr+'\t'+pos
                if cp not in chr_pos_probe:
                    chr_pos_probe[cp]=[]
                chr_pos_probe[cp].append(probe)
        return chr_pos_probe


    @staticmethod
    def split_list_for_mp(raw_list,number_of_cores):
        split_list=[]
        if len(raw_list)<number_of_cores:
            number_of_cores=len(raw_list)
        q=len(raw_list)/number_of_cores
        rem=len(raw_list)%number_of_cores
        num_split_l=[q]*number_of_cores
        for i in range(rem):
            num_split_l[i]+=1
        start_idx=0
        for num in num_split_l:
            last_idx=start_idx+num
            split_l=raw_list[start_idx:last_idx]
            split_list.append(split_l)
            start_idx=last_idx
        return split_list

    @staticmethod
    def hetero_merge(gt_ct_dict):
        merged_gt=gt_ct_dict.keys()
        merged_gt_ct_dict={}
        # merge hetero genotypes into one
        for gt in gt_ct_dict:
            gt_l=gt.split('/')
            if gt_l[0]!=gt_l[1] and gt_l[1]+'/'+gt_l[0] in merged_gt:
                merged_gt.remove(gt_l[1]+'/'+gt_l[0])
                break
        for gt in merged_gt:
            gt_l=gt.split('/')
            if gt_l[0]==gt_l[1]:
                merged_gt_ct_dict[gt]=gt_ct_dict[gt]
            else:
                another_gt=gt_l[1]+'/'+gt_l[0]
                if another_gt in gt_ct_dict.keys():
                    merged_gt_ct_dict[gt]=gt_ct_dict[gt]+gt_ct_dict[another_gt]
                else:
                    merged_gt_ct_dict[gt]=gt_ct_dict[gt]
        return merged_gt_ct_dict

    def get_probe_ngs_pon(self):
        chr_pos_probe=GenotypeCountNGS.get_chr_pos_probe(self.manifest)
        vcf_file_list=[]
        for line in self.vcf_paths:
            if line.rstrip()!='':
                vcf_file_list.append(line.rstrip())
        vcf_split_file_list=GenotypeCountNGS.split_list_for_mp(vcf_file_list,48)
        vcf_and_probe_info_list=[]
        for v_l in vcf_split_file_list:
            info_dict={'vcf':v_l,'probe':chr_pos_probe}
            vcf_and_probe_info_list.append(info_dict)
        pool=Pool(len(vcf_split_file_list))
        result_pon_l=pool.map(calcluate_ngs_pon,vcf_and_probe_info_list)
        pool.close()
        pool.join()
        merged_chr_pos_ngs_pon={}
        for pon_dict in result_pon_l:
            for cp in pon_dict:
                if cp not in merged_chr_pos_ngs_pon:
                    merged_chr_pos_ngs_pon[cp]={'Ref':pon_dict[cp]['Ref'],'genotypes':{}}
                for gt in pon_dict[cp]['genotypes']:
                    if gt not in merged_chr_pos_ngs_pon[cp]['genotypes']:
                        merged_chr_pos_ngs_pon[cp]['genotypes'][gt]=0
                    merged_chr_pos_ngs_pon[cp]['genotypes'][gt]+=pon_dict[cp]['genotypes'][gt]
        for cp in merged_chr_pos_ngs_pon:
            merged_chr_pos_ngs_pon[cp]['genotypes']=GenotypeCountNGS.hetero_merge(merged_chr_pos_ngs_pon[cp]['genotypes'])

        probe_ngs_pon={}
        for cp in chr_pos_probe:
            for probe in chr_pos_probe[cp]:
                manifest_ref_allele=self.probe_ref[probe]
                manifest_ref_homo_gt=manifest_ref_allele+'/'+manifest_ref_allele
                if cp not in merged_chr_pos_ngs_pon:
                    if len(manifest_ref_allele)==1:
                        probe_ngs_pon[probe]={manifest_ref_homo_gt:1950}
                    else:
                        probe_ngs_pon[probe]={'-':'-'}
                else:
                    sum_count=sum(merged_chr_pos_ngs_pon[cp]['genotypes'].values())
                    if sum_count<1950:
                        ngs_ref_allele=merged_chr_pos_ngs_pon[cp]['Ref']
                        ngs_ref_homo_gt=ngs_ref_allele+'/'+ngs_ref_allele
                        ref_homo_count=1950-sum_count
                        merged_chr_pos_ngs_pon[cp]['genotypes'][ngs_ref_homo_gt]=ref_homo_count
                        probe_ngs_pon[probe]=merged_chr_pos_ngs_pon[cp]['genotypes']
                    else:
                        probe_ngs_pon[probe]=merged_chr_pos_ngs_pon[cp]['genotypes']


        return probe_ngs_pon

