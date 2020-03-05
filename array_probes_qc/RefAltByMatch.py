import re
import subprocess
from multiprocessing import Pool
from QCfromManifest import *

complement_pair={'A':'T','T':'A','C':'G','G':'C',
                'a':'t','t':'a','c':'g','g':'c','':''}


def ref_alt_match_tool(probes_manifest_info,ref_genome_path,probe_design,selected_chr_probe):
    ref_genome_by_chr=RefAltMatch.get_ref_genome_by_chr(ref_genome_path)
    probes_manifest_info_selected={}
    for selected_chr in selected_chr_probe:
        if selected_chr not in probes_manifest_info_selected:
            probes_manifest_info_selected[selected_chr]={}
        for selected_probe in selected_chr_probe[selected_chr]:
            if selected_probe not in probes_manifest_info_selected[selected_chr]:
                probes_manifest_info_selected[selected_chr][selected_probe]=probes_manifest_info[selected_chr][selected_probe]
    
    chr_sorted=probes_manifest_info_selected.keys()
    chr_sorted.sort()
    probe_ref_alt={}
    for chr_num in chr_sorted:
        if chr_num=='0':
            for probe in probes_manifest_info_selected[chr_num]:
                probe_ref_alt[probe]=['Chr is zero']*2
                print(probe+'\t'+'\t'.join(probe_ref_alt[probe]))
        elif chr_num=='XY':
            rfg_genome_X=ref_genome_by_chr['X']
            rfg_genome_Y=ref_genome_by_chr['Y']
            rfg_genomes=[rfg_genome_X,rfg_genome_Y]
            for probe in probes_manifest_info_selected[chr_num]:
                infinium_design=probe_design[probe]
                ref_from_refgenome_l=[]
                pseq_a=probes_manifest_info_selected[chr_num][probe]['probe_seq_a']
                snp_l=probes_manifest_info_selected[chr_num][probe]['snp']
                if 'I' in snp_l or 'D' in snp_l:
                    probe_ref_alt[probe]=['indel']*2
                else:
                    snp_l_cp=list(map(lambda x:complement_pair[x],snp_l))
                    ref_from_refgenome_l=list(map(lambda x:RefAltMatch.ref_from_probeseq_refgenome_match(pseq_a,infinium_design,x),rfg_genomes))

                    if "probeseqA doesn't match with ref genome" in ref_from_refgenome_l and len(set(ref_from_refgenome_l))==1:
                        ref_from_refgenome="probeseqA doesn't match with ref genome"
                        alt="probeseqA doesn't match with ref genome"
                    elif "probeseqA doesn't match with ref genome" in ref_from_refgenome_l and len(set(ref_from_refgenome_l))==2:
                        ref_found_idx=1-ref_from_refgenome_l.index("probeseqA doesn't match with ref genome")
                        ref_tmp=ref_from_refgenome_l[ref_found_idx]
                        ref_from_refgenome=ref_from_refgenome_l[ref_found_idx]+'('+['X','Y'][ref_found_idx]+')'
                        if ref_tmp in snp_l:
                            ref_idx=snp_l.index(ref_tmp)
                            alt=snp_l[1-ref_idx]
                        elif ref_tmp in snp_l_cp:
                            ref_idx=snp_l_cp.index(ref_tmp)
                            alt=snp_l_cp[1-ref_idx]
                        else:
                            ref_from_refgenome='Ref is not in SNP'
                            alt='Ref is not in SNP'
                    elif len(set(ref_from_refgenome_l))==1:
                        ref_from_refgenome=ref_from_refgenome_l[0]
                        if ref_from_refgenome in snp_l:
                            ref_idx=snp_l.index(ref_from_refgenome)
                            alt=snp_l[1-ref_idx]
                        elif ref_from_refgenome in snp_l_cp:
                            ref_idx=snp_l_cp.index(ref_from_refgenome)
                            alt=snp_l_cp[1-ref_idx]
                        else:
                            ref_from_refgenome='Ref is not in SNP'
                            alt='Ref is not in SNP'
                    elif len(set(ref_from_refgenome_l))==2:
                        ref_from_refgenome="Ref doesn't match between X and Y"
                        alt="Ref doesn't match between X and Y"
                    probe_ref_alt[probe]=[ref_from_refgenome,alt]
                print(probe+'\t'+'\t'.join(probe_ref_alt[probe]))
        else:
            rfg_genome=ref_genome_by_chr[chr_num]
            for probe in probes_manifest_info_selected[chr_num]:
                infinium_design=probe_design[probe]
                pseq_a=probes_manifest_info_selected[chr_num][probe]['probe_seq_a']
                snp_l=probes_manifest_info_selected[chr_num][probe]['snp']
                
                if 'I' in snp_l or 'D' in snp_l:
                    probe_ref_alt[probe]=['indel']*2
                else:
                    snp_l_cp=list(map(lambda x:complement_pair[x],snp_l))
                    ref_from_refgenome=RefAltMatch.ref_from_probeseq_refgenome_match(pseq_a,infinium_design,rfg_genome)
                    if ref_from_refgenome=="probeseqA doesn't match with ref genome":
                        alt="probeseqA doesn't match with ref genome"
                    else:
                        if ref_from_refgenome in snp_l:
                            ref_idx=snp_l.index(ref_from_refgenome)
                            alt=snp_l[1-ref_idx]
                        elif ref_from_refgenome in snp_l_cp:
                            ref_idx=snp_l_cp.index(ref_from_refgenome)
                            alt=snp_l_cp[1-ref_idx]
                        else:
                            ref_from_refgenome='Ref is not in SNP'
                            alt='Ref is not in SNP'
                    probe_ref_alt[probe]=[ref_from_refgenome,alt]
                print(probe+'\t'+'\t'.join(probe_ref_alt[probe]))
    return probe_ref_alt

def ref_alt_match_tool_wrapper(data):
    result=ref_alt_match_tool(*data)
    return result

class RefAltMatch(object):

    def __init__(self,manifest_csv,ref_genome_file_path,num_of_cores):
        self.manifest=manifest_csv
        self.ref_genome_path=ref_genome_file_path
        self.core_num=num_of_cores
    
    @staticmethod
    def get_ref_genome_by_chr(ref_genome_file_path):
        ref_genome_by_chr={}
        cmd='ls '+ref_genome_file_path+" | grep '^chr' | grep fa"
        cmd_result=subprocess.check_output(cmd,shell=True).rstrip()
        file_l=cmd_result.split('\n')
        for f in file_l:
            chr_num=f.split('_')[0][3:]
            ref_chr_file=open(ref_genome_file_path+f,'r')
            ref_genome_by_chr[chr_num]=ref_chr_file.readlines()[1].rstrip()
            ref_chr_file.close()
        return ref_genome_by_chr

    @staticmethod
    def ref_from_probeseq_refgenome_match(probe_seq_a,infinium_design,chr_refgenome):
        if infinium_design=='infinium_1':
            match_pseq=probe_seq_a[:-1]
        else:
            match_pseq=probe_seq_a

        match_pseq_rvcp=''.join(list(map(lambda x:complement_pair[x],list(match_pseq[::-1]))))
        search_result_1=re.search(match_pseq,chr_refgenome,flags=re.I)
        search_result_2=re.search(match_pseq_rvcp,chr_refgenome,flags=re.I)
        if search_result_1:
            target_snp_idx=search_result_1.end()
            ref=chr_refgenome[target_snp_idx]
        elif search_result_2:
            target_snp_idx=search_result_2.start()-1
            ref=chr_refgenome[target_snp_idx]
        else:
            ref="probeseqA doesn't match with ref genome"
        return ref

    @staticmethod
    def split_chr_probe(probes_manifest_info,num_of_cores):
        chr_probe=[]
        for chr in probes_manifest_info:
            for p in probes_manifest_info[chr]:
                chr_probe.append((chr, p))
        quot=len(chr_probe)/num_of_cores
        rem=len(chr_probe)%num_of_cores
        probe_number_per_process_l=[quot]*num_of_cores
        for i in range(rem):
            probe_number_per_process_l[i]+=1
        first_idx=0
        splitted_chr_probe=[]
        for num in probe_number_per_process_l:
            last_idx=first_idx+num
            chr_probe_list=chr_probe[first_idx:last_idx]
            chr_probe_dict={}
            for cp in chr_probe_list:
                c,p=cp
                if c not in chr_probe_dict:
                    chr_probe_dict[c]=[]
                chr_probe_dict[c].append(p)
            splitted_chr_probe.append(chr_probe_dict)
            first_idx=last_idx
        return splitted_chr_probe
    
    
    def get_ref_alt_by_reference_genome_match(self):
        probes_manifest_info=get_info_from_manifest(self.manifest)
        design_obj=InfiniumDesign(self.manifest)
        probe_design=design_obj.get_probes_infinium_design(probes_manifest_info)
        num_of_cores=self.core_num
        input_arg_list=[]
        split_chr_pos_l=RefAltMatch.split_chr_probe(probes_manifest_info,num_of_cores)
        for cpd in split_chr_pos_l:
            print(cpd)
            arg_l=[probes_manifest_info,self.ref_genome_path,probe_design,cpd]
            input_arg_list.append(arg_l)
        pool=Pool(num_of_cores)
        ref_alt_found_result=pool.map(ref_alt_match_tool_wrapper,input_arg_list)
        print(ref_alt_found_result)
        merged_probe_ref_alt={}
        for probe_ref_alt in ref_alt_found_result:
            for p in probe_ref_alt:
                merged_probe_ref_alt[p]=probe_ref_alt[p]
        pool.close()
        return merged_probe_ref_alt
        

