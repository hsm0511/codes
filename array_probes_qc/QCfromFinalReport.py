import os, sys, re, getopt, datetime
import subprocess
import pandas as pd
from IlluminaBeadArrayFiles import ClusterFile

class IntensityR(object):
    def __init__(self,Final_Report_path_handle,egt_handle):
        self.fr_path_name=Final_Report_path_handle
        self.egt=egt_handle
    
    @staticmethod
    def collect_R_from_FianlReports(FR_paths_file): # get R value from Final_Report_Plus.txt
        fr_list=[]
        for line in FR_paths_file:
            FR_path=line.rstrip()
            fr_list.append(FR_path)
        FR_paths_file.close()
        probe_R={}
        for fr_path in fr_list:
            fr=open(fr_path,'r')
            header_found=False
            for line in fr:
                l=line.rstrip().split('\t')
                if 'Allele1 - Plus' in line: # get column index
                    header_found=True
                    probe_name_idx=l.index('SNP Name')
                    R_idx=l.index('R')
                elif header_found==False:
                    continue
                else:   # get r value
                    probe_name=l[probe_name_idx]
                    if probe_name not in probe_R:
                        probe_R[probe_name]=[]
                    r_value=float(l[R_idx])
                    probe_R[probe_name].append(r_value)
            fr.close()
        return probe_R
    
    @staticmethod
    def get_quartile(r_value_list): # calculate quartile from r values
        r_value_ser=pd.Series(r_value_list)
        quartile_tmp=list(r_value_ser.quantile([0.25,0.5,0.75]))
        quartile_l=list(map(lambda x:float(x), quartile_tmp))
        return quartile_l
    
    def get_probes_R_info(self):
        egt=ClusterFile.read_cluster_file(self.egt)
        probe_R_info=IntensityR.collect_R_from_FianlReports(self.fr_path_name)
        probe_R_result={}
        for probe_name in probe_R_info:
            r_value_l=probe_R_info[probe_name]  # get list of samples' r values
            quartile_l=IntensityR.get_quartile(r_value_l)  # get quartile of samples' r values
            probe_R_result[probe_name]={'quartile':quartile_l}
            cluster_stat=IntensityR.get_R_cluster_stat(probe_name,egt)
            probe_R_result[probe_name]['cluster_stat']=cluster_stat
        return probe_R_result

    @staticmethod
    def get_R_cluster_stat(probe_name,egt_obj):
        clst_obj=egt_obj.get_record(probe_name)
        median_std_info={'AA':[],'AB':[],'BB':[]}
        # aa stat
        median_std_info['AA'].append(clst_obj.aa_cluster_stats.r_mean)
        median_std_info['AA'].append(clst_obj.aa_cluster_stats.r_dev)
        # ab stat
        median_std_info['AB'].append(clst_obj.ab_cluster_stats.r_mean)
        median_std_info['AB'].append(clst_obj.ab_cluster_stats.r_dev)
        # bb stat
        median_std_info['BB'].append(clst_obj.bb_cluster_stats.r_mean)
        median_std_info['BB'].append(clst_obj.bb_cluster_stats.r_dev)
        return median_std_info

class GenotypeCount(object):
    def __init__(self,final_reports_path):
        self.fr_path_name=final_reports_path
    
    def get_gt_count_from_FinalReport(self):  # function to count genotype of probes from Final_Report_Plus.txt
        probe_gt_ct={}
        for line in self.fr_path_name:
            if line.rstrip()!='':
                finalreport_path=line.rstrip()
                # read genotype from Final_Report_Plus.txt
                FR=open(finalreport_path,'r')
                header_found=False
                for line in FR:
                    l=line.rstrip().split('\t')
                    if 'Allele1 - Plus' in line:
                        header_found=True
                        # get index of columns
                        probe_name_idx=l.index('SNP Name')
                        allele_1_idx=l.index('Allele1 - Plus')
                        allele_2_idx=l.index('Allele2 - Plus')
                    elif header_found==False:
                        continue
                    else:
                        # count genotype for probes
                        probe_name=l[probe_name_idx]
                        gt=l[allele_1_idx]+'/'+l[allele_2_idx]
                        if probe_name not in probe_gt_ct:
                            probe_gt_ct[probe_name]={}
                        if gt not in probe_gt_ct[probe_name]:
                            probe_gt_ct[probe_name][gt]=0
                        probe_gt_ct[probe_name][gt]+=1
                FR.close()
        final_probe_gt_ct={}
        # merge hetero genotype
        for p in probe_gt_ct:
            final_probe_gt_ct[p]=GenotypeCount.hetero_merge(p,probe_gt_ct[p])
        return final_probe_gt_ct

    @staticmethod
    def weird_gt_chk(gt_list):  # check weird probes which have more than two bases from called genotypes
        weird_flag=False
        base_list=[]
        for gt in gt_list:
            base_list.append(gt[0])
            base_list.append(gt[2])
        # get unique base from 2000egt genotypes
        base_list=list(set(base_list))
        if '-' in base_list:
            base_list.remove('-')
        if len(base_list)>2:
            weird_flag=True
        return weird_flag

    @staticmethod
    def hetero_merge(probe_name,gt_ct_dict):
        merged_gt=gt_ct_dict.keys()
        # check weird gt form
        if GenotypeCount.weird_gt_chk(merged_gt)==True:
            print('werid genotype found in -> '+probe_name+' : '+','.join(merged_gt))
        merged_gt_ct_dict={}
        # merge hetero genotypes into one
        for gt in gt_ct_dict:
            if gt[0]!=gt[2] and gt[2]+'/'+gt[0] in merged_gt:
                merged_gt.remove(gt[2]+'/'+gt[0])
                break
        for gt in merged_gt:
            if gt[0]==gt[2]:
                merged_gt_ct_dict[gt]=gt_ct_dict[gt]
            else:
                another_gt=gt[2]+'/'+gt[0]
                if another_gt in gt_ct_dict.keys():
                    merged_gt_ct_dict[gt]=gt_ct_dict[gt]+gt_ct_dict[another_gt]
                else:
                    merged_gt_ct_dict[gt]=gt_ct_dict[gt]
        return merged_gt_ct_dict
    




