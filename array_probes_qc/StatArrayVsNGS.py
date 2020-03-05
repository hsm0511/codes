import scipy.stats as stats
import numpy as np
import sys

class AlleleCount(object):
    @staticmethod
    def allele_count_from_pon(pon_info):    # allele counting for every base in genotypes
                                            # PON form example : A/A:13|A/B:98|B/B:3
        pon_info=pon_info.rstrip()
        allele_ct={}
        if pon_info!="" and '/' in pon_info:
            gt_count_l=pon_info.split('|')
            for gt_c in gt_count_l:
                gt,ct=gt_c.split(':')
                base_l=list(set(gt.split('/')))
                gt_l=gt.split('/')
                for b in base_l:
                    if b not in allele_ct:
                        allele_ct[b]=0
                    ct_in_gt=gt_l.count(b)
                    allele_ct[b]+=int(ct)*ct_in_gt
        return allele_ct
    
    @staticmethod
    def ref_alt_else_count(pon_egt,pon_ngs,ref,alt): # re-organize allele count info for ref, alt, else allele
        allele_info_dict={ref:[0,0],alt:[0,0],'else':['','']}
        pon_l=[pon_egt,pon_ngs]
        for i in range(len(pon_l)):
            p=pon_l[i]
            a_info=AlleleCount.allele_count_from_pon(p)
            if len(a_info)==0:  # no allele found from pon genotype
                for k in allele_info_dict:
                    allele_info_dict[k][i]='None'
            else:
                else_info=[]
                for base in a_info:
                    if base==ref:   # ref allele count
                        allele_info_dict[ref][i]=a_info[base]
                    elif base==alt: # alt allele count
                        allele_info_dict[alt][i]=a_info[base]
                    else:   # allele count for base which is neither ref nor alt
                        else_info.append(base+':'+str(a_info[base]))
                # merge else allele count info
                allele_info_dict['else'][i]='|'.join(else_info)
        return allele_info_dict

    @staticmethod
    def ref_alt_else_count_exception(ref,alt,indel_info):   # no Ref, Alt or indel probe
        allele_info_dict={ref:[0,0],alt:[0,0],'else':['','']}
        if indel_info=='indel':
            # indel probe, no allele counting
            ref_count=['indel']*2
            alt_count=['indel']*2
            else_count=['indel']*2
        else:
            # no ref alt info, no allele counting
            ref_count=['No ref alt']*2
            alt_count=['No ref alt']*2
            else_count=['No ref alt']*2
        allele_info_dict[ref]=ref_count
        allele_info_dict[alt]=alt_count
        allele_info_dict['else']=else_count
        return allele_info_dict
    
    def get_probe_ref_alt_count(self,probe_name,ref_alt_pon_record):
        info_l=ref_alt_pon_record[probe_name]
        ref=info_l[0]
        alt=info_l[1]
        pon_ngs=info_l[2]
        pon_egt=info_l[3]
        indel_info=info_l[4]
        if indel_info=='indel':
            allele_ct_info=AlleleCount.ref_alt_else_count_exception(ref,alt,indel_info)
        elif 'None' in [ref,alt] or '' in [ref,alt] or len(ref)!=1 or len(alt)!=1:
            allele_ct_info=AlleleCount.ref_alt_else_count_exception(ref,alt,indel_info)
        else:
            allele_ct_info=AlleleCount.ref_alt_else_count(pon_egt,pon_ngs,ref,alt)
        return allele_ct_info

class RefAltPON(object):

    @staticmethod
    def read_info_file(info_file_handle):
        probe_ref_alt_pon_indel={}
        # get pon info and count ref alt allele
        for line in info_file_handle:
            l=line.rstrip().split('\t')
            if line.rstrip()!='':
                if 'Name' in line:
                    # get index of columns
                    name_idx=l.index('Name')
                    pon_ngs_idx=l.index('PON_NGS')
                    pon_egt_idx=l.index('PON_Chip')
                    ref_idx=l.index('Ref')
                    alt_idx=l.index('Alt')
                    indel_idx=l.index('Indel_info')
                else:
                    # get ref, alt, pon, indel info
                    probe_name=l[name_idx]
                    pon_ngs=l[pon_ngs_idx]
                    pon_egt=l[pon_egt_idx]
                    ref=l[ref_idx]
                    alt=l[alt_idx]
                    if '/' in ref:
                        ref=ref.split('/')[0]
                    if '/' in alt:
                        alt=alt.split('/')[0]
                    indel_info=l[indel_idx]
                    info_l=[ref,alt,pon_ngs,pon_egt,indel_info]
                    probe_ref_alt_pon_indel[probe_name]=info_l
        return probe_ref_alt_pon_indel

class StatAlleleGTFreq(object):

    @staticmethod
    def check_stat_test_possible(allele_count_info):
        flag=True
        for k in allele_count_info:
            if k=='else':
                continue
            count_l=allele_count_info[k]
            if 'None' in count_l or 'indel' in count_l or 'No ref alt' in count_l:
                flag=False
                break
        return flag
    
    @staticmethod
    def check_chi2_possible(allele_count_info):
        flag_stat=StatAlleleGTFreq.check_stat_test_possible(allele_count_info)
        if flag_stat==False:
            return flag_stat
        flag_chi2=True
        test_d=[]
        for k in allele_count_info:
            if k=='else':
                continue
            else:
                # make 2X2 contingency table
                test_d.append(allele_count_info[k])
        test_d=np.array(test_d)
        if 0 in test_d.sum(axis=0) or 0 in test_d.sum(axis=1):
            flag_chi2=False
        return flag_stat and flag_chi2

    @staticmethod
    def bon_ferroni_adjust(p_val_l,alpha):
        qc_result=np.array(['pass']*len(p_val_l))
        adjust_p_val_arr=np.array(p_val_l)*len(p_val_l)
        adjust_p_val_arr[adjust_p_val_arr>1]=1
        adjust_p_val_l=list(adjust_p_val_arr)
        qc_tmp=adjust_p_val_arr<alpha
        qc_result[qc_tmp]='fail'
        qc_result=list(qc_result)
        return adjust_p_val_l,qc_result

    def probes_allele_freq_fisher_test(self,ref_alt_pon_record,alpha=0.05,p_adjust_method='Bon Ferroni'):
        probes_fisher_result={}
        probes_fisher_result_value_only=[]
        for probe in ref_alt_pon_record:
            ac=AlleleCount()
            allele_count=ac.get_probe_ref_alt_count(probe,ref_alt_pon_record)
            fisher_possible=StatAlleleGTFreq.check_stat_test_possible(allele_count)
            if fisher_possible:
                allele_ct_chip_ngs=[]
                for k in allele_count:
                    if k=='else':
                        continue
                    else:
                        # make 2X2 contingency table
                        allele_ct_chip_ngs.append(allele_count[k])
                # execute fisher test
                test_d=np.array(allele_ct_chip_ngs)
                print(test_d)
                p_value=stats.fisher_exact(test_d)[1]
                fisher_result=['',p_value]
                probes_fisher_result_value_only.append([probe,p_value])
            else:
                fisher_result=['No test','No test']
            probes_fisher_result[probe]=fisher_result
        # adjusting p_value
        if p_adjust_method=='Bon Ferroni':
            p_val_l=[pf[1] for pf in probes_fisher_result_value_only]
            p_adjusted_l,qc_result=StatAlleleGTFreq.bon_ferroni_adjust(p_val_l,alpha)
            for i in range(len(p_adjusted_l)):
                probe=probes_fisher_result_value_only[i][0]
                probes_fisher_result[probe]=[qc_result[i],p_adjusted_l[i]]
        return probes_fisher_result
    
    def probes_allele_freq_chi_square(self,ref_alt_pon_record,alpha=0.05,p_adjust_method='Bon Ferroni'):
        probes_chi2_result={}
        probes_chi2_result_value_only=[]
        for probe in ref_alt_pon_record:
            ac=AlleleCount()
            allele_count=ac.get_probe_ref_alt_count(probe,ref_alt_pon_record)
            chi2_possible=StatAlleleGTFreq.check_chi2_possible(allele_count)
            if chi2_possible:
                allele_ct_chip_ngs=[]
                for k in allele_count:
                    if k=='else':
                        continue
                    else:
                        # make 2X2 contingency table
                        allele_ct_chip_ngs.append(allele_count[k])
                # execute chi2 test
                try:
                    test_d=np.array(allele_ct_chip_ngs)
                    print(test_d)
                    p_value=stats.chi2_contingency(test_d,correction=True)[1]
                    chi2_result=['',p_value]
                    probes_chi2_result_value_only.append([probe,p_value])
                except ValueError as e:
                    print(probe+'\t'+'Error : '+str(e))
                    sys.exit(0)
            else:
                chi2_result=['No test','No test']
            probes_chi2_result[probe]=chi2_result
        # adjusting p_value
        if p_adjust_method=='Bon Ferroni':
            p_val_l=[pf[1] for pf in probes_chi2_result_value_only]
            p_adjusted_l,qc_result=StatAlleleGTFreq.bon_ferroni_adjust(p_val_l,alpha)
            for i in range(len(p_adjusted_l)):
                probe=probes_chi2_result_value_only[i][0]
                probes_chi2_result[probe]=[qc_result[i],p_adjusted_l[i]]
        return probes_chi2_result


    



