# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 18:53:30 2019

@author: zky
"""
from sys import argv
from math import log
from scipy import stats
input_bin25_file = argv[1]
ip_bin25_file = argv[2]
input_total_reads_count = int(argv[3])
ip_total_reads_count = int(argv[4])
peak_windows_number = int(argv[5])
output_ip_file = argv[6]
def windows_fisher_test(input_count, ip_count, input_total_reads_count, ip_total_reads_count):
    """fisher test for the PeakCalling of meyer"""
    site_input_rest_reads_count = input_total_reads_count - int(input_count)
    site_ip_rest_reads_count = ip_total_reads_count - int(ip_count)
    ip_oddsratio, ip_pvalue = stats.fisher_exact([[input_count, ip_count], [input_total_reads_count, ip_total_reads_count]], 'less')
    input_oddsratio, input_pvalue = stats.fisher_exact([[input_count, ip_count], [site_input_rest_reads_count, site_ip_rest_reads_count]], 'greater')
    return input_pvalue,ip_pvalue

def cluster_bin( bonferroni_filter_list ):
    bonferroni_peak = []
    peak_line = []
    idx = 0
    pre_end_position = 0
    for data in bonferroni_filter_list:
        distance = data[1] - pre_end_position
        if pre_end_position == 0 or distance > 0 :
            if peak_line :
                peak_region = peak_line[2] - peak_line[1]
                if peak_region >= 100 :
                    bonferroni_peak.append([])
                    bonferroni_peak[idx] = peak_line
                    idx += 1
            peak_line = []
            peak_line = data[:]
            pre_end_position = data[2]
        else:
            peak_line[2] = data[2]
            pre_end_position = data[2]
            peak_line.append(data[3])
    for data in bonferroni_peak:
        statistic, pval = stats.combine_pvalues(data[3:len(data)], method='fisher', weights=None)
        data[3] = pval
        del data[4:len(data)]
    return bonferroni_peak

with open (input_bin25_file) as input_bin25,open (ip_bin25_file) as ip_bin25:
    """Generate the list of bonferroni_filter_windows"""
    ip_bonferroni_filter_list = []
    ip_index = 0
    print ("Generate the list of bonferroni_filter_windows")
    while True:
        input_line = input_bin25.readline().rstrip("\n")
        ip_line = ip_bin25.readline().rstrip("\n")
        if input_line == '':
            break
        input_line_list = input_line.split("\t")
        ip_line_list = ip_line.split("\t")
        input_pvalue,ip_pvalue = windows_fisher_test(input_line_list[-1],ip_line_list[-1],input_total_reads_count,ip_total_reads_count)
        if (ip_pvalue < 0.05/peak_windows_number ):
            del ip_line_list[-1]
            ip_line_list.append(ip_pvalue)
            ip_line_list[1] = int(ip_line_list[1])
            ip_line_list[2] = int(ip_line_list[2])
            ip_bonferroni_filter_list.append([])
            ip_bonferroni_filter_list[ip_index] = ip_line_list
            ip_index += 1
"""Generate the list of bonferroni_filter_peaks"""
print ("Generate the list of bonferroni_filter_peaks")
ip_bonferroni_peak = cluster_bin(ip_bonferroni_filter_list[:])
"""Write the list of bonferroni_filter_peaks"""
print ("Write the list of bonferroni_filter_peaks")
with open(output_ip_file,'w') as output_file:
    for data in ip_bonferroni_peak:  
        output_file.write('\t'.join(str(i) for i in data))
        output_file.write('\n')