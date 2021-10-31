# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 16:44:32 2019

@author: zky
"""
from math import log
import sys
import numpy as np
#from scipy import stats
if len(sys.argv) <= 2:
    print("This script need two parameters. For example,\
          python normalize_peaks.py <input_bed> <ouput_bed>")
    sys.exit()
input_bed_file = sys.argv[1]
output_bed_file = sys.argv[2]

def MaxMinNormalization(x,Max,Min):
    if Max != Min :
        x = 1e-20 + (x - Min)*(1-1e-20) / (Max - Min)
    return x
#def Z_ScoreNormalization(x,mu,sigma):
#    x = (x - mu) / sigma;
#    return x
with open(input_bed_file) as peaks_bed:
    pvalue_array = []
    normalized_peaks = []
    max_pvalue = min_pvalue = 0
    for line in peaks_bed:
        data = line.replace('\n','').replace('\r','').split('\t')
        pvalue = float(data[4])
        pvalue_array.append(pvalue)
        normalized_peaks.append(data)
    if pvalue_array :
        max_pvalue = np.max(pvalue_array)
        min_pvalue = np.min(pvalue_array)
#    mu = np.average(pvalue_array)
#    sigma = np.std(pvalue_array)cd 
    for data in normalized_peaks:
        data[4] = MaxMinNormalization(float(data[4]),max_pvalue,min_pvalue)
        data[4] = -log(data[4],10)
with open(output_bed_file,'w') as output_file:
    for data in normalized_peaks:
        output_file.write('\t'.join(str(i) for i in data))
        output_file.write('\n')
