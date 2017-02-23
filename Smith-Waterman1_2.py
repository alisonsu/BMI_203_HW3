#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This code tests all of the substitution matrices at the gap opening and extension
penalties determined in Part 1 question 1, and plots a ROC curve for each matrix.
This answers Part 1, question 2.

@author: alisonsu
"""

import numpy as np
from hw2skeleton import *
import matplotlib.pyplot as plt
from sklearn import metrics

# Use gap open and extension penalties determined in part 1, q 1
gap_open = 8
gap_extend = 2

color_map = {0:'r',1:'b',2:'c',3:'m',4:'k'}   

fpr_dict = {}

# Read in substitution matrix files
filenames = ["BLOSUM50","BLOSUM62","MATIO","PAM100","PAM250"]
for ix,filename in enumerate(filenames):
    print(filename) # print statement lets me know code is still running
    score_matrix, key = align.read_score_matrix(filename)
    
    negpairs, pospairs = align.read_pos_and_neg() # read in negpairs and pospairs
    
    np_threshold = {}
    # Align negative pairs
    neg_score_list = np.empty(50)
    total_neg_pairs = 0
    for index,pair in enumerate(negpairs):
        total_neg_pairs += 1
        s1 = pair[0].sequence
        s2 = pair[1].sequence
        #print("Aligning:", pair)
        maxscoren=align.SWalign(s1,s2,gap_open, gap_extend,score_matrix, key)
        neg_score_list[index]=maxscoren
                   
    # Align positive pairs
    pos_score_list = np.empty(50)
    total_pos_pairs = 0
    for index,pair in enumerate(pospairs):
        total_pos_pairs +=1
        s1 = pair[0].sequence
        s2 = pair[1].sequence
        #print("Aligning:", pair)
        maxscorep=align.SWalign(s1,s2,gap_open, gap_extend,score_matrix,key)
        pos_score_list[index]=maxscorep
    
    pos_ID = np.ones(50) 
    neg_ID = np.zeros(50)
    
    ROC_ID = np.append(pos_ID, neg_ID)
    ROC_values = np.append(pos_score_list,neg_score_list)
    
    fpr, tpr, thresholds = metrics.roc_curve(ROC_ID, ROC_values,pos_label=1)
    
    #find closest value to tpr = 0.7
    idx = (np.abs(tpr-0.7)).argmin()
    fpr_at_tpr70 = fpr[idx]
    n = filename + " tpr:" + str(tpr[idx])
    # store these values in dictionary that can be printed so actual values can be seen
    fpr_dict[n] = float(fpr_at_tpr70)
    
    # Plot ROC curves
    lw = 2
    plt.plot(fpr, tpr, color=color_map[ix],lw=lw,label=filename)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.axhline(y=.7, xmin=0, xmax=1, linewidth=lw, color = 'navy')
    plt.legend(loc='lower right')
    plt.title("ROC curves for all scoring matrices")
plt.savefig("test")
print(fpr_dict)
