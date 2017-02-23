#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This code tests the effect of normalizing the alignment score to the length of 
the shorter sequence in the pair of sequences being aligned. This answers Part 1,
question 3.

@author: alisonsu
"""

import numpy as np
from hw2skeleton import *
import matplotlib.pyplot as plt
from sklearn import metrics

gap_open = 8
gap_extend = 2

color_map = {0:'r',1:'b',2:'c',3:'m',4:'k'}   

fpr_dict = {}

filenames = ["BLOSUM62"]
for ix,filename in enumerate(filenames):
    
    score_matrix, key = align.read_score_matrix(filename) # Read in score matrix file
    
    negpairs, pospairs = align.read_pos_and_neg() # read in negpairs and pospairs
    
    np_threshold = {}
    # Align negative pairs
    neg_score_list = np.empty(50)
    norm_neg_score_list = np.empty(50)
    total_neg_pairs = 0
    for index,pair in enumerate(negpairs):
        total_neg_pairs += 1
        s1 = pair[0].sequence
        s2 = pair[1].sequence
        maxscoren=align.SWalign(s1,s2,gap_open, gap_extend,score_matrix, key)
        norm_maxscoren = maxscoren/min(len(s1),len(s2)) # normalize to shorter seq len
        neg_score_list[index]=maxscoren
        norm_neg_score_list[index] = norm_maxscoren
                   
    # Align positive pairs
    pos_score_list = np.empty(50)
    norm_pos_score_list = np.empty(50)
    total_pos_pairs = 0
    for index,pair in enumerate(pospairs):
        total_pos_pairs +=1
        s1 = pair[0].sequence
        s2 = pair[1].sequence
        maxscorep=align.SWalign(s1,s2,gap_open, gap_extend,score_matrix,key)
        norm_maxscorep = maxscorep/min(len(s1),len(s2)) # normalize to shorter seq len
        pos_score_list[index]=maxscorep
        norm_pos_score_list[index] = norm_maxscorep
    
    # Format for ROC curve
    pos_ID = np.ones(50) 
    neg_ID = np.zeros(50)
    
    ROC_ID = np.append(pos_ID, neg_ID)
    ROC_values = np.append(pos_score_list,neg_score_list)
    ROC_values_norm = np.append(norm_pos_score_list,norm_neg_score_list)
    
    # Calculate ROC curves for both non-normalized and normalized
    fpr, tpr, thresholds = metrics.roc_curve(ROC_ID, ROC_values,pos_label=1)
    fpr_norm, tpr_norm, thresholds_norm = metrics.roc_curve(ROC_ID, ROC_values_norm,pos_label=1)
    
    # Plot ROC curves and save figure
    lw = 2
    plt.plot(fpr, tpr, color="r",lw=lw,label="BLOSUM62 not normalized")
    plt.plot(fpr_norm, tpr_norm, color="b",lw=lw,label="BLOSUM62 normalized")
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.axhline(y=.7, xmin=0, xmax=1, linewidth=lw, color = 'navy')
    plt.legend(loc='lower right')
    plt.title("ROC curves for BLOSUM62 with and without score normalization")
plt.savefig("normalized_BLOSUM62")