#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is very similar to the code for optimizing the BLOSUM62 matrix, except
that it uses MATIO as the matrix to optimize. Still start with static alignments
generated using BLOSUM62, as specified in the HW.

@author: alisonsu
"""

import numpy as np
from hw2skeleton import *
import matplotlib.pyplot as plt
from sklearn import metrics
import random
import math

def static_align(s,pospairs,negpairs):
    
    #Rebuild score matrix
    score_matrix = np.zeros((24,24))
    inds = np.triu_indices_from(score_matrix)
    score_matrix[inds] = s
    score_matrix[(inds[1], inds[0])] = s
            
    #Calculate score in static alignment    
    pos_scores = []
    for pair in pospairs:
        se1 = pair[0]
        se2 = pair[1]
        score = 0
        gap_penalty1 = gap_open
        gap_penalty2 = gap_open
        for i in range(len(se1)):
            if se1[i] == "-":
                score -= gap_penalty1
                gap_penalty1 = gap_extend
        
            elif se2[i] == "-":
                score -= gap_penalty2
                gap_penalty2 = gap_extend
    
            else:
                loc1 = key.index(se1[i])
                loc2 = key.index(se2[i])
                score += score_matrix[loc1,loc2]
                gap_penalty1 = gap_open
                gap_penalty2 = gap_open
        pos_scores.append(score)
       
    neg_scores = []
    for pair in negpairs:
        se1 = pair[0]
        se2 = pair[1]
        score = 0
        gap_penalty1 = gap_open
        gap_penalty2 = gap_open
        for i in range(len(se1)):
            if se1[i] == "-":
                score -= gap_penalty1
                gap_penalty1 = gap_extend
        
            elif se2[i] == "-":
                score -= gap_penalty2
                gap_penalty2 = gap_extend
    
            else:
                loc1 = key.index(se1[i].upper())
                loc2 = key.index(se2[i].upper())
                score += score_matrix[loc1,loc2]
                gap_penalty1 = gap_open
                gap_penalty2 = gap_open
        neg_scores.append(score)
    
    return(pos_scores,neg_scores)

def objective_function(pos_scores,neg_scores):  
    # Generate ROC scores
    pos_ID = np.ones(50) 
    neg_ID = np.zeros(50)
    
    ROC_ID = np.append(pos_ID, neg_ID)
    ROC_values = np.append(pos_scores,neg_scores)
    
    fpr, tpr, thresholds = metrics.roc_curve(ROC_ID, ROC_values,pos_label=1)
    
    rev_fpr = np.fliplr([fpr])[0]
    rev_tpr = np.fliplr([tpr])[0]
        
    #find indices of values for fpr = 0.0, 0.1, 0.2, and 0.3
    id1 = (np.abs(rev_fpr-0.0)).argmin()
    id2 =(np.abs(rev_fpr-0.1)).argmin()
    id3 = (np.abs(rev_fpr-0.2)).argmin()
    id4 = (np.abs(rev_fpr-0.3)).argmin()
    
    # To optimize this sum using a minimization optimization, multiply by -1
    tpr_sum = (float(rev_tpr[id1])+float(rev_tpr[id2])+float(rev_tpr[id3])+float(rev_tpr[id4])) * -1
    print("value:", tpr_sum)
    return(tpr_sum)

"""
Begin code
"""
color_map = {0:'r',1:'b',2:'c',3:'m',4:'k'}   
lw = 2
gap_open = 8
gap_extend = 2

scorem = "BLOSUM62"
s_matrix, key = align.read_score_matrix(scorem)

negpairs, pospairs = align.read_pos_and_neg()

# First, generate static alignments using BLOSUM62
# Align negative pairs
neg_score_list = []
total_neg_pairs = 0
neg_pair_static = []
for pair in negpairs:
    total_neg_pairs += 1
    s1 = pair[0].sequence
    s2 = pair[1].sequence
    maxscoren,aseq1,aseq2=align.SWalign_return_static(s1,s2,gap_open,gap_extend,s_matrix,key)
    neg_score_list.append(maxscoren)
    apair = aseq1,aseq2
    neg_pair_static.append(apair)
               
# Align positive pairs
pos_score_list = []
pos_pair_static = []
total_pos_pairs = 0
for pair in pospairs:
    total_pos_pairs +=1
    s1 = pair[0].sequence
    s2 = pair[1].sequence
    maxscorep,aseq1,aseq2=align.SWalign_return_static(s1,s2,gap_open,gap_extend,s_matrix,key)
    pos_score_list.append(maxscorep)  
    apair = aseq1,aseq2
    pos_pair_static.append(apair)

#Calculate MATIO alignment from static alignment
scorem = "MATIO"
s_matrix, key = align.read_score_matrix(scorem) 
pos_scores = []
for pair in pos_pair_static:
    se1 = pair[0]
    se2 = pair[1]
    score = 0
    gap_penalty1 = gap_open
    gap_penalty2 = gap_open
    for i in range(len(se1)):
        if se1[i] == "-":
            score -= gap_penalty1
            gap_penalty1 = gap_extend
    
        elif se2[i] == "-":
            score -= gap_penalty2
            gap_penalty2 = gap_extend

        else:
            loc1 = key.index(se1[i])
            loc2 = key.index(se2[i])
            score += s_matrix[loc1,loc2]
            gap_penalty1 = gap_open
            gap_penalty2 = gap_open
    pos_scores.append(score)
   
neg_scores = []
for pair in neg_pair_static:
    se1 = pair[0]
    se2 = pair[1]
    score = 0
    gap_penalty1 = gap_open
    gap_penalty2 = gap_open
    for i in range(len(se1)):
        if se1[i] == "-":
            score -= gap_penalty1
            gap_penalty1 = gap_extend
    
        elif se2[i] == "-":
            score -= gap_penalty2
            gap_penalty2 = gap_extend

        else:
            loc1 = key.index(se1[i].upper())
            loc2 = key.index(se2[i].upper())
            score += s_matrix[loc1,loc2]
            gap_penalty1 = gap_open
            gap_penalty2 = gap_open
    neg_scores.append(score)
    
# Generate ROC scores
pos_ID = np.ones(50) 
neg_ID = np.zeros(50)

ROC_ID = np.append(pos_ID, neg_ID)
ROC_values = np.append(pos_scores,neg_scores)

fpr, tpr, thresholds = metrics.roc_curve(ROC_ID, ROC_values,pos_label=1)
plt.plot(fpr, tpr, color=color_map[0],lw=lw,label="MATIO")
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.0])
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.legend(loc='lower right')

# This is same optimization code as used in BLOSUM. For more detailed comments, see that code
# Vectorize top half of substitution_matrix so can be optimized over
s = s_matrix[np.triu_indices(24)]

posscore,negscore = static_align(s,pos_pair_static, neg_pair_static)
score = objective_function(posscore,negscore)
start_score = score

# Set paramaters for simulated annealing
T0 = 200 #know that maximum spread is 4, so if want to start by accepting ~98% of scores, this equates to ~200 starting temp
alpha = 0.8
print("Starting optimization")
for k in range(500):
    # Set temperature to follow cooling schedule
    T = T0 * math.pow(alpha,k)
    for iterations in range(50):
        new_s = np.empty(300)
        # Perturb substitution matrix using Guassian distributions centered around starting matrix value
        for i in range(len(new_s)):
            new_s[i] = np.random.normal(s[i],1)
        
        # Evaluate objective function    
        new_posscore,new_negscore = static_align(new_s,pos_pair_static, neg_pair_static)
        new_score = objective_function(new_posscore,new_negscore)
        #always accept new matrix if score is less than previous score
        if new_score < score:
            s = new_s
            score=new_score
        # if score > previous score, accept with probability based on cooling schedule
        else:
            P = math.exp(-abs(new_score-score)/T)
            print(P)
            if random.random()<P:
                s = new_s
                score=new_score
                print("Accepted new sub matrix from probability")
            else:
                print("Did not accept new sub matrix")
                
# Build final matrix and write to new file that can be read back in
final_score_matrix = np.zeros((24,24))
inds = np.triu_indices_from(final_score_matrix)
final_score_matrix[inds] = s
final_score_matrix[(inds[1], inds[0])] = s

outfile = open("MATIO_optimized","w")
outfile.write("# Optimized MATIO substitution matrix\n")
for letter in key:
    outfile.write(letter)
    outfile.write(" ")
outfile.write("\n")
for i in range(len(final_score_matrix)):
    for j in range(len(final_score_matrix)):
        outfile.write(str(final_score_matrix[i][j]))
        outfile.write(" ")
    outfile.write("\n")

outfile.close()

# Read back in optimized matrix
test_s_matrix, test_key = align.read_score_matrix("MATIO_optimized")

#Calculate score in static alignment using new matrix    
pos_scores = []
for pair in pos_pair_static:
    se1 = pair[0]
    se2 = pair[1]
    score = 0
    gap_penalty1 = gap_open
    gap_penalty2 = gap_open
    for i in range(len(se1)):
        if se1[i] == "-":
            score -= gap_penalty1
            gap_penalty1 = gap_extend
    
        elif se2[i] == "-":
            score -= gap_penalty2
            gap_penalty2 = gap_extend

        else:
            loc1 = test_key.index(se1[i])
            loc2 = test_key.index(se2[i])
            score += test_s_matrix[loc1,loc2]
            gap_penalty1 = gap_open
            gap_penalty2 = gap_open
    pos_scores.append(score)
   
neg_scores = []
for pair in neg_pair_static:
    se1 = pair[0]
    se2 = pair[1]
    score = 0
    gap_penalty1 = gap_open
    gap_penalty2 = gap_open
    for i in range(len(se1)):
        if se1[i] == "-":
            score -= gap_penalty1
            gap_penalty1 = gap_extend
    
        elif se2[i] == "-":
            score -= gap_penalty2
            gap_penalty2 = gap_extend

        else:
            loc1 = test_key.index(se1[i].upper())
            loc2 = test_key.index(se2[i].upper())
            score += test_s_matrix[loc1,loc2]
            gap_penalty1 = gap_open
            gap_penalty2 = gap_open
    neg_scores.append(score)

   
# Generate ROC scores and plot
pos_ID = np.ones(50) 
neg_ID = np.zeros(50)

ROC_ID = np.append(pos_ID, neg_ID)
ROC_values = np.append(pos_scores,neg_scores)

fpr_new, tpr_new, thresholds = metrics.roc_curve(ROC_ID, ROC_values,pos_label=1)
plt.plot(fpr_new, tpr_new, color=color_map[1],lw=lw,label="MATIO_optimized")
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.0])
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.legend(loc='lower right')

#Realign using new optimized matrix:
    
# Align negative pairs
neg_score_list = []
total_neg_pairs = 0
neg_pair_static = []
for pair in negpairs:
    total_neg_pairs += 1
    s1 = pair[0].sequence
    s2 = pair[1].sequence
    maxscoren,aseq1,aseq2=align.SWalign_return_static(s1,s2,gap_open,gap_extend,test_s_matrix,test_key)
    neg_score_list.append(maxscoren)
    apair = aseq1,aseq2
    neg_pair_static.append(apair)
               
# Align positive pairs
pos_score_list = []
pos_pair_static = []
total_pos_pairs = 0
for pair in pospairs:
    total_pos_pairs +=1
    s1 = pair[0].sequence
    s2 = pair[1].sequence
    maxscorep,aseq1,aseq2=align.SWalign_return_static(s1,s2,gap_open,gap_extend,test_s_matrix,test_key)
    pos_score_list.append(maxscorep)  
    apair = aseq1,aseq2
    pos_pair_static.append(apair)

pos_ID = np.ones(50) 
neg_ID = np.zeros(50)
    
ROC_ID = np.append(pos_ID, neg_ID)
ROC_values = np.append(pos_score_list,neg_score_list)
    
fpr, tpr, thresholds = metrics.roc_curve(ROC_ID, ROC_values,pos_label=1)

plt.plot(fpr, tpr, color=color_map[2],lw=lw,label="Realignment using MATIO_optimized")
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.0])
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.legend(loc='lower right')
plt.savefig("MATIO_optimization")
