#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This code uses BLOSUM50 to test a range of gap opening and gap extension
penalties, as asked in Part 1 question 1.

@author: alisonsu
"""

import numpy as np
from hw2skeleton import *
import matplotlib.pyplot as plt
from operator import itemgetter


# Read in score matrix file
score_matrix, key = align.read_score_matrix("BLOSUM50")

# REad in positive and negative pairs
negpairs, pospairs = align.read_pos_and_neg()

# Conduct alignments
np_threshold = {}
for i in range(1,21):
    print(i) # helps let me know that code is still running ok
    for j in range(1,6):
        # Align negative pairs
        neg_score_list = []
        total_neg_pairs = 0
        for pair in negpairs:
            total_neg_pairs += 1
            s1 = pair[0].sequence
            s2 = pair[1].sequence
            #print("Aligning:", pair)
            maxscoren=align.SWalign(s1,s2,i,j,score_matrix,key)
            neg_score_list.append(maxscoren)
                       
        # Align positive pairs
        pos_score_list = []
        total_pos_pairs = 0
        for pair in pospairs:
            total_pos_pairs +=1
            s1 = pair[0].sequence
            s2 = pair[1].sequence
            #print("Aligning:", pair)
            maxscorep=align.SWalign(s1,s2,i,j,score_matrix,key)
            pos_score_list.append(maxscorep)  
        
        # Since there are a total of 50 positive pairs, the value of the 15th lowest score
        # defines the threshold where 35 scores (or 70%) are higher than that value
        threshold = sorted(pos_score_list)[15]-1
        sorted_neg_pairs = sorted(neg_score_list)
        np_above_threshold = len([x for x in sorted_neg_pairs if x>threshold])/total_neg_pairs
        np_threshold[i,j] = np_above_threshold
sorted_np_threshold = sorted(np_threshold.items(),key=itemgetter(1),reverse=True)
print(sorted_np_threshold) # this gives actual values of fpr at tpr=0.7, since x axis of graph is hard to read
y_pos = np.arange(len(sorted_np_threshold))

# Construct bar plot of results.
x_val = []
y_val = []
for sublist in sorted_np_threshold:
    x_val.append(sublist[0])
    y_val.append(sublist[1])
    
plt.bar(y_pos,y_val,align='center')
plt.xticks(y_pos, x_val)
plt.ylabel('False positive rate')
plt.xlabel('(gap opening, gap extension)')
plt.savefig('Determining best false positive rate')