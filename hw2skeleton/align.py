import numpy as np
from hw2skeleton import *
import math

def read_score_matrix(score_matrix_name):
    """
    Reads in substitution matrix file and returns matrix and key
    """
    score_matrix_file = open(score_matrix_name, "r")
    sm = []
    # Save score matrix file as a list of lists
    templetter = []
    i = 0
    for index, line in enumerate(score_matrix_file):
        line = line.strip("\n")
        if line[0] == "#":
            continue
        elif i == 0: # Save letter key to list "key"
            line = line.replace(" ","")
            for letter in line:
                templetter.append(letter)
            key = templetter
            i+=1
        else:
            templine = []
            line = line.split()
            for number in line:
                templine.append(number)
            sm.append(templine)

    score_matrix = np.zeros([len(sm), len(sm[0])])
    for index,sub_list in enumerate(sm):
        score_matrix[index,:] = sm[index]
    return(score_matrix,key)


def read_pos_and_neg():
    """
    Reads in positive and negative sequence pairs, saves them as Sequence
    class, and returns lists of pairs
    """
    negpair_file = open("Negpairs.txt","r")
    negpairs = []
    for line in negpair_file:
        ID = line.split()
        ID1 = ID[0]
        ID2 = ID[1]
        s1 = read_FASTA(ID1)
        s2 = read_FASTA(ID2)
        pair = (s1,s2)
        negpairs.append(pair)
        
    pospair_file = open("Pospairs.txt","r")
    pospairs = []
    for line in pospair_file:
        ID = line.split()
        ID1 = ID[0]
        ID2 = ID[1]
        s1 = read_FASTA(ID1)
        s2 = read_FASTA(ID2)
        pair = (s1,s2)
        pospairs.append(pair)
    
    negpair_file.close()
    pospair_file.close()

    return(negpairs,pospairs)


def SWalign(seq1,seq2,gap_open,gap_extend,score_matrix,key):
    """
    Performs Smith-Waterman alignment with affine gap penalty
    """
    l1 = len(seq1)
    l2 = len(seq2)
    # Initialize matrices    
    M = np.zeros([l1+1, l2+1], dtype=np.int) # score matrix
    P = np.zeros([l1+1, l2+1], dtype=np.int) # Pointer matrix
    
    gap_penalty1 = gap_open
    gap_penalty2 = gap_open
    for i in range(1,l1+1):
        for j in range(1,l2+1):
            if P[i,j-1] == 1: # check if coming from a gap, in which case, use extension penalty
                gap_penalty1 = gap_extend
            elif P[i,j-1] == -1: # if P[i,j-1] came from vertical gap, do not allow horizontal gap
                gap_penalty1 = math.inf
            else: # otherwise, not coming from gap, so use gap_open penalty
                gap_penalty1 = gap_open
                
            if P[i-1,j] == -1:# check if coming from a gap, in which case, use extension penalty
                gap_penalty2 = gap_extend
            elif P[i-1,j] == 1: #if P[i-1,j] came from horizontal gap, do not allow vertical gap
                gap_penalty2 = math.inf
            else: # otherwise, not coming from gap, so use gap_open penalty
                gap_penalty2 = gap_open 
            # compute score of aligning to a gap for both sequences    
            gap_seq1 = M[i,j-1] - gap_penalty1
            gap_seq2 = M[i-1,j] - gap_penalty2
            
            # Look up location in "key" of letters being aligned
            loc1 = key.index(seq1[i-1].upper())
            loc2 = key.index(seq2[j-1].upper())
            # Determine score of those letters being aligned
            score = score_matrix[loc1, loc2]
            # Fill in main matrix
            M[i,j] = max(M[i-1,j-1]+score,gap_seq1,gap_seq2,0)

            # build pointer matrix
            if M[i,j] == gap_seq1:
                P[i,j] = 1 #gap in seq1
            elif M[i,j]== gap_seq2:
                P[i,j] = -1 #gap in seq2
            else:
                P[i,j]=0 #match
    # Find max score in score matrix = alignment score            
    mx = M.max()

    #Backtrace - note, this is not actually used in this assignment since, in this function,
    # I only return the max score. However, it is necessary for visualizing alignments and
    # troubleshooting, which I did, so I've left it here. Uncomment to use.
    """
    aseq1 = []
    aseq2 = []
    m,n = np.unravel_index(M.argmax(),M.shape)
    m = int(m)
    n = int(n)
    
    while M[m,n] > 0: # stop when hit 0 = end of local alignment
        # Use pointer matrix to backtrack through score matrix
        if P[m,n] == 1: # gap in seq1
            aseq1.append("-")
            aseq2.append(seq2[n-1])
            n -= 1
        elif P[m,n] == -1: # gap in seq2
            aseq1.append(seq1[m-1])
            aseq2.append("-")
            m -= 1
        else: # match
            aseq1.append(seq1[m-1])
            aseq2.append(seq2[n-1])
            m -=1
            n -=1
    # Currently have alignment in 2 lists, both reversed. So reverse and join into strings
    aseq1.reverse()
    final1 = "".join(aseq1)
    aseq2.reverse()
    final2 = "".join(aseq2)
    """
    
    return(mx)



def SWalign_return_static(seq1,seq2,gap_open,gap_extend,score_matrix,key):
    """
    This is the same code as SWalign, except that it returns the static alignments
    as well as the alignment score
    """
    l1 = len(seq1)
    l2 = len(seq2)
    
    M = np.zeros([l1+1, l2+1], dtype=np.int)
    P = np.zeros([l1+1, l2+1], dtype=np.int) # Pointer matrix
    
    gap_penalty1 = gap_open
    gap_penalty2 = gap_open
    for i in range(1,l1+1):
        for j in range(1,l2+1):
            if P[i,j-1] == 1:
                gap_penalty1 = gap_extend
            elif P[i,j-1] == -1: #if P[i,j-1] came from vertical gap, do not allow horizontal gap
                gap_penalty1 = math.inf
            else:
                gap_penalty1 = gap_open
            if P[i-1,j] == -1:
                gap_penalty2 = gap_extend
            elif P[i-1,j] == 1: #if P[i-1,j] came from horizontal gap, do not allow vertical gap
                gap_penalty2 = math.inf
            else:
                gap_penalty2 = gap_open
                
            gap_seq1 = M[i,j-1] - gap_penalty1
            gap_seq2 = M[i-1,j] - gap_penalty2
            
            # Look up location in "key" of letters being aligned
            loc1 = key.index(seq1[i-1].upper())
            loc2 = key.index(seq2[j-1].upper())
            # Determine score of those letters being aligned
            score = score_matrix[loc1, loc2]
            # Fill in main matrix
            M[i,j] = max(M[i-1,j-1]+score,gap_seq1,gap_seq2,0)
            
            if M[i,j] == gap_seq1:
                P[i,j] = 1 #deletion in seq1
            elif M[i,j]== gap_seq2:
                P[i,j] = -1 #insertion in seq1
            else:
                P[i,j]=0 #match
                
    mx = M.max()


    #Backtrace
    aseq1 = []
    aseq2 = []
    m,n = np.unravel_index(M.argmax(),M.shape)
    m = int(m)
    n = int(n)
    
    while M[m,n] > 0:
        if P[m,n] == 1: #gap in seq1
            aseq1.append("-")
            aseq2.append(seq2[n-1])
            n -= 1
        elif P[m,n] == -1: # gap in seq2
            aseq1.append(seq1[m-1])
            aseq2.append("-")
            m -= 1
        else:
            aseq1.append(seq1[m-1])
            aseq2.append(seq2[n-1])
            m -=1
            n -=1
    # Currently have alignment in 2 lists, both reversed. So reverse and join into strings.
    aseq1.reverse()
    final1 = "".join(aseq1)
    aseq2.reverse()
    final2 = "".join(aseq2)

    return(mx,final1,final2)

