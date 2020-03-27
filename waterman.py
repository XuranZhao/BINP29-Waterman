#!/usr/bin/env python3
import sys
import numpy
import string
import argparse
########## Parsing arguments ###########
desc="""This program which could do simple local alignment"""
parser=argparse.ArgumentParser(description=desc)
parser.add_argument('-seq1',metavar="Sequence1_File",help="Input sequence1 file",type=argparse.FileType('r'),required=1)
parser.add_argument('-seq2',metavar="Sequence2_File",help="Input sequence2 file",type=argparse.FileType('r'),required=1)
args=parser.parse_args()

def load_sequence(sequence_file):
    sequence_file = sequence_file.readlines()
    sequence = []; line_number = 0
    while line_number in range(len(sequence_file)):
        if sequence_file[line_number].startswith('>'):
            line_number +=1
        else:
            sequence.append(sequence_file[line_number].rstrip())
            line_number += 1
    del sequence_file
    return ''.join(sequence).upper()



## Score matrix for nucleotide alignment
NUC44=numpy.array([[5,-4,-4,-4,-2],\
                   [-4,5,-4,-4,-2],\
                   [-4,-4,5,-4,-2],\
                   [-4,-4,-4,5,-2],\
                   [-2,-2,-2,-2,-1]])

NBET='ATGCN'

## define the function for calculating score matrix and arrow matrix:
def scoreMat(NUC44,NBET,seq1,seq2,gap=-8):
    len1,len2=len(seq1),len(seq2)
    scorMat=numpy.zeros((len1+1,len2+1),int)
    arrowMat=numpy.zeros((len1+1,len2+1),int)

#    scorMat[0,:]=numpy.arange(len2+1)*gap
#    scorMat[:,0]=numpy.arange(len1+1)*gap
    arrowMat[0,:]=numpy.ones(len2+1)
    arrowMat[1:,0]=numpy.zeros(len1)
    for i in range(1,len1+1):
        for j in range(1,len2+1):
            s_mat=numpy.zeros(4)
            s_mat[0]=scorMat[i-1,j]+gap
            s_mat[1]=scorMat[i,j-1]+gap
            n1,n2=NBET.index(seq1[i-1]),NBET.index(seq2[j-1])
            s_mat[2]=scorMat[i-1,j-1]+NUC44[n1,n2]

            scorMat[i,j]=numpy.max(s_mat)
            arrowMat[i,j]=numpy.argmax(s_mat)
    return scorMat,arrowMat

## obtain the alignment of the sequences
def DynAlign(scorMat,arrow,seq1,seq2):
    aln_seq1,aln_seq2='',''
    flat_scorMat=numpy.ravel(scorMat)
    v,h=divmod(numpy.argmax(flat_scorMat),len(seq2)+1)
    # print (v,h)
    while True:
        if arrow[v,h]==0:
            aln_seq1+=seq1[v-1]
            aln_seq2+='-'
            v-=1
        elif arrow[v,h]==1:
            aln_seq1+='-'
            aln_seq2+=seq2[h-1]
            h-=1
        elif arrow[v,h]==2:
            aln_seq1+=seq1[v-1]
            aln_seq2+=seq2[h-1]
            v-=1
            h-=1
        elif arrow[v,h]==3:
            break
        if (v==0 and h==0) or scorMat[v,h]==0:
            break

    aln1=aln_seq1[::-1]
    aln2=aln_seq2[::-1]


    return aln1,aln2


seq1_file = args.seq1
seq2_file = args.seq2

seq1 = load_sequence(seq1_file)
seq2 = load_sequence(seq2_file)

scoreMatrix,arrowMatrix=scoreMat(NUC44,NBET,seq1,seq2,gap=-8)
aln1,aln2=DynAlign(scoreMatrix,arrowMatrix,seq1,seq2)

print ('\nThe aligned sequences are:\n')
print (aln1,"\n")
print (aln2)
