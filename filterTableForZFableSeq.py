#!/bin/env python

from sys import *

ZFTripletString="GAAGCAGGAGTAGACGCCGGCGTCGAGGCGGGGGTGGATGCTGGTGTTAAAAACAAGAATACAACCACGACTAGAAGCAGGAGTATAATGATTTGATGGTAGCAACAGCATCCACCCCCGCCTCGACGCCGGCGTCTACTGCACCTT"

ZFTriplets=[]

for i in range(0,len(ZFTripletString),3):
    ZFTriplets.append(ZFTripletString[i:i+3])

#print(ZFTriplets,file=stderr)

def isSequenceZFable(seq):
    if len(seq)%3!=0: #not multiples of 3
        return False

    for i in range(0,len(seq),3):
        if seq[i:i+3] not in ZFTriplets:
            return False

    return True

def printUsageAndExit(programName):
    print(programName,"filename col1(,sep(,com1))",file=stderr)
    exit(1)

if __name__=='__main__':
    programName=argv[0]
    splitsep=""
    com=0
    try:
        filename,col1=argv[1:]
        col1=col1.split(",")
        col0=int(col1[0])-1
        if len(col1)>1:
            splitsep=col1[1]
        if len(col1)>2:
            com=int(col1[2])-1
    except:
        printUsageAndExit(programName)

    fil=open(filename)
    for lin in fil:
        lin=lin.rstrip("\r\n")
        fields=lin.split("\t")
        targetField=fields[col0]
        if splitsep!="":
            targetField=targetField.split(splitsep)[com]
        if(isSequenceZFable(targetField)):
            print(lin,file=stdout)
    fil.close()
    #print(seq,file=stdout,end=' ')
    #if isSequenceZFable(seq):
    #    print("is ZFable",file=stdout)
    #else:
    #    print("is NOT ZFable",file=stdout)
    
    
