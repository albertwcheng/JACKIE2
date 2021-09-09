#!/bin/env python

from sys import *


def within(TestPoint,TestRange):
    return TestPoint>=TestRange[0] and TestPoint<=TestRange[1]

def printUsageAndExit(programName):
    print >> stderr,programName,"file1 file2 > overlappedFile"
    exit(1)

if __name__=='__main__':
    programName=argv[0]
    try:
        file1,file2=argv[1:]
    except:
        printUsageAndExit(programName)


    file1Dict=dict()
    
    fil=open(file1)
    for lin in fil:
        lin=lin.rstrip("\r\n")
        fields=lin.split("\t")
        chrom=fields[0]
        start0=int(fields[1])
        end1=int(fields[2])
        try:
            chromIntervals=file1Dict[chrom]
        except:
            chromIntervals=[]
            file1Dict[chrom]=chromIntervals
        
        chromIntervals.append((start0,end1,lin))

    fil.close()

    fil=open(file2)
    for lin in fil:
        lin=lin.rstrip("\r\n")
        fields=lin.split("\t")
        chrom=fields[0]
        start0=int(fields[1])
        end1=int(fields[2])
        try:
            chromIntervals=file1Dict[chrom]
        except KeyError:
            continue #not request, skip

        for chromInterval in chromIntervals:
            if within(start0,chromInterval) or within(end1,chromInterval):
                print >> stdout,chromInterval[2]+"\t"+lin

    fil.close()