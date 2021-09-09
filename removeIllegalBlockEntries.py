#!/bin/env python

from sys import *

def printUsageAndExit(programName):
    print >> stderr,programName,"filename passthru filtered"
    exit(1)

if __name__=='__main__':
    programName=argv[0]
    args=argv[1:]
    try:
        filename,passthru,filtered=args
    except:
        printUsageAndExit(programName)

    passthru=open(passthru,mode="w")
    filtered=open(filtered,mode="w")

    fil=open(filename)
    for lin in fil:
        lin=lin.rstrip("\r\n")
        fields=lin.split("\t")
        
        blockSizes=[int(x) for x in fields[10].split(",")]
        blockStarts=[int(x) for x in fields[11].split(",")]
        legal=True

        if(len(blockSizes)>1):
            for i in range(1,len(blockSizes)):
                if blockStarts[i-1]+blockSizes[i-1]>blockStarts[i]:
                    legal=False
                    break
         
        if legal:
            print >> passthru,lin
        else:
            print >> filtered,lin

    fil.close()
    passthru.close()
    filtered.close()
