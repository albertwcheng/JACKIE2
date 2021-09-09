#!/bin/env python

from sys import *
from operator import itemgetter

def printUsageAndExit(programName):
    print >> stderr,programName,"filename spanMax > ofilename"
    exit(1)


def outMemLines(ostream,memLines,chroms,chromMin,chromMax,spanMax,exons):

    if(len(memLines)<1):
        return

    if len(chroms)>1:
        return
    
    thisChromSpan=chromMax-chromMin

    if spanMax>0 and thisChromSpan>spanMax:
        return

    exons.sort(key=itemgetter(0))

    memLineFields=memLines[0].split("\t")
    outputFields=[memLineFields[0],chromMin,chromMax,memLineFields[3]+"/k="+str(len(exons))+"/s="+str(chromMax-chromMin),memLineFields[4],'+',chromMin,chromMax,"0,0,0",len(exons)]
    
    blockSizes=[]
    blockStarts=[]
    for exon in exons:
        blockSizes.append(str(exon[1]-exon[0]))
        blockStarts.append(str(exon[0]-chromMin))

    outputFields.append(",".join(blockSizes))
    outputFields.append(",".join(blockStarts))

    print >> ostream,"\t".join([str(x) for x in outputFields])


if __name__=='__main__':
    programName=argv[0]
    args=argv[1:]
    try:
        filename,spanMax=args
    except:
        printUsageAndExit(programName)

    fil=open(filename)
    memLines=[]
    itemName=""
    chroms=set()
    chromMin=100000000000
    chromMax=0
    lino=0
    exons=[]
    for lin in fil:
        lino+=1
        if lino%1000000==1:
            print >> stderr,"processing line",lino
        lin=lin.rstrip("\r\n")
        fields=lin.split("\t")
        thisItemName=fields[3]
        thisChrom=fields[0]
        thisChromStart0=int(fields[1])
        thisChromEnd1=int(fields[2])

        if thisItemName!=itemName:
            outMemLines(stdout,memLines,chroms,chromMin,chromMax,spanMax,exons)
            memLines=[lin]
            itemName=thisItemName
            chroms=set()
            chroms.add(thisChrom)
            chromMin=thisChromStart0
            chromMax=thisChromEnd1
            exons=[[thisChromStart0,thisChromEnd1]]
        else:
            chroms.add(thisChrom)
            chromMin=min(chromMin,thisChromStart0)
            chromMax=max(chromMax,thisChromEnd1)
            memLines.append(lin)
            exons.append([thisChromStart0,thisChromEnd1])
            
    fil.close()

    outMemLines(stdout,memLines,chroms,chromMin,chromMax,spanMax,exons)

