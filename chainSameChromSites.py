#!/usr/bin/env python2.7

from sys import *
from operator import itemgetter

def printUsageAndExit(programName):
    print >> stderr,programName,"filename spanMax majorityThreshold (range from 0-1, use -1 or 2 for print-all and print-only-one-chrom-exclusive) > ofilename"
    exit(1)


def outMemLines(ostream,chroms,spanMax,itemName,itemScore,majorityThres):

    if majorityThres>=1 and len(chroms)>1: #want only in one chromosome
        return
    
    for chromName,chromIntervals in chroms.items():
        chromIntervals.sort(key=itemgetter(0)) #sort by start coordinates
        chromMin=chromIntervals[0][0]
        chromMax=chromIntervals[-1][1]


        thisChromSpan=chromMax-chromMin

        if spanMax>0 and thisChromSpan>spanMax:
            return

        lastEnd=-1

        blockSizes=[]
        blockStarts=[]

        for thisStart,thisEnd in chromIntervals:
            if thisStart<lastEnd: #this has overlapped with last block, skip
                continue
            thisBlockSize=thisEnd-thisStart
            blockSizes.append(str(thisBlockSize))
            blockStarts.append(str(thisStart-chromMin))
            lastEnd=thisEnd
        
        fractionClusterOnThisChrom=float(len(blockStarts))/float(itemScore)


        if majorityThres>0 and fractionClusterOnThisChrom<majorityThres: #not majority, skip this chrom cluster
            continue



        #memLineFields=memLines[0].split("\t")
        outputFields=[chromName,chromMin,chromMax,itemName+"/k="+str(len(blockStarts))+"/s="+str(thisChromSpan)+"/f="+str(fractionClusterOnThisChrom),itemScore,'+',chromMin,chromMax,"0,0,0",len(blockStarts)]
        #print >> stderr,outputFields

        outputFields.append(",".join(blockSizes))
        outputFields.append(",".join(blockStarts))



        print >> ostream,"\t".join([str(x) for x in outputFields])


if __name__=='__main__':
    programName=argv[0]
    args=argv[1:]
    try:
        filename,spanMax,majorityThres=args
        spanMax=int(spanMax)
        majorityThres=float(majorityThres)
    except:
        printUsageAndExit(programName)

    fil=open(filename)
    memLines=[]
    itemName=""
    chroms=dict()
    chromMin=100000000000
    chromMax=0
    itemScore=0
    lino=0
    exons=[]
    for lin in fil:
        lino+=1
        if lino%1000000==1:
            print >> stderr,"processing line",lino
        lin=lin.rstrip("\r\n")
        fields=lin.split("\t")
        thisItemScore=fields[4]
        thisItemName=fields[3]
        thisChrom=fields[0].strip()
        thisChromStart0=int(fields[1])
        thisChromEnd1=int(fields[2])
        if(len(thisChrom)==0):
            continue

        if thisItemName!=itemName:
            outMemLines(stdout,chroms,spanMax,itemName,itemScore,majorityThres)
            itemName=thisItemName
            itemScore=thisItemScore
            chroms=dict()
            
        try:
            chromIntervals=chroms[thisChrom]
        except KeyError:
            chromIntervals=[]
            chroms[thisChrom]=chromIntervals
        
        chromIntervals.append([thisChromStart0,thisChromEnd1])
        #else:
            #chroms.add(thisChrom)
            #if len(chroms)==1: #no need to add if len(chrom)>1 to save time and memory.
            #    chromMin=min(chromMin,thisChromStart0)
            #    chromMax=max(chromMax,thisChromEnd1)
            #    memLines.append(lin)
            #    exons.append([thisChromStart0,thisChromEnd1])
            
    fil.close()

    outMemLines(stdout,chroms,spanMax,itemName,itemScore,majorityThres)

