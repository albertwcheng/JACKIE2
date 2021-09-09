from sys import *
from operator import itemgetter

def printUsageAndExit(programName):
    print >> stderr,programName,"filename minSpan maxSpan minNumBlocks maxNumBlocks nameFormer(col,sep,item)"
    exit(1)

def output(chrom,chromEntries,u,nameFormer):
    gStart0=chromEntries[u[0]][0]
    gEnd1=chromEntries[u[-1]][1]
    nBlocks=len(u)
    blockLens=[]
    blockStarts=[]
    names=[]
    for i in u:
        start0,end1,name,strand=chromEntries[i]
        blockLens.append(str(end1-start0))
        blockStarts.append(str(start0-gStart0))
        names.append(name+"/"+strand)
    
    fields=[chrom,str(gStart0),str(gEnd1),":".join(names),"0",".",str(gStart0),str(gEnd1),"0,0,0",str(nBlocks),",".join(blockLens),",".join(blockStarts)]
    print >> stdout,"\t".join(fields)


if __name__=='__main__':
    programName=argv[0]
    try:
        filename,minSpan,maxSpan,minNumBlocks,maxNumBlocks,nameFormer=argv[1:]
        minSpan=int(minSpan)
        maxSpan=int(maxSpan)
        minNumBlocks=int(minNumBlocks)
        maxNumBlocks=int(maxNumBlocks)
        nameFormer=nameFormer.split(",")
        nameFormer[0]=int(nameFormer[0])-1
        if len(nameFormer)==3:
            nameFormer[2]=int(nameFormer[2])-1
    except:
        printUsageAndExit(programName)

    bedItems=dict()

    fil=open(filename)
    for lin in fil:
        lin=lin.rstrip("\r\n")
        fields=lin.split("\t")
        chrom=fields[0]
        start0=int(fields[1])
        end1=int(fields[2])
        name=fields[nameFormer[0]]
        strand=fields[5]
        if len(nameFormer)==3:
            name=name.split(nameFormer[1])[nameFormer[2]]

        try:
            chromEntries=bedItems[chrom]
        except KeyError:
            chromEntries=[]
            bedItems[chrom]=chromEntries
        
        chromEntries.append((start0,end1,name,strand))

    fil.close()

    for chrom,chromEntries in bedItems.items():
        print >> stderr,"process",chrom
        chromEntries.sort(key=itemgetter(0))
        Q=[]
        for i in range(0,len(chromEntries)):
            Q.append([i])
        
        while len(Q)>0:
            u=Q.pop(0)
            #print >> stderr,u
            lastI=u[-1]
            lastStart,lastEnd=chromEntries[lastI][:2]

            hasChildren=False
            
            if len(u)<maxNumBlocks:
                for i in range(lastI+1,len(chromEntries)):
                    nextStart,nextEnd=chromEntries[i][:2]
                    #if nextStart<=lastEnd:
                    #    continue
                
                    if lastEnd+minSpan<=nextStart and lastEnd+maxSpan>=nextStart:
                        hasChildren=True
                        #print >> stderr,u+[i]
                        Q.append(u+[i])
                
            if not hasChildren:
                if len(u)>=minNumBlocks and len(u)<=maxNumBlocks:
                    output(chrom,chromEntries,u,nameFormer)

    