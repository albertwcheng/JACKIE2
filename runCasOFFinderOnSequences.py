#!/usr/bin/env python

from sys import *
from tempfile import mkdtemp
import os

def printUsageAndExit(programName):
    print >> stderr,programName,"filename col1,sep,elementIdx mismatches genomeDir casIODir"
    exit(1)

def registerMatch(patternList,mismatches):
    while len(patternList)<=mismatches:
        patternList.append(0)
    
    #print >> stderr,patternList
    patternList[mismatches]+=1


if __name__=='__main__':
  
    try:
        programName=argv[0]
        args=argv[1:]
        filename,col1,mismatches,genomeDir,tmpDir=args
    except:
        printUsageAndExit(programName)

    if "," in col1:
        col,sep,elementIdx=col1.split(",")
        col=int(col)-1
        elementIdx=int(elementIdx)-1
    else:
        sep=""
        col=int(col1)-1
        
    #tmpDir=mkdtemp()
    print >> stderr,"casOffinder output dir:",tmpDir

    try:
        os.mkdir(tmpDir)
    except:
        print >> stderr,"cannot make directory",tmpDir

        
    CasOF_inname=tmpDir+"/in.txt"
    CasOF_outname=tmpDir+"/out.txt"
    
    fout=open(CasOF_inname,"w")
    print >> fout,genomeDir
    print >> fout,"NNNNNNNNNNNNNNNNNNNNNGG"

    fil=open(filename)

    for lin in fil:
        lin=lin.rstrip("\r\n")
        fields=lin.split("\t")
        if sep=="":
            seq=fields[col]
        else:
            seq=fields[col].split(sep)[elementIdx]
        print >> fout,seq+"NNN "+str(mismatches)

    fil.close()

    fout.close()

    
    
    os.system("cas-offinder "+CasOF_inname+" G "+CasOF_outname)
    
    
    fil=open(CasOF_outname)
    

    #fil=open("/private/var/folders/2b/jtscc7_x1tbdtjhfbmlsk9nm9f7b5j/T/tmpPPogVQ/out.txt")

    matchDict=dict()
    
    for lin in fil:
        lin=lin.rstrip("\r\n")
        #print >> stdout,lin
        pattern,chrom,location,seq,strand,thisMismatches=lin.split("\t")
        try:
            patternList=matchDict[pattern]
        except:
            patternList=[]
            matchDict[pattern]=patternList

        thisMismatches=int(thisMismatches)

        registerMatch(patternList,thisMismatches)

    fil.close()

    #for matchPattern,patternList in matchDict.items():
    #    print >> stdout,matchPattern[:-3]+"\t"+"/".join([str(x) for x in patternList])

    fil=open(filename)
    for lin in fil:
        lin=lin.rstrip("\r\n")
        fields=lin.split("\t")
        
        seq=fields[col]
        
    
        patternList=matchDict[seq.upper()+"NNN"]
        print >> stdout,lin+"\t"+"/".join([str(x) for x in patternList])

    
