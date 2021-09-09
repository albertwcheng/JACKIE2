from sys import *
from getopt import getopt

def printUsageAndExit(programName):
    print >> stderr,programName,"[options] filename > ofilename"
    exit(1)

if __name__=='__main__':
    programName=argv[0]
    #default values:
    gcRange=[0.4,0.6]
    noLowercase=True
    maxTandemT=4
    col0=3
    copyNumCol0=4
    sep="/"
    splitCom0=1
    startsWith=""
    showProgress=0
    cpRange=[1,-1]
    opts,args=getopt(argv[1:],'',['cp-range=','gc-range=','allow-lowercase','max-tandem-t=','field=','starts-with=','show-progress='])



    try:
        
        for o,v in opts:
            if o=='--gc-range':
                gcRange=[float(x) for x in v.split(",")]
            elif o=='--cp-range':
                cpRange=[int(x) for x in v.split(",")]
            elif o=='--allow-lowercase':
                noLowercase=False
            elif o=='--max-tandem-t':
                maxTandemT=int(v)
            elif o=='--field':
                col0,sep,splitCom0=v.split(",")
                col0=int(col0)-1
                splitCom0=int(splitCom0)-1
            elif o=='--starts-with':
                startsWith=v.upper()
            elif o=='--show-progress':
                showProgress=int(v)
        
        
        filename,=args
    except:
	#print(exc_info())        
	printUsageAndExit(programName)

    #print >> stderr,str(gcRange[0])+"<=gcRange<="+str(gcRange[1])
    lino=0
    fil=open(filename)
    for lin in fil:
        lino+=1
        if showProgress>0 and lino%showProgress==1:
            print >> stderr, "processing line",lino
        
        lin=lin.rstrip("\r\n")
        fields=lin.split("\t")
        copyNum=int(fields[copyNumCol0])
	itemName=fields[col0]
        itemNameSplits=itemName.split(sep)
        seq=itemNameSplits[splitCom0]
        if copyNum<cpRange[0]:
            continue
        if cpRange[1]>0 and copyNum>cpRange[1]:
            continue
        if noLowercase and seq.count("a")+seq.count("g")+seq.count("t")+seq.count("c")>0:
            continue

        seqUpper=seq.upper()

        if seqUpper.count("T")>maxTandemT:
            continue

        GCCount=seqUpper.count("G")+seqUpper.count("C")
        LenSeq=len(seqUpper)
        gcFrac=float(GCCount)/LenSeq

        #print >> stderr,seqUpper,gcFrac
        
        if len(startsWith)>0 and seqUpper[:len(startsWith)]!=startsWith:
            continue

        if gcFrac<gcRange[0]:
            continue

        if gcFrac>gcRange[1]:
            continue

        print >>  stdout,lin

    fil.close()
