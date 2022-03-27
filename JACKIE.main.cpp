/***************************************************************************
 Copyright Wu Albert Cheng <albertwcheng@gmail.com>
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 *******************************************************************************/ 



#include <iostream>
#include <fstream>
#include <limits.h>
#include <vector>
#include <deque>
#include <map>
#include <set>
#include <list>
#include <queue>
#include <algorithm>
#include "Commonf.h"
#include "StringUtil.h"
#include "FastaFile.h"
#include "NucleicBit.h"
#include "KeyedPosition.h"





void printFGHelp(const char* programname)
{

	
	
	#ifdef __BIN_MAIN__
	cerr<<"Usage:"<<programname<<" ";
	cerr<<"fasta output_prefix output_suffix prefixLength outchrRef prefix2 ignoreLowercaseSeq kmerPatternMatch(e.g.,ATNWXXXXXNVG)"<<endl;
	cerr<<"\tGenerete binary files of simulated reads binned by prefix (only those with prefix2)"<<endl;
	#endif

	#ifdef __SORTTOBED_MAIN__
	cerr<<"Usage:"<<programname<<" ";
	cerr<<"chrRef binaryUnfold outBedfile thresholdLow thresholdHigh(=0 for no upperbound)"<<endl;
	cerr<<"\tOutput bed files for reads occuring [thresholdLow,thresholdHigh] inclusive times "<<endl;
   	#endif	
		
	#ifdef __JACKIE_MAP__
	cerr<<"Usage:"<<programname<<" ";
	cerr<<"fasta prefix2 ignoreLowercaseSeq(y/n) (singles/clusters) kmerPatternMatch(e.g.,ATNWXXXXXNVG)"<<endl;
	cerr<<"\tsort genome by k-mer sequences"<<endl;
	#endif


	#ifdef __JACKIE_VECTOR__
	cerr<<"Usage:"<<programname<<" ";
	cerr<<"fasta prefix2 ignoreLowercaseSeq(y/n) (singles/clusters) kmerPatternMatch(e.g.,ATNWXXXXXNVG)"<<endl;
	cerr<<"\tsort genome by k-mer sequences"<<endl;
	#endif
	
}

void foldGenomics_generateBinary2(int argc,const char** argv)
{
	if(argc<9)
	{
		printFGHelp(argv[0]);
		return;
	}
	
	GenomeNmersEncoder encoder(argv[3],argv[4],StringUtil::atoi(argv[5]),argv[6],23);
	encoder.transferFromFastaFile(argv[2],argv[7],argv[8][0]=='N' || argv[8][0]=='n' );
	
}

void foldGenomics_generateBinary3(int argc,const char** argv)
{
	if(argc<10)
	{
		printFGHelp(argv[0]);
		return;
	}
	
	GenomeNmersEncoder encoder(argv[3],argv[4],StringUtil::atoi(argv[5]),argv[6],StringUtil::atoi(argv[9]));
	encoder.transferFromFastaFile(argv[2],argv[7],argv[8][0]=='N' || argv[8][0]=='n' );
	
}

void foldGenomics_generateBinary4(int argc,const char** argv)
{
	if(argc<10)
	{
		printFGHelp(argv[0]);
		return;
	}
	
	GenomeNmersEncoder encoder(argv[3],argv[4],StringUtil::atoi(argv[5]),argv[6],argv[9]);
	encoder.transferFromFastaFile(argv[2],argv[7],argv[8][0]=='N' || argv[8][0]=='n' );
	
}

void foldGenomics_generateBinary4_new(int argc,const char** argv)
{
	if(argc<9)
	{
		printFGHelp(argv[0]);
		return;
	}
	
	GenomeNmersEncoder encoder(argv[2],argv[3],StringUtil::atoi(argv[4]),argv[5],argv[8]);
	encoder.transferFromFastaFile(argv[1],argv[6],argv[7][0]=='N' || argv[7][0]=='n' );
	
}


void JACKIE_mapSeq(int argc,const char** argv)
{
	if(argc<6)
	{
		printFGHelp(argv[0]);
		return;
	}
	
	GenomeNmersMap mapper(argv[5]);
	mapper.transferFromFastaFile(argv[1],argv[2],argv[3][0]=='N' || argv[3][0]=='n' , !strcmp(argv[4],"clusters"));
	
}

void JACKIE_vectorSeq(int argc,const char** argv)
{
	if(argc<6)
	{
		printFGHelp(argv[0]);
		return;
	}
	
	GenomeNmersVector vect(argv[5]);
	vect.transferFromFastaFile(argv[1],argv[2],argv[3][0]=='N' || argv[3][0]=='n' , !strcmp(argv[4],"clusters") );
	
}


void foldGenomics_generateBinary(int argc,const char** argv)
{
	if(argc<8)
	{
		printFGHelp(argv[0]);
		return;
	}

	GenomeNmersEncoder encoder(argv[3],argv[4],StringUtil::atoi(argv[5]),argv[6],StringUtil::atoi(argv[7]));
	encoder.transferFromFastaFile(argv[2]);

}


void foldGenomics_fold_stdsort(int argc,const char** argv)
{
	if(argc<6)
	{
		printFGHelp(argv[0]);
		return;
	}

	GenomeNmersDecoder decoder(argv[2]);
	int numEntries=decoder.numEntriesPending(KEYEDPOSITION);

	KeyedPosition* kps=new KeyedPosition[numEntries];

	ofstream ffileOut(argv[3],ios::binary|ios::out);

	int nEntries=0;
	int fEntries=0;

	int thr=StringUtil::atoi(argv[4]);
	
	
	int i=0;
	
	while(decoder.readEntry())
	{
		//decoder.kpos.printAsText(cerr);
		//cerr<<endl;
		kps[i++]=decoder.getKeyedPosition();

		nEntries++;

	}
	
	cerr<<nEntries<<" of sim reads read"<<endl;
	
	//cerr<<"*** Now fold ***"<<endl;
	//sort entries now

	std::sort(kps,kps+numEntries);

	bool outputPlainPosition=!strcmp(argv[5],"p");

	//now entries happening for a particular number of time and output

	KeyedPosition prePos=kps[0];
	int freq=1;
	
	for(i=1;i<numEntries;i++)
	{
		const KeyedPosition& thisPos=kps[i];

		if(thisPos==prePos)
		{
			freq++;
		}
		else
		{
			if(freq<=thr)
			{
				if(outputPlainPosition)
					ffileOut<<(Position)prePos;
				else
					ffileOut<<prePos;
				fEntries++;
			}

			prePos=thisPos;
			freq=1;
		}
	}

	if(freq<=thr)
	{
		if(outputPlainPosition)
			ffileOut<<(Position)prePos;
		else
			ffileOut<<prePos;
		fEntries++;
	}

	delete[] kps;


	cout<<nEntries<<" folded to "<<fEntries<<" with at most "<<thr<<" ocurrence(s)"<<endl;
	cerr<<fEntries<<" of sim reads after fold"<<endl;

	ffileOut.close();

}

void outKPtoBedFile(ofstream& ffileOut,SmartChrMap & indexedChrMap, int freq,vector<KeyedPosition> &store)
{

    string seq="";
    bool firstItem=true;
    for(vector<KeyedPosition>::iterator i=store.begin();i!=store.end();i++){
        
        ffileOut<<indexedChrMap.getChrMapInfoFromPChrID(i->chrID).k1; //chrom
        ffileOut<<"\t";


		if(firstItem){
			seq=Key3b2Nuc(i->b);
			firstItem=false;
		}

        int start0=i->getPos()-1;
        int end1=start0+seq.length();
        char strand=(i->isForward()?'+':'-');
        ffileOut<<start0<<"\t"<<end1<<"\t";

        ffileOut<<i->b<<"."<<freq<<"/"<<seq<<"\t";
        ffileOut<<freq<<"\t";
        ffileOut<<strand<<endl;
    }
    
}

void outKPtoBedFile_old(ofstream& ffileOut,SmartChrMap & indexedChrMap, int freq,vector<KeyedPosition> &store)
{


    for(vector<KeyedPosition>::iterator i=store.begin();i!=store.end();i++){
        
        ffileOut<<indexedChrMap.getChrMapInfoFromPChrID(i->chrID).k1; //chrom
        ffileOut<<"\t";
        int start0=i->getPos()-1;
        int end1=start0+20;
        char strand=(i->isForward()?'+':'-');
        ffileOut<<start0<<"\t"<<end1<<"\t";
        ffileOut<<i->b<<"."<<freq<<"\t";
        ffileOut<<freq<<"\t";
        ffileOut<<strand<<endl;
    }
    
}

void foldGenomics_foldToBed_stdsort_new(int argc,const char** argv)
{
    //-f3 chrRef binaryUnfold outBedfile thresholdLow thresholdHigh
    
	if(argc<6)
	{
		printFGHelp(argv[0]);
		return;
	}
    
    SmartChrMap indexedChrMap(argv[1]);
    
	GenomeNmersDecoder decoder(argv[2]);
    
	int numEntries=decoder.numEntriesPending(KEYEDPOSITION);
    
	KeyedPosition* kps=new KeyedPosition[numEntries];
    vector<KeyedPosition> store;
    
    
	ofstream ffileOut(argv[3],ios::out);
    
	int nEntries=0;
	int fEntries=0;
    
	int thrLo=StringUtil::atoi(argv[4]);
	int thrHi=StringUtil::atoi(argv[5]);


    
	
	int i=0;
	
	while(decoder.readEntry())
	{
		//decoder.kpos.printAsText(cerr);
		//cerr<<endl;
		kps[i++]=decoder.getKeyedPosition();
        
		nEntries++;
        
	}
	
	cerr<<nEntries<<" of sim reads read"<<endl;
	
	//cerr<<"*** Now fold ***"<<endl;
	//sort entries now
    
	std::sort(kps,kps+numEntries);
    
	//bool outputPlainPosition=!strcmp(argv[5],"p");
    
	//now entries happening for a particular number of time and output
    
	KeyedPosition prePos=kps[0];
    
    store.push_back(prePos);
	int freq=1;
	
	for(i=1;i<numEntries;i++)
	{
		const KeyedPosition& thisPos=kps[i];
        
		if(thisPos==prePos)
		{
            
			freq++;
		}
		else
		{
			if(freq>=thrLo && (thrHi==0 || freq<=thrHi))
			{
				
                outKPtoBedFile(ffileOut,indexedChrMap,freq,store);
                
				fEntries++;
			}
            
			prePos=thisPos;
			freq=1;
            store.clear();
            
		}
        
        store.push_back(thisPos);
	}
    
	if(freq>=thrLo && (thrHi==0 || freq<=thrHi))
	{
		outKPtoBedFile(ffileOut,indexedChrMap,freq,store);
        
		fEntries++;
	}
    
	delete[] kps;
    
    
	cout<<nEntries<<" folded to "<<fEntries<<" with ["<<thrLo<<","<<thrHi<<"] ocurrence(s)"<<endl;
	cerr<<fEntries<<" of sim reads after fold"<<endl;
    
	ffileOut.close();
    
}


void foldGenomics_foldToBed_stdsort(int argc,const char** argv)
{
    //-f3 chrRef binaryUnfold outBedfile thresholdLow thresholdHigh
    
	if(argc<7)
	{
		printFGHelp(argv[0]);
		return;
	}
    
    SmartChrMap indexedChrMap(argv[2]);
    
	GenomeNmersDecoder decoder(argv[3]);
    
	int numEntries=decoder.numEntriesPending(KEYEDPOSITION);
    
	KeyedPosition* kps=new KeyedPosition[numEntries];
    vector<KeyedPosition> store;
    
    
	ofstream ffileOut(argv[4],ios::out);
    
	int nEntries=0;
	int fEntries=0;
    
	int thrLo=StringUtil::atoi(argv[5]);
	int thrHi=StringUtil::atoi(argv[6]);


    
	
	int i=0;
	
	while(decoder.readEntry())
	{
		//decoder.kpos.printAsText(cerr);
		//cerr<<endl;
		kps[i++]=decoder.getKeyedPosition();
        
		nEntries++;
        
	}
	
	cerr<<nEntries<<" of sim reads read"<<endl;
	
	//cerr<<"*** Now fold ***"<<endl;
	//sort entries now
    
	std::sort(kps,kps+numEntries);
    
	//bool outputPlainPosition=!strcmp(argv[5],"p");
    
	//now entries happening for a particular number of time and output
    
	KeyedPosition prePos=kps[0];
    
    store.push_back(prePos);
	int freq=1;
	
	for(i=1;i<numEntries;i++)
	{
		const KeyedPosition& thisPos=kps[i];
        
		if(thisPos==prePos)
		{
            
			freq++;
		}
		else
		{
			if(freq>=thrLo && (thrHi==0 || freq<=thrHi))
			{
				
                outKPtoBedFile(ffileOut,indexedChrMap,freq,store);
                
				fEntries++;
			}
            
			prePos=thisPos;
			freq=1;
            store.clear();
            
		}
        
        store.push_back(thisPos);
	}
    
	if(freq>=thrLo && (thrHi==0 || freq<=thrHi))
	{
		outKPtoBedFile(ffileOut,indexedChrMap,freq,store);
        
		fEntries++;
	}
    
	delete[] kps;
    
    
	cout<<nEntries<<" folded to "<<fEntries<<" with ["<<thrLo<<","<<thrHi<<"] ocurrence(s)"<<endl;
	cerr<<fEntries<<" of sim reads after fold"<<endl;
    
	ffileOut.close();
    
}





void foldGenomics_fold(int argc,const char** argv)
{
	if(argc<5)
	{
		printFGHelp(argv[0]);
		return;
	}

	GenomeNmersDecoder decoder(argv[2]);
	
	map<KeyedPosition,uint> kpf;
	typedef map<KeyedPosition,uint>::iterator I;
	typedef pair<I,bool> IStat;
	typedef map<KeyedPosition,uint>::value_type V;
	
	ofstream ffileOut(argv[3],ios::binary|ios::out);

	int nEntries=0;
	int fEntries=0;

	unsigned int thr=StringUtil::atoi(argv[4]); //changed unsigned-signed

	while(decoder.readEntry())
	{
		//decoder.kpos.printAsText(cerr);
		//cerr<<endl;
		IStat stat=kpf.insert(V(decoder.kpos,1));
		if(!stat.second) //more than once!!
		{
			stat.first->second++;
		}

		nEntries++;

	}

	cerr<<nEntries<<" of sim reads read"<<endl;

	//cerr<<"*** Now fold ***"<<endl;
	for(I i=kpf.begin();i!=kpf.end();i++)
	{
		if(i->second<=thr)
		{
			fEntries++;

			ffileOut<<(i->first);

			//i->first.printAsText(cerr);
			//cerr<<endl;
		}
	}
	cout<<nEntries<<" folded to "<<fEntries<<" with at most "<<thr<<" ocurrence(s)"<<endl;
	cerr<<fEntries<<" of sim reads after fold"<<endl;
	
	ffileOut.close();

}

void foldGenomics_print(int argc,const char** argv,int formatBin)
{
	if(argc<3)
	{
		printFGHelp(argv[0]);
		return;
	}

	GenomeNmersDecoder decoder(argv[2]);



	int nEntries=0;

	
	while(decoder.readEntry(formatBin))
	{
		switch(formatBin)
		{

		case KEYEDPOSITION:
			decoder.getKeyedPosition().printAsText(cout);
			break;
		case COMPACTPOSITION:
			decoder.getCompactPosition().printAsText(cout);
			break;
		default:
			decoder.getPosition().printAsText(cout);
		}
		cout<<endl;

		nEntries++;

	}
	
	cerr<<nEntries<<" of entries read and printed"<<endl;

}


void foldGenomics_partitionChr(int argc,const char**argv)
{
	if(argc<7)
	{
		printFGHelp(argv[0]);
		return;
	}
	
	PositionChrPartitioner pcp(argv[2],argv[3],argv[4]);
	pcp.partition(argv[5],(!strcmp(argv[6],"k")?KEYEDPOSITION:POSITION));
}

void foldGenomics_sort(int argc,const char**argv)
{
	if(argc<3)
	{
		printFGHelp(argv[0]);
		return;
	}
	
	PositionSorter::sort(argv[2],argv[3]);
}

void foldGenomics_sortCompact(int argc,const char**argv)
{
	if(argc<3)
	{
		printFGHelp(argv[0]);
		return;
	}
	
	PositionSorter::sortCompact(argv[2],argv[3]);
}



int main(int argc, const char **argv)
{

	cerr<<"JACKIE (Jackie & Albert's CRISPR k-mer instance enumerator) v2.0"<<endl;
	cerr<<"[Built:"<<__DATE__<<" "<<__TIME__<<"]"<<endl;	
	#ifdef __BIN_MAIN__
	cerr<<"Subroutine JACKIE.bin"<<endl;
	foldGenomics_generateBinary4_new(argc,argv);
	#endif

	#ifdef __SORTTOBED_MAIN__
	cerr<<"Subroutine JACKIE.sortToBed"<<endl;
	foldGenomics_foldToBed_stdsort_new(argc,argv);
	#endif

	#ifdef __JACKIE_MAP__
	cerr<<"Subroutine JACKIE.map"<<endl;
	JACKIE_mapSeq(argc,argv);
	#endif

	#ifdef __JACKIE_VECTOR__
	cerr<<"Subroutine JACKIE.vector"<<endl;
	JACKIE_vectorSeq(argc,argv);
	#endif

	cerr<<"<Done>"<<endl;
	return 0;

}


