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


#ifndef KEYEDPOSITION_H_
#define KEYEDPOSITION_H_


#include "Nucleic.h"
#include "NucleicBit.h"
#include "FastaFile.h"

#define uint unsigned int

class KeyedPosition;

class Position
{
public:
	uint chrID;
	int pos;
	Position(int _chrID=-1,int _pos=INT_MIN) :chrID(_chrID), pos(_pos)
	{
		
	}
	inline int getPos() const
	{
		if(isForward())
			return pos;
		else
			return -1*pos;
	}
	inline bool isValid() const
	{
		return pos!=INT_MIN;
	}
	inline bool isForward() const
	{
		return pos>0;
	}
	inline string getStrand() const
	{
		return (isForward()?"+":"-");
	}

	inline bool operator > (const Position& right) const
	{
		if(chrID==right.chrID)
		{
			return getPos() > right.getPos();
		}
		else
			return chrID > right.chrID;
	}
	inline bool operator < (const Position& right) const
	{
		if(chrID==right.chrID)
		{
			return getPos() < right.getPos();
		}
		else
			return chrID < right.chrID;
	}
	
	inline void readFromKeyedPosition(istream& is);
	inline void printAsText(ostream& os) const
	{
		os<<chrID;
		os<<"\t";
		
		if(pos>0)
		{
		os<<"F:";
		os<<pos;
		}
		else
		{
		os<<"R:";
		os<<(-1*pos);
		}
	}
	
	
};

inline ostream& operator << (ostream& os, const Position& pos)
{
	os.write((char*)&pos.chrID,sizeof(uint));
	os.write((char*)&pos.pos,sizeof(int));
	return os;
}

inline istream& operator >> (istream& is, Position& pos)
{
	is.read((char*)&pos.chrID,sizeof(uint));
	
	if(is.eof())
	{
		pos.pos=INT_MIN;
	}

	is.read((char*)&pos.pos,sizeof(int));
	
	return is;
}

class CompactPosition: public Position
{
public:
	
	int end0;
	inline CompactPosition(uint _chrID=-1,int _start0=INT_MIN,int _end0=INT_MIN): Position(_chrID,_start0),end0(_end0)
	{
		
	}
	inline bool operator > (const Position& right) const
	{
		if(chrID==right.chrID)
		{
			return getPos() > right.getPos();
		}
		else
			return chrID > right.chrID;
	}
	inline bool operator < (const Position& right) const
	{
		if(chrID==right.chrID)
		{
			return getPos() < right.getPos();
		}
		else
			return chrID < right.chrID;
	}
	
	inline int getStart0() const
	{
		return this->getPos();
	}
	
	inline int getEnd0() const
	{
		if(isForward())
			return end0;
		else
			return -1*end0;
	}
	
	inline int getStart1() const
	{
		return getStart0()+1;
	}
	
	inline int getEnd1() const
	{
		return getEnd0()+1;
	}
	
	inline KeyPair<int,int> getBound11() const
	{
		return KeyPair<int,int>(getStart1(),getEnd1());
	}
	
	inline KeyPair<int,int> getBound01() const
	{
		return KeyPair<int,int>(getStart0(),getEnd1());
	}

	inline int getLength() const
	{
		return getEnd0()-getStart0()+1;
	}
	
	
	inline void printAsText(ostream& os) const
	{
		os<<chrID;
		os<<"\t";
		
		if(pos>0)
		{
		os<<"F:";
		os<<pos;
		}
		else
		{
		os<<"R:";
		os<<(-1*pos);
		}
		os<<"\t";
		if(end0>0)
		{
		os<<"F:";
		os<<end0;
		}
		else
		{
		os<<"R:";
		os<<(-1*end0);
		}		
		
	}
	
};

inline ostream& operator << (ostream& os, const CompactPosition& pos)
{
	os.write((char*)&pos.chrID,sizeof(uint));
	os.write((char*)&pos.pos,sizeof(int));
	os.write((char*)&pos.end0,sizeof(int));	
	return os;
}

inline istream& operator >> (istream& is, CompactPosition& pos)
{
	is.read((char*)&pos.chrID,sizeof(uint));
	
	if(is.eof())
	{
		pos.pos=INT_MIN;
	}

	is.read((char*)&pos.pos,sizeof(int));
	is.read((char*)&pos.end0,sizeof(int));	
	
	return is;
}


class KeyedPosition: public Position
{
public:
	
	NucKey3b a;
	NucKey3b b;
	
	KeyedPosition()
	{
		a=NULL_KEY;
		b=NULL_KEY;
	}


	inline void writeAsPosition(ostream& os) const
	{
		os<<(const Position&)(*this);
	}
	
	inline bool operator > (const KeyedPosition& right) const
	{
		if(b==right.b)
		{
			return a<right.a;
		}
		else
			return b<right.b;
	}


	inline bool absEquals(const KeyedPosition& right) const
	{
		return a==right.a && b==right.b && chrID==right.chrID && pos==right.pos;
	}
	inline bool operator < (const KeyedPosition& right) const
	{
		if(b==right.b)
		{
			return a>right.a;
		}
		else
			return b>right.b;
	}
	
	inline bool operator == (const KeyedPosition& right) const
	{
		return a==right.a && b==right.b;
	}
	
	inline bool operator !=(const KeyedPosition& right) const
	{
		return !(*this==right);
	}
	
	inline bool operator >=(const KeyedPosition& right) const
	{
		return !(*this<right);
		
	}
	inline bool operator <=(const KeyedPosition& right) const
	{
		return !(*this>right);
	}

	inline KeyedPosition(const char* seq,int length,uint _chrID,int _pos,bool forward):Position(_chrID,_pos*(forward?1:-1))
	{
		a=NULL_KEY;
		b=NULL_KEY;
		
		if(length>KEY_CAPACITY*2)
		{
			cerr<<"error try to encode nuc > "<<(KEY_CAPACITY*2)<<endl;
			return;
		}
		
		if(length>KEY_CAPACITY)
		{
			b=Nuc2Key3b(seq,KEY_CAPACITY);
			a=Nuc2Key3b(seq+KEY_CAPACITY,length-KEY_CAPACITY);
		}
		else
		{
			b=Nuc2Key3b(seq,length);
		}
	}
	inline void printAsText(ostream& os) const
	{
		bin(a,os);
		bin(b,os);
		os<<"\t";
		os<<chrID;
		os<<"\t";
		
		if(pos>0)
		{
		os<<"F:";
		os<<pos;
		}
		else
		{
		os<<"R:";
		os<<(-1*pos);
		}
	}
};

inline ostream& operator << (ostream& os, const KeyedPosition& kpos)
{
	os.write((char*)&kpos.a,sizeof(NucKey3b));
	os.write((char*)&kpos.b,sizeof(NucKey3b));
	os.write((char*)&kpos.chrID,sizeof(uint));
	os.write((char*)&kpos.pos,sizeof(int));
	return os;
}

inline istream& operator >> (istream& is,  KeyedPosition& kpos)

{
	is.read((char*)&kpos.a,sizeof(NucKey3b));
	if(is.eof())
	{
		kpos.pos=INT_MIN;
	}
	is.read((char*)&kpos.b,sizeof(NucKey3b));
	is.read((char*)&kpos.chrID,sizeof(uint));
	is.read((char*)&kpos.pos,sizeof(int));
	return is;
}


#define UNREAD 0
#define KEYEDPOSITION 1
#define POSITION 2
#define COMPACTPOSITION 3



class GenomeNmersDecoder
{
public:
	ifstream fin;
	
	int typeRead;
	
	KeyedPosition kpos;
	Position pos;
	CompactPosition cpos;
	
	inline GenomeNmersDecoder(string filename):fin(filename.c_str(),ios::in|ios::binary),typeRead(UNREAD)
	{
		if(!fin.good())
		{
			cerr<<"file "<<filename<<" cannot be opened"<<endl;
			return;
		}
	}
	
	inline int numEntriesPending(int typeToRead=KEYEDPOSITION)
	{
		std::streampos origpos=fin.tellg();
		fin.seekg(0,ios::end);
		std::streampos endpos=fin.tellg();
		std::streampos bytesPending=endpos-origpos;
		

		int sizeItem;
		switch(typeToRead)
		{
	
		case KEYEDPOSITION:
			sizeItem=sizeof(KeyedPosition);
			break;
		case COMPACTPOSITION:
			sizeItem=sizeof(CompactPosition);
			break;
		default:
			sizeItem=sizeof(Position);
			
		}

		
		fin.seekg(origpos,ios::beg);
		
		if(bytesPending%sizeItem!=0)
		{
			cerr<<"Error: the size doesn't match Format"<<endl;
		}
		
		return bytesPending/sizeItem;
	}
	
	
	inline bool readEntry(int typeToRead=KEYEDPOSITION)
	{
		if(fin.eof())
			return false;
		
		//cerr<<"format="<<typeToRead<<endl;
		
		switch(typeToRead)
		{
		case KEYEDPOSITION:
			fin>>kpos;
			typeRead=KEYEDPOSITION;
			return kpos.isValid();
		case COMPACTPOSITION:
			fin>>cpos;
			typeRead=COMPACTPOSITION;
			return cpos.isValid();
		default:
			//cerr<<"readPosition"<<endl;
			fin>>pos;	
			//pos.printAsText(cerr);
			typeRead=POSITION;
			return pos.isValid();
		}
	}
	
	inline const Position& getPosition()
	{
		switch(typeRead)
		{	
		case KEYEDPOSITION:
			return (const Position&)kpos;
		case COMPACTPOSITION:
			return (const Position&)cpos;
		default:
			return pos;
		}
	}
	
	inline const CompactPosition& getCompactPosition()
	{
		return cpos;
	}
	
	inline const KeyedPosition& getKeyedPosition()
	{
		return kpos;
	}
	
	inline ~GenomeNmersDecoder()
	{
		fin.close();
	}
};

inline void Position::readFromKeyedPosition(istream& is)
{
	KeyedPosition kpos;
	is>>kpos;
	(*this)=(Position)kpos;
}

class PositionSorter
{
public:
	inline static void sort(string filebinPos,string sortedOut)
	{
		GenomeNmersDecoder dec(filebinPos);
		
		int numEntries=dec.numEntriesPending(POSITION);
		
		Position *poss=new Position[numEntries];
		
		int ind=0;
		while(dec.readEntry(POSITION))
		{
			poss[ind++]=dec.getPosition();
		}
		
		std::sort(poss,poss+numEntries);
		
		ofstream fout(sortedOut.c_str(),ios::binary|ios::out);
		
		for(int i=0;i<numEntries;i++)
		{
			fout<<poss[i];
		}
		
		cout<<numEntries<<" entries sorted"<<endl;
		
		delete[] poss;
		fout.close();
		
	}

	inline static void sortCompact(string filebinPos,string sortedOut)
	{
		GenomeNmersDecoder dec(filebinPos);
		
		int numEntries=dec.numEntriesPending(POSITION);
		
		Position *poss=new Position[numEntries];
		
		int ind=0;
		while(dec.readEntry(POSITION))
		{
			poss[ind++]=dec.getPosition();
		}
		
		std::sort(poss,poss+numEntries);
		
		ofstream fout(sortedOut.c_str(),ios::binary|ios::out);
		
		CompactPosition cpp;
		
		int nCompact=0;
		
		for(int i=0;i<numEntries;i++)
		{
			Position& curPos=poss[i];
			
			//cerr<<"compare"<<endl;
			//cerr<<"\t";
			//cpp.printAsText(cerr);
			//cerr<<endl;
			//cerr<<"\t";
			//curPos.printAsText(cerr);
			//cerr<<endl;			
			
			if(cpp.chrID!=curPos.chrID || curPos.getPos()>cpp.getEnd0()+1)
			{
				if(cpp.isValid())
				{	
					fout<<cpp;
					nCompact++;
				}
				
				cpp.chrID=curPos.chrID;
				cpp.pos=curPos.pos;
				cpp.end0=curPos.pos;
			}else
				cpp.end0=curPos.pos;
			
		}
		
		if(cpp.isValid())
		{
			
			fout<<cpp;
			nCompact++;
		}
		
		cout<<numEntries<<" entries sorted to "<<nCompact<<" contiguous positions"<<endl;
		delete[] poss;
		fout.close();
		
	}
	
};


class SmartChrMap : public ChrMap
{
public:
	inline SmartChrMap(string filename): ChrMap(filename)
	{
		
	}
	
	typedef KeyPair<string,int> ChrMapInfo;
	ChrMapInfo getChrMapInfoFromPChrID(int id)
	{
		string chrName=this->_map[id];
		vector<string> chrPart;
		StringUtil::split(chrName,"|",chrPart);
		if(chrPart.size()>1)
		{
			vector<string> chrPart1;
			StringUtil::split(chrPart[1],":",chrPart1);
			return KeyPair<string,int>(chrPart[0],StringUtil::atoi(chrPart1[1]));
		}else
			return KeyPair<string,int>(chrPart[0],-1);
		
		
		//return KeyPair<string,int>(chrPart[0],(chrPart.size()>1?StringUtil::atoi(chrPart[1]):-1));
	}
};


class PositionChrPartitioner
{
public:
	
	ChrMap indexedChrMap;
	
	string outputPrefix;
	string outputSuffix;
	bool ignoreReverse;
	typedef buffered_ofstream<Position,vector<Position>,vector<Position>::iterator> BufOut;
	map<string,BufOut*> bfOuts; //binned by chr
	typedef map<string,BufOut*>::iterator I;
	typedef map<string,BufOut*>::value_type V;
	
	inline PositionChrPartitioner(string _indexedChrFile,string _outputPrefix,string _outputSuffix,bool _ignoreReverse=true): indexedChrMap(_indexedChrFile),outputPrefix(_outputPrefix),outputSuffix(_outputSuffix),ignoreReverse(_ignoreReverse)
	{
		/*for(map<int,string>::iterator i=indexedChrMap._map.begin();i!=indexedChrMap._map.end();i++)
		{
			cerr<<i->first<<"=>"<<i->second<<endl;
		}*/
	}
	
	inline void partition(string filebin,int formatBin=KEYEDPOSITION)
	{
		GenomeNmersDecoder dec(filebin);
		
		//cerr<<"filebin="<<filebin<<endl;
		//cerr<<"format="<<formatBin<<endl;
		
		int nPart=0;
		int nPass=0;
		
		while(dec.readEntry(formatBin))
		{
			nPass++;
			//cerr<<"a"<<nPass<<endl;
			const Position& pos=dec.getPosition();
			
			//if(nPass%100)
			//cerr<<"folding "<<nPass<<endl;
			
			if(ignoreReverse && !pos.isForward())
				continue;
			//cerr<<"b"<<nPass<<endl;
			nPart++;
			uint chr=pos.chrID;
			vector<string> chrPart;
			//cerr<<"c"<<nPass<<"chr="<<chr<<endl;
			string chrName=indexedChrMap._map[chr];
			//cerr<<"d"<<nPass<<"got "<<chrName<<endl;
			StringUtil::split(chrName,"|",chrPart);
			//cerr<<"e"<<nPass<<endl;
			string realName=chrPart[0];
			//cerr<<"f"<<nPass<<endl;
			BufOut* bfout=NULL;
			I i=bfOuts.find(realName);
			//cerr<<"g"<<nPass<<endl;
			if(i==bfOuts.end())
			{
				bfout=new BufOut(outputPrefix+realName+outputSuffix,false);
				bfOuts.insert(V(realName,bfout));
			}else
				bfout=i->second;
			//cerr<<"h"<<nPass<<endl;
			bfout->push(pos,ios::binary|ios::out|ios::app,false);
			//cerr<<"i"<<nPass<<endl;
			
			
		}
		
		cerr<<"finish loop"<<endl;
		
		for(I i=bfOuts.begin();i!=bfOuts.end();i++)
		{
			BufOut* bfout=i->second;
			bfout->flush(ios::binary|ios::out|ios::app,false);
			delete bfout;
		}
		
		cout<<nPass<<" entries passed among which "<<nPart<<" were partitioned"<<endl;
		
	}
	
	inline ~PositionChrPartitioner()
	{
		
	}
	
	
};


class GenomeNmersEncoder
{
public:

	ofstream fmap;
	int nmersize;
	int prefixLength;
	string outputPrefix;
	string outputSuffix;
	string nmerPattern;
	
	int keyStart;
	int keyEnd;

	typedef buffered_ofstream<KeyedPosition,vector<KeyedPosition>,vector<KeyedPosition>::iterator> BufOut;
	map<string,BufOut*> bfouts;

	
	inline ~GenomeNmersEncoder()
	{	
		for(map<string,BufOut*>::iterator i=bfouts.begin();i!=bfouts.end();i++)
		{
			BufOut* fout=(*i).second;
			fout->flush(ios::out|ios::app|ios::binary,false);
			delete fout;
		}
		fmap.close();
	}
	inline GenomeNmersEncoder(string _outputPrefix,string _outputSuffix,int _prefixLength,string filenameMap,int _nmersize):outputPrefix(_outputPrefix),outputSuffix(_outputSuffix),prefixLength(_prefixLength),fmap(filenameMap.c_str()),nmersize(_nmersize)
	{
		keyStart=0;
		keyEnd=nmersize;
	}

	inline GenomeNmersEncoder(string _outputPrefix,string _outputSuffix,int _prefixLength,string filenameMap,string _nmerPattern):outputPrefix(_outputPrefix),outputSuffix(_outputSuffix),prefixLength(_prefixLength),fmap(filenameMap.c_str()),nmerPattern(_nmerPattern)
	{
		nmersize=nmerPattern.length();

		for(int i=0;i<nmersize;i++){
			if(nmerPattern[i]=='X'){
				keyStart=i;
				break;
			}
		}

		for(int i=nmersize-1;i>=0;i--){
			if(nmerPattern[i]=='X'){
				keyEnd=i+1;
				break;
			}
		}
	}

	inline const char* getStrPtr(int i,const string& str)
	{
		return str.c_str()+i;
	}
	inline string getPrefix(int i,const string& str)
	{
		return str.substr(i,prefixLength);
	}

	inline bool acceptPAM(int i,const string& str)
    {
        
		const char* _tmpNmerPattern=this->nmerPattern.c_str();
        
        for(int j=i;j<i+nmersize;j++){
		
			char thisBase=str[j];

            if(thisBase>'T' || thisBase=='N'){

                return false;
            }
            
			char thisPatternBase=*_tmpNmerPattern;

			if(thisPatternBase!='X' && thisPatternBase!='N'){
				switch(thisPatternBase){
					case 'A':
					case 'C':
					case 'G':
					case 'T':
						if(thisBase!=thisPatternBase){
							return false;
						}
						break;
					case 'W':
						if(thisBase!='A' && thisBase!='T'){
							return false;
						}
						break;
					case 'S':
						if(thisBase!='C' && thisBase!='G'){
							return false;
						}
						break;
					case 'M':
						if(thisBase!='A' && thisBase!='C'){
							return false;
						}
						break;
					case 'K':
						if(thisBase!='G' && thisBase!='T'){
							return false;
						}
						break;							
					case 'R':
						if(thisBase!='A' && thisBase!='G'){
							return false;
						}
						break;
					case 'Y':
						if(thisBase!='C' && thisBase!='T'){
							return false;
						}
						break;
					case 'B':
						if(thisBase=='A'){
							return false;
						}
						break;
					case 'D':
						if(thisBase=='C'){
							return false;
						}
						break;
					case 'H':
						if(thisBase=='G'){
							return false;
						}
						break;
					case 'V':
						if(thisBase=='T'){
							return false;
						}
						break;
					default:
						return false;
						break;
				}

			}

			_tmpNmerPattern++;
        }

        return true;
    }


    inline bool acceptPAM_NGG(int i,const string& str)
    {
        //NGG = 22G 23G
        
        for(int j=i;j<i+nmersize;j++){
            //cerr<<"check "<<str[j]<<" " <<(str[j]>'T' || str[j]=='N')<<endl;
            if(str[j]>'T' || str[j]=='N'){
                //lowercase or N
                //cerr<<"******reject "<<str.substr(i,23)<<endl;
                return false;
            }
            
            /*f(str[j]=='n' || str[j]=='N'){
                //N
                //cerr<<"******reject "<<str.substr(i,23)<<endl;
                return false;
            }*/
            
            
        }
        
        
        return (str[i+nmersize-2]=='G' || str[i+nmersize-2]=='g') && (str[i+nmersize-1]=='G' || str[i+nmersize-1]=='g');
    }
    
	inline void transferFromFastaFile(string fastaFileName,string prefixConstraint="",bool autoUpperCase=true)
	{
		int curChrID=0;
		FastaFile ffile(fastaFileName,autoUpperCase);
		int nReads=0;
		
		int constraintLength=prefixConstraint.length();
		
		while(ffile.readEntry())
		{
			int nCurReads=0;
			curChrID++;
			
			fmap<<curChrID<<"\t"<<ffile.seqName<<endl;
			cout<<"Encoding "<<curChrID<<":"<<ffile.seqName<<" of length "<<ffile.seq.length()<<endl;
			
			int seqLength=ffile.seq.length();
			if(seqLength<nmersize)
			{
				cerr<<"Ignored: sequence of "<<ffile.seqName<<" has length smaller than nmersize"<<endl;
				cout<<": ignored: sequence has length smaller then nmersize "<<endl;
				continue;
			}
			
                
			for(int i=0;i<=seqLength-nmersize;i++)
			{
                
                
                if(!acceptPAM(i,ffile.seq)){
                    continue;
                }
                
                //cerr<<"accept "<<ffile.seq.substr(i,23)<<endl;
                
				BufOut* bout=NULL;
				string prefix=getPrefix(i+keyStart,ffile.seq);
				
				if(prefixConstraint!="" && prefix.substr(0,constraintLength)!=prefixConstraint)
				{
					continue;
				}
				
				map<string,BufOut*>::iterator bi=bfouts.find(prefix);
				if(bi==bfouts.end())
				{
					bout=new BufOut(outputPrefix+prefix+outputSuffix);
					bfouts.insert(map<string,BufOut*>::value_type(prefix,bout));
				}else
					bout=(*bi).second;
				
				KeyedPosition kpos(getStrPtr(i+keyStart,ffile.seq),keyEnd-keyStart,curChrID,i+keyStart+1,true); //nmersize-3 to discount the PAM sequence
                
				//kpos.printAsText(cerr);
				//cerr<<endl;
				
                bout->push(kpos,ios::out|ios::binary|ios::app,false);
				nCurReads++;
				
				
				
				
				
			}
			
			
			
			string rseq=reverse_complement(ffile.seq);
			
			
			for(int i=0;i<=seqLength-nmersize;i++)
			{

                if(!acceptPAM(i,rseq)){
                    continue;
                }
                
                //cerr<<"accept "<<rseq.substr(i,23)<<endl;
                
                
				BufOut* bout=NULL;
				string prefix=getPrefix(i+keyStart,rseq);
				
				if(prefixConstraint!="" && prefix.substr(0,constraintLength)!=prefixConstraint)
				{
					continue;
				}
				
				map<string,BufOut*>::iterator bi=bfouts.find(prefix);
				if(bi==bfouts.end())
				{
					bout=new BufOut(outputPrefix+prefix+outputSuffix);
					bfouts.insert(map<string,BufOut*>::value_type(prefix,bout));
				}else
					bout=(*bi).second;
                
				KeyedPosition kpos(getStrPtr(i+keyStart,rseq),keyEnd-keyStart,curChrID,seqLength-i-keyEnd+1/*seqLength-nmersize-i+1+3*/,false); //nmersize-3  seqLength-nmersize-i+1+3
				
                //kpos.printAsText(cerr);
				//cerr<<endl;
                
				bout->push(kpos,ios::out|ios::binary|ios::app,false);
				nCurReads++;
				
				
				
				
			}			
			
			cout<<" outputing "<<nCurReads<<" reads with prefixConstraint "<<prefixConstraint<<endl;
			nReads+=nCurReads;
			
		}
		
		cerr<<curChrID<<" sequence(s) processed."<<nReads<<" simulated read(s) outputed "<<endl;
		//cout<<curChrID<<" sequence(s) processed."<<nReads<<" simulated read(s) outputed "<<endl;
		
	}
	
	
};



class GenomeNmersMap
{
public:

	int nmersize;


	string nmerPattern;
	
	int keyStart;
	int keyEnd;
	map<NucKey3b/*NucKey*/,vector<Position> > NucKey_Position_Map;
	map<int/*chrID*/,string/*chromName*/> chrID2chrName_Map;
	
	inline ~GenomeNmersMap()
	{	

	}
	inline GenomeNmersMap(int _nmersize):nmersize(_nmersize)
	{
		keyStart=0;
		keyEnd=nmersize;
	}

	inline GenomeNmersMap(string _nmerPattern):nmerPattern(_nmerPattern)
	{
		nmersize=nmerPattern.length();

		for(int i=0;i<nmersize;i++){
			if(nmerPattern[i]=='X'){
				keyStart=i;
				break;
			}
		}

		for(int i=nmersize-1;i>=0;i--){
			if(nmerPattern[i]=='X'){
				keyEnd=i+1;
				break;
			}
		}
	}

	inline const char* getStrPtr(int i,const string& str)
	{
		return str.c_str()+i;
	}


	inline bool acceptPAM(int i,const string& str)
    {
        
		const char* _tmpNmerPattern=this->nmerPattern.c_str();
        
        for(int j=i;j<i+nmersize;j++){
		
			char thisBase=str[j];

            if(thisBase>'T' || thisBase=='N'){

                return false;
            }
            
			char thisPatternBase=*_tmpNmerPattern;

			if(thisPatternBase!='X' && thisPatternBase!='N'){
				switch(thisPatternBase){
					case 'A':
					case 'C':
					case 'G':
					case 'T':
						if(thisBase!=thisPatternBase){
							return false;
						}
						break;
					case 'W':
						if(thisBase!='A' && thisBase!='T'){
							return false;
						}
						break;
					case 'S':
						if(thisBase!='C' && thisBase!='G'){
							return false;
						}
						break;
					case 'M':
						if(thisBase!='A' && thisBase!='C'){
							return false;
						}
						break;
					case 'K':
						if(thisBase!='G' && thisBase!='T'){
							return false;
						}
						break;							
					case 'R':
						if(thisBase!='A' && thisBase!='G'){
							return false;
						}
						break;
					case 'Y':
						if(thisBase!='C' && thisBase!='T'){
							return false;
						}
						break;
					case 'B':
						if(thisBase=='A'){
							return false;
						}
						break;
					case 'D':
						if(thisBase=='C'){
							return false;
						}
						break;
					case 'H':
						if(thisBase=='G'){
							return false;
						}
						break;
					case 'V':
						if(thisBase=='T'){
							return false;
						}
						break;
					default:
						return false;
						break;
				}

			}

			_tmpNmerPattern++;
        }

        return true;
    }


    inline bool acceptPAM_NGG(int i,const string& str)
    {
        //NGG = 22G 23G
        
        for(int j=i;j<i+nmersize;j++){
            //cerr<<"check "<<str[j]<<" " <<(str[j]>'T' || str[j]=='N')<<endl;
            if(str[j]>'T' || str[j]=='N'){
                //lowercase or N
                //cerr<<"******reject "<<str.substr(i,23)<<endl;
                return false;
            }
            
            /*f(str[j]=='n' || str[j]=='N'){
                //N
                //cerr<<"******reject "<<str.substr(i,23)<<endl;
                return false;
            }*/
            
            
        }
        
        
        return (str[i+nmersize-2]=='G' || str[i+nmersize-2]=='g') && (str[i+nmersize-1]=='G' || str[i+nmersize-1]=='g');
    }
    
	inline void transferFromFastaFile(string fastaFileName,string prefixConstraint,bool autoUpperCase,bool collapseSameChromosomeToExtendedBed)
	{
		int curChrID=0;
		FastaFile ffile(fastaFileName,autoUpperCase);
		int nReads=0;
		
		int constraintLength=prefixConstraint.length();
		
		while(ffile.readEntry())
		{
			int nCurReads=0;
			curChrID++;
			
			
			cerr<<"Encoding "<<curChrID<<":"<<ffile.seqName<<" of length "<<ffile.seq.length()<<endl;
			chrID2chrName_Map.insert(map<int,string>::value_type(curChrID,ffile.seqName));
			
			int seqLength=ffile.seq.length();
			if(seqLength<nmersize)
			{
				cerr<<"Ignored: sequence of "<<ffile.seqName<<" has length smaller than nmersize"<<endl;
				cerr<<": ignored: sequence has length smaller then nmersize "<<endl;
				continue;
			}
			
                
			for(int i=0;i<=seqLength-nmersize;i++)
			{
                
                
                if(!acceptPAM(i,ffile.seq)){
                    continue;
                }
                
                //cerr<<"accept "<<ffile.seq.substr(i,23)<<endl;
                
				
				if(prefixConstraint!="" && ffile.seq.substr(i+keyStart,constraintLength)!=prefixConstraint)
				{
					continue;
				}
				
				
				KeyedPosition kpos(getStrPtr(i+keyStart,ffile.seq),keyEnd-keyStart,curChrID,i+keyStart+1,true); //nmersize-3 to discount the PAM sequence
                
				map<NucKey3b/*NucKey*/,vector<Position> >::iterator keyed_vector_I=NucKey_Position_Map.find(kpos.b);
				if(keyed_vector_I==NucKey_Position_Map.end()){
					NucKey_Position_Map.insert(map<NucKey3b/*NucKey*/,vector<Position> >::value_type(kpos.b,vector<Position>(1,(Position&)kpos)));
				}else{
					keyed_vector_I->second.push_back((Position&)kpos);
				}
				
				


				nCurReads++;
				
				
				
				
				
			}
			
			
			
			string rseq=reverse_complement(ffile.seq);
			
			
			for(int i=0;i<=seqLength-nmersize;i++)
			{

                if(!acceptPAM(i,rseq)){
                    continue;
                }
                
                //cerr<<"accept "<<rseq.substr(i,23)<<endl;
                

				
				if(prefixConstraint!="" && rseq.substr(i+keyStart,constraintLength)!=prefixConstraint)
				{
					continue;
				}
				
                
				KeyedPosition kpos(getStrPtr(i+keyStart,rseq),keyEnd-keyStart,curChrID,seqLength-i-keyEnd+1/*seqLength-nmersize-i+1+3*/,false); //nmersize-3  seqLength-nmersize-i+1+3
				
				map<NucKey3b/*NucKey*/,vector<Position> >::iterator keyed_vector_I=NucKey_Position_Map.find(kpos.b);
				if(keyed_vector_I==NucKey_Position_Map.end()){
					NucKey_Position_Map.insert(map<NucKey3b/*NucKey*/,vector<Position> >::value_type(kpos.b,vector<Position>(1,(Position&)kpos)));
				}else{
					keyed_vector_I->second.push_back((Position&)kpos);
				}
                

				nCurReads++;
				
				
				
				
			}			
			
			cerr<<" outputing "<<nCurReads<<" reads with prefixConstraint "<<prefixConstraint<<endl;
			nReads+=nCurReads;
			
		}


		for(map<NucKey3b/*NucKey*/,vector<Position> >::iterator i=NucKey_Position_Map.begin();i!=NucKey_Position_Map.end();i++){
			string seq=Key3b2Nuc(i->first);
			int seqLength=seq.length();
			if(collapseSameChromosomeToExtendedBed){

				bool sameChr=true;

				std:sort(i->second.begin(),i->second.end());
				vector<Position>&sortedPositions=i->second;
				int firstChrID=sortedPositions[0].chrID;
				string firstStrand=sortedPositions[0].getStrand();
				int firstStart0=sortedPositions[0].getPos()-1;

				string blockStarts="0";
				string blockLen=StringUtil::str(seqLength);
				string blockLengths=blockLen;

				for(int x=1;x<sortedPositions.size();x++){
					Position& thisPos=sortedPositions[x];
					if(thisPos.chrID!=firstChrID){
						sameChr=false;
						break;
					}
					if(thisPos.getStrand()!=firstStrand){
						firstStrand=".";
					}

					if(sortedPositions[x].getPos()-sortedPositions[x-1].getPos()<seqLength){
						sameChr=false; //illegal block, overlap!
						break;
					}
					blockStarts+=","+StringUtil::str(thisPos.getPos()-1-firstStart0);
					blockLengths+=","+blockLen;
				}

				if(sameChr){
					int lastEnd1=(seq.length()+sortedPositions[sortedPositions.size()-1].getPos()-1);
					cout<<chrID2chrName_Map[firstChrID]<<"\t"<<firstStart0<<"\t"<<lastEnd1
						<<"\t"<<i->first<<"."<<sortedPositions.size()<<"/"<<seq<<"/k="<<sortedPositions.size()<<"/s="<<(lastEnd1-firstStart0)
						<<"\t"<<sortedPositions.size()<<"\t"<<firstStrand<<"\t"<<firstStart0<<"\t"<<lastEnd1<<"\t0,0,0\t"<<sortedPositions.size()<<"\t"<<blockLengths<<"\t"<<blockStarts<<endl;
				}


			}else{
				for(vector<Position>::iterator j=i->second.begin();j!=i->second.end();j++){
					
					int start0=j->getPos()-1;
					int end1=start0+seq.length();
					char strand=(j->isForward()?'+':'-');
					int freq=i->second.size();
					cout<<chrID2chrName_Map[j->chrID]<<"\t"<<start0<<"\t"<<end1<<"\t"<<i->first<<"."<<freq<<"/"<<seq<<"\t"<<freq<<"\t"<<strand<<endl;
					
				}
			}
		}


		
		cerr<<curChrID<<" sequence(s) processed."<<nReads<<" simulated read(s) outputed "<<endl;
		//cout<<curChrID<<" sequence(s) processed."<<nReads<<" simulated read(s) outputed "<<endl;

		
		
	}
	
	
};


class GenomeNmersVector
{
public:

	int nmersize;


	string nmerPattern;
	
	int keyStart;
	int keyEnd;
	vector<KeyedPosition> kps;
	map<int/*chrID*/,string/*chromName*/> chrID2chrName_Map;
	
	inline ~GenomeNmersVector()
	{	

	}
	inline GenomeNmersVector(int _nmersize):nmersize(_nmersize)
	{
		keyStart=0;
		keyEnd=nmersize;
	}

	inline GenomeNmersVector(string _nmerPattern):nmerPattern(_nmerPattern)
	{
		nmersize=nmerPattern.length();

		for(int i=0;i<nmersize;i++){
			if(nmerPattern[i]=='X'){
				keyStart=i;
				break;
			}
		}

		for(int i=nmersize-1;i>=0;i--){
			if(nmerPattern[i]=='X'){
				keyEnd=i+1;
				break;
			}
		}
	}

	inline const char* getStrPtr(int i,const string& str)
	{
		return str.c_str()+i;
	}


	inline bool acceptPAM(int i,const string& str)
    {
        
		const char* _tmpNmerPattern=this->nmerPattern.c_str();
        
        for(int j=i;j<i+nmersize;j++){
		
			char thisBase=str[j];

            if(thisBase>'T' || thisBase=='N'){

                return false;
            }
            
			char thisPatternBase=*_tmpNmerPattern;

			if(thisPatternBase!='X' && thisPatternBase!='N'){
				switch(thisPatternBase){
					case 'A':
					case 'C':
					case 'G':
					case 'T':
						if(thisBase!=thisPatternBase){
							return false;
						}
						break;
					case 'W':
						if(thisBase!='A' && thisBase!='T'){
							return false;
						}
						break;
					case 'S':
						if(thisBase!='C' && thisBase!='G'){
							return false;
						}
						break;
					case 'M':
						if(thisBase!='A' && thisBase!='C'){
							return false;
						}
						break;
					case 'K':
						if(thisBase!='G' && thisBase!='T'){
							return false;
						}
						break;							
					case 'R':
						if(thisBase!='A' && thisBase!='G'){
							return false;
						}
						break;
					case 'Y':
						if(thisBase!='C' && thisBase!='T'){
							return false;
						}
						break;
					case 'B':
						if(thisBase=='A'){
							return false;
						}
						break;
					case 'D':
						if(thisBase=='C'){
							return false;
						}
						break;
					case 'H':
						if(thisBase=='G'){
							return false;
						}
						break;
					case 'V':
						if(thisBase=='T'){
							return false;
						}
						break;
					default:
						return false;
						break;
				}

			}

			_tmpNmerPattern++;
        }

        return true;
    }


    inline bool acceptPAM_NGG(int i,const string& str)
    {
        //NGG = 22G 23G
        
        for(int j=i;j<i+nmersize;j++){
            //cerr<<"check "<<str[j]<<" " <<(str[j]>'T' || str[j]=='N')<<endl;
            if(str[j]>'T' || str[j]=='N'){
                //lowercase or N
                //cerr<<"******reject "<<str.substr(i,23)<<endl;
                return false;
            }
            
            /*f(str[j]=='n' || str[j]=='N'){
                //N
                //cerr<<"******reject "<<str.substr(i,23)<<endl;
                return false;
            }*/
            
            
        }
        
        
        return (str[i+nmersize-2]=='G' || str[i+nmersize-2]=='g') && (str[i+nmersize-1]=='G' || str[i+nmersize-1]=='g');
    }
    
	void outputKP(int freq,vector<Position>& store,const KeyedPosition&preKP){
		string seq=Key3b2Nuc(preKP.b);
		int seqLength=seq.length();
		
		for(int i=0;i<store.size();i++){
	
			int start0=store[i].getPos()-1;
			int end1=start0+seqLength;
			char strand=(store[i].isForward()?'+':'-');
			
			cout<<chrID2chrName_Map[store[i].chrID]<<"\t"<<start0<<"\t"<<end1<<"\t"<<preKP.b<<"."<<freq<<"/"<<seq<<"\t"<<freq<<"\t"<<strand<<endl;
		}
	}

	void outputKPcluster(int freq,vector<Position>& store,const KeyedPosition&preKP){
		string seq=Key3b2Nuc(preKP.b);
		int seqLength=seq.length();
		bool sameChr=true;
		
		std:sort(store.begin(),store.end());
		vector<Position>&sortedPositions=store;
		int firstChrID=sortedPositions[0].chrID;
		string firstStrand=sortedPositions[0].getStrand();
		int firstStart0=sortedPositions[0].getPos()-1;

		string blockStarts="0";
		string blockLen=StringUtil::str(seqLength);
		string blockLengths=blockLen;

		for(int x=1;x<sortedPositions.size();x++){
			Position& thisPos=sortedPositions[x];
			if(thisPos.chrID!=firstChrID){
				sameChr=false;
				break;
			}
			if(thisPos.getStrand()!=firstStrand){
				firstStrand=".";
			}

			if(sortedPositions[x].getPos()-sortedPositions[x-1].getPos()<seqLength){
				sameChr=false; //illegal block, overlap!
				break;
			}
			blockStarts+=","+StringUtil::str(thisPos.getPos()-1-firstStart0);
			blockLengths+=","+blockLen;
		}

		if(sameChr){
			int lastEnd1=(seq.length()+sortedPositions[sortedPositions.size()-1].getPos()-1);
			cout<<chrID2chrName_Map[firstChrID]<<"\t"<<firstStart0<<"\t"<<lastEnd1
				<<"\t"<<preKP.b<<"."<<sortedPositions.size()<<"/"<<seq<<"/k="<<sortedPositions.size()<<"/s="<<(lastEnd1-firstStart0)
				<<"\t"<<sortedPositions.size()<<"\t"<<firstStrand<<"\t"<<firstStart0<<"\t"<<lastEnd1<<"\t0,0,0\t"<<sortedPositions.size()<<"\t"<<blockLengths<<"\t"<<blockStarts<<endl;
		}

	}

	inline void transferFromFastaFile(string fastaFileName,string prefixConstraint,bool autoUpperCase,bool collapseSameChromosomeToExtendedBed)
	{
		int curChrID=0;
		FastaFile ffile(fastaFileName,autoUpperCase);
		int nReads=0;
		
		int constraintLength=prefixConstraint.length();
		
		while(ffile.readEntry())
		{
			int nCurReads=0;
			curChrID++;
			
			
			cerr<<"Encoding "<<curChrID<<":"<<ffile.seqName<<" of length "<<ffile.seq.length()<<endl;
			chrID2chrName_Map.insert(map<int,string>::value_type(curChrID,ffile.seqName));
			
			int seqLength=ffile.seq.length();
			if(seqLength<nmersize)
			{
				cerr<<"Ignored: sequence of "<<ffile.seqName<<" has length smaller than nmersize"<<endl;
				cerr<<": ignored: sequence has length smaller then nmersize "<<endl;
				continue;
			}
			
                
			for(int i=0;i<=seqLength-nmersize;i++)
			{
                
                
                if(!acceptPAM(i,ffile.seq)){
                    continue;
                }
                
                //cerr<<"accept "<<ffile.seq.substr(i,23)<<endl;
                
				
				if(prefixConstraint!="" && ffile.seq.substr(i+keyStart,constraintLength)!=prefixConstraint)
				{
					continue;
				}
				
				
				KeyedPosition kpos(getStrPtr(i+keyStart,ffile.seq),keyEnd-keyStart,curChrID,i+keyStart+1,true); //nmersize-3 to discount the PAM sequence
                
				kps.push_back(kpos);
				
				


				nCurReads++;
				
				
				
				
				
			}
			
			
			
			string rseq=reverse_complement(ffile.seq);
			
			
			for(int i=0;i<=seqLength-nmersize;i++)
			{

                if(!acceptPAM(i,rseq)){
                    continue;
                }
                
                //cerr<<"accept "<<rseq.substr(i,23)<<endl;
                

				
				if(prefixConstraint!="" && rseq.substr(i+keyStart,constraintLength)!=prefixConstraint)
				{
					continue;
				}
				
                
				KeyedPosition kpos(getStrPtr(i+keyStart,rseq),keyEnd-keyStart,curChrID,seqLength-i-keyEnd+1/*seqLength-nmersize-i+1+3*/,false); //nmersize-3  seqLength-nmersize-i+1+3
				
				kps.push_back(kpos);
                

				nCurReads++;
				
				
				
				
			}			
			
			cerr<<" outputing "<<nCurReads<<" reads with prefixConstraint "<<prefixConstraint<<endl;
			nReads+=nCurReads;
			
		}

		std::sort(kps.begin(),kps.end());


		KeyedPosition prePos=kps[0];
		
		vector<Position> store;
		
		store.push_back(prePos);
		int freq=1;
		int fEntries=0;

		for(int i=1;i<kps.size();i++)
		{
			const KeyedPosition& thisPos=kps[i];
			
			if(thisPos==prePos)
			{
				
				freq++;
			}
			else
			{

				if(collapseSameChromosomeToExtendedBed){
					outputKPcluster(freq,store,prePos);
				}else{
					outputKP(freq,store,prePos);
				}	
				fEntries++;

				
				prePos=thisPos;
				freq=1;
				store.clear();
				
			}
			
			store.push_back(thisPos);
		}
		

		if(collapseSameChromosomeToExtendedBed){
			outputKPcluster(freq,store,prePos);
		}else{
			outputKP(freq,store,prePos);
		}	
			
		fEntries++;
		


		
		cerr<<curChrID<<" sequence(s) processed."<<nReads<<" simulated read(s) outputed "<<endl;
		//cout<<curChrID<<" sequence(s) processed."<<nReads<<" simulated read(s) outputed "<<endl;

		
		
	}
	
	
};

class GenomeNmersVector2
{
public:

	int nmersize;


	string nmerPattern;
	
	int keyStart;
	int keyEnd;
	vector<KeyedPosition> kps;
	map<int/*chrID*/,string/*chromName*/> chrID2chrName_Map;
	
	inline ~GenomeNmersVector2()
	{	

	}
	inline GenomeNmersVector2(int _nmersize):nmersize(_nmersize)
	{
		keyStart=0;
		keyEnd=nmersize;
	}

	inline GenomeNmersVector2(string _nmerPattern):nmerPattern(_nmerPattern)
	{
		nmersize=nmerPattern.length();

		for(int i=0;i<nmersize;i++){
			if(nmerPattern[i]=='X'){
				keyStart=i;
				break;
			}
		}

		for(int i=nmersize-1;i>=0;i--){
			if(nmerPattern[i]=='X'){
				keyEnd=i+1;
				break;
			}
		}
	}

	inline const char* getStrPtr(int i,const string& str)
	{
		return str.c_str()+i;
	}


	inline bool acceptPAM(int i,const string& str)
    {
        
		const char* _tmpNmerPattern=this->nmerPattern.c_str();
        
        for(int j=i;j<i+nmersize;j++){
		
			char thisBase=str[j];

            if(thisBase>'T' || thisBase=='N'){

                return false;
            }
            
			char thisPatternBase=*_tmpNmerPattern;

			if(thisPatternBase!='X' && thisPatternBase!='N'){
				switch(thisPatternBase){
					case 'A':
					case 'C':
					case 'G':
					case 'T':
						if(thisBase!=thisPatternBase){
							return false;
						}
						break;
					case 'W':
						if(thisBase!='A' && thisBase!='T'){
							return false;
						}
						break;
					case 'S':
						if(thisBase!='C' && thisBase!='G'){
							return false;
						}
						break;
					case 'M':
						if(thisBase!='A' && thisBase!='C'){
							return false;
						}
						break;
					case 'K':
						if(thisBase!='G' && thisBase!='T'){
							return false;
						}
						break;							
					case 'R':
						if(thisBase!='A' && thisBase!='G'){
							return false;
						}
						break;
					case 'Y':
						if(thisBase!='C' && thisBase!='T'){
							return false;
						}
						break;
					case 'B':
						if(thisBase=='A'){
							return false;
						}
						break;
					case 'D':
						if(thisBase=='C'){
							return false;
						}
						break;
					case 'H':
						if(thisBase=='G'){
							return false;
						}
						break;
					case 'V':
						if(thisBase=='T'){
							return false;
						}
						break;
					default:
						return false;
						break;
				}

			}

			_tmpNmerPattern++;
        }

        return true;
    }


    inline bool acceptPAM_NGG(int i,const string& str)
    {
        //NGG = 22G 23G
        
        for(int j=i;j<i+nmersize;j++){
            //cerr<<"check "<<str[j]<<" " <<(str[j]>'T' || str[j]=='N')<<endl;
            if(str[j]>'T' || str[j]=='N'){
                //lowercase or N
                //cerr<<"******reject "<<str.substr(i,23)<<endl;
                return false;
            }
            
            /*f(str[j]=='n' || str[j]=='N'){
                //N
                //cerr<<"******reject "<<str.substr(i,23)<<endl;
                return false;
            }*/
            
            
        }
        
        
        return (str[i+nmersize-2]=='G' || str[i+nmersize-2]=='g') && (str[i+nmersize-1]=='G' || str[i+nmersize-1]=='g');
    }
    
	void outputKP(int freq,vector<Position>& store,const KeyedPosition&preKP){
		string seq=Key3b2Nuc(preKP.b);
		int seqLength=seq.length();
		
		for(int i=0;i<store.size();i++){
	
			int start0=store[i].getPos()-1;
			int end1=start0+seqLength;
			char strand=(store[i].isForward()?'+':'-');
			
			cout<<chrID2chrName_Map[store[i].chrID]<<"\t"<<start0<<"\t"<<end1<<"\t"<<preKP.b<<"."<<freq<<"/"<<seq<<"\t"<<freq<<"\t"<<strand<<endl;
		}
	}

	void outputKPcluster(int freq,vector<Position>& store,const KeyedPosition&preKP){
		string seq=Key3b2Nuc(preKP.b);
		int seqLength=seq.length();
		bool sameChr=true;
		
		std:sort(store.begin(),store.end());
		vector<Position>&sortedPositions=store;
		int firstChrID=sortedPositions[0].chrID;
		string firstStrand=sortedPositions[0].getStrand();
		int firstStart0=sortedPositions[0].getPos()-1;

		string blockStarts="0";
		string blockLen=StringUtil::str(seqLength);
		string blockLengths=blockLen;

		for(int x=1;x<sortedPositions.size();x++){
			Position& thisPos=sortedPositions[x];
			if(thisPos.chrID!=firstChrID){
				sameChr=false;
				break;
			}
			if(thisPos.getStrand()!=firstStrand){
				firstStrand=".";
			}

			if(sortedPositions[x].getPos()-sortedPositions[x-1].getPos()<seqLength){
				sameChr=false; //illegal block, overlap!
				break;
			}
			blockStarts+=","+StringUtil::str(thisPos.getPos()-1-firstStart0);
			blockLengths+=","+blockLen;
		}

		if(sameChr){
			int lastEnd1=(seq.length()+sortedPositions[sortedPositions.size()-1].getPos()-1);
			cout<<chrID2chrName_Map[firstChrID]<<"\t"<<firstStart0<<"\t"<<lastEnd1
				<<"\t"<<preKP.b<<"."<<sortedPositions.size()<<"/"<<seq<<"/k="<<sortedPositions.size()<<"/s="<<(lastEnd1-firstStart0)
				<<"\t"<<sortedPositions.size()<<"\t"<<firstStrand<<"\t"<<firstStart0<<"\t"<<lastEnd1<<"\t0,0,0\t"<<sortedPositions.size()<<"\t"<<blockLengths<<"\t"<<blockStarts<<endl;
		}

	}

	inline void transferFromFastaFile(string fastaFileName,string prefixConstraint,bool autoUpperCase,bool collapseSameChromosomeToExtendedBed,int maxFreq)
	{
		int curChrID=0;
		FastaFile ffile(fastaFileName,autoUpperCase);
		int nReads=0;
		
		int constraintLength=prefixConstraint.length();
		
		while(ffile.readEntry())
		{
			int nCurReads=0;
			curChrID++;
			
			
			cerr<<"Encoding "<<curChrID<<":"<<ffile.seqName<<" of length "<<ffile.seq.length()<<endl;
			chrID2chrName_Map.insert(map<int,string>::value_type(curChrID,ffile.seqName));
			
			int seqLength=ffile.seq.length();
			if(seqLength<nmersize)
			{
				cerr<<"Ignored: sequence of "<<ffile.seqName<<" has length smaller than nmersize"<<endl;
				cerr<<": ignored: sequence has length smaller then nmersize "<<endl;
				continue;
			}
			
                
			for(int i=0;i<=seqLength-nmersize;i++)
			{
                
                
                if(!acceptPAM(i,ffile.seq)){
                    continue;
                }
                
                //cerr<<"accept "<<ffile.seq.substr(i,23)<<endl;
                
				
				if(prefixConstraint!="" && ffile.seq.substr(i+keyStart,constraintLength)!=prefixConstraint)
				{
					continue;
				}
				
				
				KeyedPosition kpos(getStrPtr(i+keyStart,ffile.seq),keyEnd-keyStart,curChrID,i+keyStart+1,true); //nmersize-3 to discount the PAM sequence
                
				kps.push_back(kpos);
				
				


				nCurReads++;
				
				
				
				
				
			}
			
			
			
			string rseq=reverse_complement(ffile.seq);
			
			
			for(int i=0;i<=seqLength-nmersize;i++)
			{

                if(!acceptPAM(i,rseq)){
                    continue;
                }
                
                //cerr<<"accept "<<rseq.substr(i,23)<<endl;
                

				
				if(prefixConstraint!="" && rseq.substr(i+keyStart,constraintLength)!=prefixConstraint)
				{
					continue;
				}
				
                
				KeyedPosition kpos(getStrPtr(i+keyStart,rseq),keyEnd-keyStart,curChrID,seqLength-i-keyEnd+1/*seqLength-nmersize-i+1+3*/,false); //nmersize-3  seqLength-nmersize-i+1+3
				
				kps.push_back(kpos);
                

				nCurReads++;
				
				
				
				
			}			
			
			cerr<<" outputing "<<nCurReads<<" reads with prefixConstraint "<<prefixConstraint<<endl;
			nReads+=nCurReads;
			
		}

		std::sort(kps.begin(),kps.end());


		KeyedPosition prePos=kps[0];
		
		vector<Position> store;
		
		store.push_back(prePos);
		int freq=1;
		int fEntries=0;

		for(int i=1;i<kps.size();i++)
		{
			const KeyedPosition& thisPos=kps[i];
			
			if(thisPos==prePos)
			{
				
				freq++;
			}
			else
			{
				if(freq<=maxFreq){
					if(collapseSameChromosomeToExtendedBed){
						outputKPcluster(freq,store,prePos);
					}else{
						outputKP(freq,store,prePos);
					}	
					fEntries++;
				}
				
				prePos=thisPos;
				freq=1;
				store.clear();
				
			}
			if(freq<=maxFreq){
				store.push_back(thisPos);
			}
		}
		
		if(freq<=maxFreq){
			if(collapseSameChromosomeToExtendedBed){
				outputKPcluster(freq,store,prePos);
			}else{
				outputKP(freq,store,prePos);
			}	
				
			fEntries++;
		}
		


		
		cerr<<curChrID<<" sequence(s) processed."<<nReads<<" simulated read(s) outputed "<<endl;
		//cout<<curChrID<<" sequence(s) processed."<<nReads<<" simulated read(s) outputed "<<endl;

		
		
	}
	
	
};

#endif /*KEYEDPOSITION_H_*/
