#ifndef FASTAFILE_H_
#define FASTAFILE_H_

#include <stdint.h>
#include <fstream>

using namespace std;


class FastaFile
{
public:
	ifstream fin;
	string seqName;
	string seq;
	string buffered_line;
	int buffered_outCur;
	uint64_t read_cur;
	std::streampos last;
	bool upperCase;
	
	inline FastaFile(string filename,bool _upperCase=true): upperCase(_upperCase), buffered_outCur(0), read_cur(0)
	{
		
		fin.open(filename.c_str());
		if(!fin.good())
		{
			fin.close();
			cerr<<"file "<<filename<<" cannot be opened"<<endl;
		}
	}
	
	inline bool readEntry()
	{
		seq="";
		seqName="";
		if(fin.eof())
			return false;
		
		fin>>seqName;

		if(seqName[0]!='>')
		{
			cerr<<"bad format!"<<endl;
			return false;
		}
		seqName=seqName.substr(1);
		last=fin.tellg();
		string line=" ";
		
		
		while(!fin.eof())
		{
			line="";
			fin>>line;
			//if(fin.eof())
			//	break;
			
			if(line[0]=='>')
			{
				fin.seekg(last);
				break;
				
			}
			
			if(upperCase)
				seq+=toUpperAndNoSpecialDNA(line);
			else
				seq+=line;
			last=fin.tellg();
		}
		this->read_cur=0;
		return seq!="";
	}

	inline char readNextNtFromLoadedSeq(){
		if(this->read_cur>=this->seq.length())
			return '\0';

		char ret=this->seq[this->read_cur];

		this->read_cur++;
		return ret;
	}
	

	/*void rewindLastLine(){
		fin.seekg(last);
	}*/
	inline char readNextNt()
	{
		//last=fin.tellg();

		if(this->buffered_outCur+1>=this->buffered_line.length()){

			if(fin.eof())
				return '\0';

			this->buffered_line="";
			fin>>this->buffered_line;

			if(this->buffered_line[0]=='>'){
				this->seqName=this->buffered_line.substr(1);
				this->buffered_line="";
				this->buffered_outCur=0;	
				return '>';
			}
			

			this->buffered_outCur=-1;		
		}

		this->buffered_outCur++;
		return this->buffered_line[this->buffered_outCur];
	}

	/*inline string readKmer(int k){
		string kmer;
		for(int i=0;i<k;i++){
			char c=this->readNextNt();
			switch(c){
				case '\0':
					return kmer;
				case '>':
					//have to rewind;
					this->rewindLastLine();
					return kmer;
				default:
					kmer+=c;
					break;
			}
		}

		return kmer;

	}*/

	inline ~FastaFile()
	{
		fin.close();
	}
	
};

#endif /*FASTAFILE_H_*/
