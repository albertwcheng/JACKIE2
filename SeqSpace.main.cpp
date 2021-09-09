#include <string>
#include <iostream>
#include <stdexcept>
#include <ctime>

#include "BitString.h"
#include "Commonf.h"
#include "FastaFile.h"
#include <queue>

using namespace std;

int printUsageAndExit(string programName)
{
    cerr<<programName<<" kmer fasta1 fasta2 ... fastaN"<<endl;
    return 1;
}

int printUsageAndExit_Mode3Offaware(string programName)
{
    cerr<<programName<<" kmer mthres mthress(?/?/?/?) seqPattern maxNumToOutput fasta1 fasta2 ... fastaN"<<endl;
    return 1;
}

int printUsageAndExit_Mode3Offaware2(string programName)
{
    cerr<<programName<<" kmer mthres mthress(?/?/?/?) seqPrefix maxNumToOutput fasta1 fasta2 ... fastaN"<<endl;
    return 1;
}

int printUsageAndExit_writer(string programName)
{
    cerr<<programName<<" outputSeqBits(.gz) kmer fasta1 fasta2 ... fastaN"<<endl;
    return 1;
}

int printUsageAndExit_reader(string programName)
{
    cerr<<programName<<" bitStringFile kmer mthreshold searchStringFile(.gz)"<<endl;
    return 1;
}


#define A 0
#define C 1
#define G 2
#define T 3

inline unsigned char NtToIdx(char c)
{
    switch(c){
            case 'A':
            return 0;
            case 'C':
            return 1;
            case 'G':
            return 2;
            case 'T':
            return 3;
        default:
            throw std::overflow_error(string("undefined nt ")+c);
    }
}

inline char IdxToNt(uint64_t idx)
{
    switch(idx){
        case 0:
            return 'A';
        case 1:
            return 'C';
        case 2:
            return 'G';
        case 3:
            return 'T';
        default:
            throw std::overflow_error(string("undefined idx ")+StringUtil::str(int(idx)));
    }
}

class SeqSpace
{
protected:
    
    int kmer;
    uint64_t jumpAheadBits;
public:
    BitString bits;

    
    SeqSpace(int _kmer,string inBitStringFileName, bool _useGzip=false):kmer(_kmer),bits(inBitStringFileName,(_useGzip?(uint64_t(1)<<(2*_kmer)):0)){

    }

    SeqSpace(int _kmer):kmer(_kmer),bits(uint64_t(1)<<(2*_kmer))
    {
        this->jumpAheadBits=~uint64_t(0);
        cerr<<"bitstring with "<<this->bits.numBits<<" created for k="<<kmer<<endl;
    }
    
    string indexToSeq64(uint64_t idx)
    {
        string seq;
        for(int i=0;i<32;i++){
            seq=IdxToNt(idx & 3)+seq;
            idx>>=2;
        }
        
        return seq;
        
    }

    uint64_t seqToIndex(string _seq){
        uint64_t idx=NtToIdx(_seq[0]);
        for(int i=1;i<this->kmer;i++)
        {
            idx<<=2;
            idx|=NtToIdx(_seq[i]);
        }
        
        return idx; 
    }

    static uint64_t seqToIndexStatic(string _seq){
        int _kmer=_seq.length();
        uint64_t idx=NtToIdx(_seq[0]);
        for(int i=1;i<_kmer;i++)
        {
            idx<<=2;
            idx|=NtToIdx(_seq[i]);
        }
        
        return idx; 
    }

    string indexToSeq(uint64_t idx)
    {
        string seq;
        for(int i=0;i<this->kmer;i++){
            seq=IdxToNt(idx & 3)+seq;
            idx>>=2;
        }
        
        return seq;
        
    }
    
    void registerSeq(string seq)
    {
        uint64_t idx=this->seqToIndex(seq);
        this->bits.setBit(idx,1);
    }
    
    void enumRegisterSeq(string _seq){
        uint64_t seqL=_seq.length();
        if(seqL<this->kmer){
            return;
        }
        for(uint64_t i=0;i<=seqL-this->kmer+1;i++)
        {
            /*if(i%1000000==1){
                cerr<<"processing pos "<<i<<endl;
            }*/
            string subSeq=_seq.substr(i,this->kmer);
            
            try{
                this->registerSeq(subSeq);
                this->registerSeq(reverse_complement(subSeq));
            }catch(std::exception e){
                
            }
        }
    }


    void enumRegisterSeqNGG(string _seq){
        uint64_t seqL=_seq.length();
        if(seqL<this->kmer){
            return;
        }
        for(uint64_t i=0;i<=seqL-this->kmer+1-3;i++)
        {
            /*if(i%1000000==1){
                cerr<<"processing pos "<<i<<endl;
            }*/
            string subSeq=_seq.substr(i,this->kmer);

            try{
                if(_seq.substr(i+this->kmer+1,2)=="GG"){
                    this->registerSeq(subSeq);
                }

                if(i>=3 && _seq.substr(i-3,2)=="CC"){
                    this->registerSeq(reverse_complement(subSeq));
                }
            }catch(std::exception e){
                
            }
        }
    }

    void printAbsentKmers(ostream &os){

        uint64_t output_count=0;

        for(uint64_t i=0;i<bits.numBits;i++){

            if(!bits.getBit(i)){
                os<<this->indexToSeq(i)<<"\n";
                output_count++;
            }

        }

        cerr<<output_count<<" VaKmers outputted"<<endl;
    }

    uint64_t getRandomUInt64()
    {

        uint64_t _randNum=(uint64_t)(rand()&255);
        for(int i=1;i<8;i++)
        {
            _randNum<<=8;
            _randNum|=(uint64_t)(rand()&255);
            
        }
        //cerr<<"got randNum "<<_randNum<<endl;
        return _randNum;

    }

    void sampleAbsentKmers(ostream &os,uint64_t samplingCount){

        uint64_t output_count=0;

        while(true){

            uint64_t randomNum=getRandomUInt64()%bits.numBits;

            if(!bits.getBit(randomNum)){
                os<<this->indexToSeq(randomNum)<<"\n";
                output_count++;

                //now set this bit to one to avoid double-printing:
                bits.setBit(randomNum,1);
            }

            if(output_count>=samplingCount)
                break; //got enough

        }

        cerr<<output_count<<" VaKmers outputted"<<endl;
    }


    void printAbsentKmers2(ostream &os){

        uint64_t output_count=0;

        Byte *_byte_ptr=this->bits._data;
        uint64_t _bit_idx=0;

        

        for(uint64_t _byte_idx=0;_byte_idx<bits.numBytes;_byte_idx++){
            Byte thisByte=*_byte_ptr;
            if(thisByte==255){
                //the whole byte is 11111111, skip, but advance idx
                _bit_idx+=8;
            }else{
                for(register int _bit_shift=0;_bit_shift<8;_bit_shift++)
                {
                    if(!(thisByte&(1<<_bit_shift)))
                    {
                        if(_bit_idx>=bits.numBits)
                            break;

                        os<<this->indexToSeq(_bit_idx)<<"\n";
                        output_count++;
                    }
                    _bit_idx++;
                }
                
            }
            

            //if(!bits.getBit(i)){
            //    os<<this->indexToSeq(i)<<endl;
            //}

            _byte_ptr++;
        }

        cerr<<output_count<<" VaKmers outputted"<<endl;
    } 


    void countAbsentKmers2(){

        uint64_t output_count=0;

        Byte *_byte_ptr=this->bits._data;
        uint64_t _bit_idx=0;

        

        for(uint64_t _byte_idx=0;_byte_idx<bits.numBytes;_byte_idx++){
            Byte thisByte=*_byte_ptr;
            if(thisByte==255){
                //the whole byte is 11111111, skip, but advance idx
                _bit_idx+=8;
            }else{
                for(register int _bit_shift=0;_bit_shift<8;_bit_shift++)
                {
                    if(!(thisByte&(1<<_bit_shift)))
                    {
                        if(_bit_idx>=bits.numBits)
                            break;

                        
                        output_count++;
                    }
                    _bit_idx++;
                }
                
            }
            

            //if(!bits.getBit(i)){
            //    os<<this->indexToSeq(i)<<endl;
            //}

            _byte_ptr++;
        }

        cerr<<output_count<<" VaKmers counted"<<endl;
    } 

    bool diff(SeqSpace& right)
    {
        bool same=true;
        if(this->bits.numBytes!=right.bits.numBytes){
            cout<<"numBytes different:"<<this->bits.numBytes<<" "<<right.bits.numBytes<<endl;
            return false;
        }


        for(uint64_t i=0;i<this->bits.numBytes;i++){
            if(this->bits.getBit(i)!=right.bits.getBit(i)){
                cout<<indexToSeq(i)<<"\t"<<(this->bits.getBit(i)?1:0)<<"\t"<<(right.bits.getBit(i)?1:0)<<endl;
                same=false;
            }    
        }


        return same;

    }

    ~SeqSpace(){
        
    }
};

class StepwiseSeqSpace:public SeqSpace
{   
    private:
    int idx_i;

    uint64_t fwd_idx;
    uint64_t rev_idx;
    uint64_t mask;
    uint64_t rev_A_setter;
    uint64_t rev_C_setter;
    uint64_t rev_G_setter;
    uint64_t rev_T_setter;
    
    public:

    void init_StepwiseSeqSpace(int _kmer){
        this->idx_i=0;
        this->fwd_idx=0;
        this->rev_idx=0;
  
        this->mask=3;

        
        this->rev_A_setter=0;
        this->rev_C_setter=1;
        this->rev_G_setter=2;
        this->rev_T_setter=3;

        for(int i=1;i<_kmer;i++)
        {
            this->mask<<=2;
            this->mask|=3;

            this->rev_A_setter<<=2;
            this->rev_C_setter<<=2;
            this->rev_G_setter<<=2;
            this->rev_T_setter<<=2;
        }
    }

    StepwiseSeqSpace(int _kmer):SeqSpace(_kmer)
    {

        init_StepwiseSeqSpace(_kmer);

    }

    StepwiseSeqSpace(int _kmer,string inBitStringFileName, bool _useGzip=false):SeqSpace(_kmer,inBitStringFileName,_useGzip){
        init_StepwiseSeqSpace(_kmer);
    }    

    

    string Uint64ToString(uint64_t idxer)
    {
        string s;
        for(int i=0;i<64;i++)
        {
            s=(((idxer&1)==1)?'1':'0')+s;
            idxer>>=1;
        }

        return s;
    }

    inline void see(char c)
    {
        switch(c)
        {
            case 'A':

                this->idx_i++; 
                if(this->idx_i>1){
                    this->fwd_idx<<=2;
                    this->rev_idx>>=2;
                }

                this->rev_idx|=this->rev_T_setter;

                break;
            case 'C':

                this->idx_i++;
                if(this->idx_i>1){
                    this->fwd_idx<<=2;
                    this->rev_idx>>=2;
                }

                this->fwd_idx|=1;
                this->rev_idx|=this->rev_G_setter;

                break;
            case 'G':

                this->idx_i++;
                if(this->idx_i>1){
                    this->fwd_idx<<=2;
                    this->rev_idx>>=2;
                }

                this->fwd_idx|=2;
                this->rev_idx|=this->rev_C_setter;

                break;
            case 'T':

                this->idx_i++;
                if(this->idx_i>1){
                    this->fwd_idx<<=2;
                    this->rev_idx>>=2;
                }


                this->fwd_idx|=3;
                this->rev_idx|=this->rev_A_setter;

                break;
            default:
                this->idx_i=0;
                this->fwd_idx=0;
                this->rev_idx=0;

                break;
        }

        if(this->idx_i==this->kmer)
        {

            //cerr<<"\tFwd "<<Uint64ToString(this->fwd_idx & this->mask)<<" ~ "<<this->indexToSeq(this->fwd_idx & this->mask)<<endl;
            //cerr<<"\tRev "<<Uint64ToString(this->rev_idx & this->mask)<<" ~ "<<this->indexToSeq(this->rev_idx & this->mask)<<endl;
            this->bits.setBit(this->fwd_idx & this->mask,1);
            this->bits.setBit(this->rev_idx & this->mask,1);
            this->idx_i--;
        }
    }

};





class SeqBitMask{
    public:
        uint64_t AND_MASK;
        uint64_t OR_MASK;
        int nmismatches;

        vector<uint64_t> positionMasks;
    
        uint64_t getPositionMask(int len,int pos){
            uint64_t mask;
            if(pos==0){
                mask=3;
            }else{
                mask=0;
            }
            for(int i=1;i<len;i++){
                mask<<=2;
                if(pos==i){
                    mask|=3;
                }
            }

            return mask;
        }

        SeqBitMask(string _seq):AND_MASK(0),OR_MASK(0),nmismatches(0)
        {
            for(int i=0;i<_seq.length();i++){
                char c=_seq[i];
                switch(c){
                    case 'A':
                        this->AND_MASK|=3;
                        //this->OR_MASK|=0;
                        this->nmismatches++;
                        positionMasks.push_back(getPositionMask(_seq.length(),i));
                        break;
                    case 'C':
                        this->AND_MASK|=3;
                        this->OR_MASK|=1;
                        this->nmismatches++;
                        positionMasks.push_back(getPositionMask(_seq.length(),i));
                        break;
                    case 'G':
                        this->AND_MASK|=3;
                        this->OR_MASK|=2;
                        this->nmismatches++;
                        positionMasks.push_back(getPositionMask(_seq.length(),i));
                        break;
                    case 'T':
                        this->AND_MASK|=3;
                        this->OR_MASK|=3;
                        this->nmismatches++;
                        positionMasks.push_back(getPositionMask(_seq.length(),i));
                        break;
                    default:
                        break;
                }
                if(i<_seq.length()-1){
                    this->OR_MASK<<=2;
                    this->AND_MASK<<=2;
                }

                

                
            }

            this->AND_MASK=~(this->AND_MASK);



        }

        string bitStringToString(uint64_t x){
            string ret;
            for(int i=0;i<64;i++){
                ret=((x&1)?"1":"0")+ret;
                x>>=1;
            }

            return ret;
        }

        string toString(){
            return bitStringToString(this->AND_MASK)+"\t"+bitStringToString(this->OR_MASK);
        }
        string ANDToString(){
            return bitStringToString(this->AND_MASK);
        }
        string ORToString(){
            return bitStringToString(this->OR_MASK);
        }
        string spacedString(string s){
            string ret=s.substr(0,1);
            for(int i=1;i<s.length();i++){
                ret+=" "+s.substr(i,1);
            }

            int l=s.length();

            for(int i=0;i<64-l*2;i++){
                ret=" "+ret;
            }
            return ret;
        }

};


class StepwiseSeqSpaceWithMotifs:public SeqSpace
{   
    private:
    int idx_i;

    uint64_t N_idx;

    uint64_t fwd_idx;
    uint64_t rev_idx;
    uint64_t mask;
    uint64_t rev_A_setter;
    uint64_t rev_C_setter;
    uint64_t rev_G_setter;
    uint64_t rev_T_setter;
    
    SeqBitMask fwd_mask;

    int motifLength;

    uint64_t nucI;

    public:

    void init_StepwiseSeqSpace(int _kmer)
    {
        this->idx_i=0;
        this->fwd_idx=0;
        this->rev_idx=0;

        this->N_idx=0;
  
        this->mask=3;

        
        this->rev_A_setter=0;
        this->rev_C_setter=1;
        this->rev_G_setter=2;
        this->rev_T_setter=3;

        for(int i=1;i<_kmer;i++)
        {
            this->mask<<=2;
            this->mask|=3;

            this->rev_A_setter<<=2;
            this->rev_C_setter<<=2;
            this->rev_G_setter<<=2;
            this->rev_T_setter<<=2;
        }
    }

    StepwiseSeqSpaceWithMotifs(int _kmer, string _fwd_motif):SeqSpace(_kmer),fwd_mask(_fwd_motif)
    {

        init_StepwiseSeqSpace(_fwd_motif.length());
        this->motifLength=_fwd_motif.length();
        //nucI=0;
    }

    //StepwiseSeqSpace(int _kmer,string inBitStringFileName, bool _useGzip=false):SeqSpace(_kmer,inBitStringFileName,_useGzip){
    //    init_StepwiseSeqSpace(_kmer);
    //}    

    

    string Uint64ToString(uint64_t idxer)
    {
        string s;
        for(int i=0;i<64;i++)
        {
            s=(((idxer&1)==1)?'1':'0')+s;
            idxer>>=1;
        }

        return s;
    }

    inline void see(char c)
    {
        nucI++;

        switch(c)
        {
            case 'A':

                this->idx_i++; 
                if(this->idx_i>1){
                    this->fwd_idx<<=2;
                    this->rev_idx>>=2;

                    this->N_idx<<=2;
                }

                this->rev_idx|=this->rev_T_setter;

                

                break;
            case 'C':

                this->idx_i++;
                if(this->idx_i>1){
                    this->fwd_idx<<=2;
                    this->rev_idx>>=2;

                    this->N_idx<<=2;
                }

                this->fwd_idx|=1;
                this->rev_idx|=this->rev_G_setter;

                break;
            case 'G':

                this->idx_i++;
                if(this->idx_i>1){
                    this->fwd_idx<<=2;
                    this->rev_idx>>=2;

                    this->N_idx<<=2;
                }

                this->fwd_idx|=2;
                this->rev_idx|=this->rev_C_setter;

                break;
            case 'T':

                this->idx_i++;
                if(this->idx_i>1){
                    this->fwd_idx<<=2;
                    this->rev_idx>>=2;

                    this->N_idx<<=2;
                }


                this->fwd_idx|=3;
                this->rev_idx|=this->rev_A_setter;

                break;
            case '>':
                //cerr<<"see ["<<c<<"]"<<endl;
                this->nucI=0;
                this->idx_i=0;
                this->fwd_idx=0;
                this->rev_idx=0;
                this->N_idx=0;

            case 'N':
            default:

                this->idx_i++;
                if(this->idx_i>1){
                    this->fwd_idx<<=2;
                    this->rev_idx>>=2;

                    this->N_idx<<=2;
                }

                this->N_idx|=3;

                break;
        }

        //string checkSeq("AAAATATGACTGGG");

        if(this->idx_i==this->motifLength)//this->kmer)
        {
            //now match fwd

            uint64_t masked_fwd_idx=this->fwd_idx & this->mask;
            if((masked_fwd_idx&(~this->fwd_mask.AND_MASK))==this->fwd_mask.OR_MASK){
                uint64_t actual_fwd_idx=masked_fwd_idx>>(2*(this->motifLength-this->kmer));

                //check N
                uint64_t masked_N_idx=this->N_idx & this->mask;
                uint64_t actual_N_idx=masked_N_idx>>(2*(this->motifLength-this->kmer));

                

                if(actual_N_idx==0)
                {
                    this->bits.setBit(actual_fwd_idx,1);
                
                    /*if(indexToSeq(actual_fwd_idx)==checkSeq){
                        cerr<<"SET FWD\t"<<StringUtil::str(int(nucI+1))<<"\t"<<indexToSeq64(actual_fwd_idx)<<endl;
                        //cerr<<"set Fwd bit"<<endl;
                        cerr<<"\tthis->fwd_idx\t"<<Uint64ToString(this->fwd_idx)<<"\t"<<indexToSeq64(this->fwd_idx)<<endl;
                        cerr<<"\tthis->mask\t"<<Uint64ToString(this->mask)<<"\t"<<indexToSeq64(this->mask)<<endl;
                        cerr<<"\tmasked_fwd_idx\t"<<Uint64ToString(masked_fwd_idx)<<"\t"<<indexToSeq64(masked_fwd_idx)<<endl;
                        cerr<<"\tactual_fwd_idx\t"<<Uint64ToString(actual_fwd_idx)<<"\t"<<indexToSeq64(actual_fwd_idx)<<endl;
                        cerr<<"\tget_bit at\t"<<Uint64ToString(actual_fwd_idx)<<"\t"<<indexToSeq(actual_fwd_idx)<<endl;
                    }*/

                }
                else
                {
                    /*cerr<<"has N FWD\t"<<StringUtil::str(int(nucI+1))<<"\t"<<indexToSeq64(actual_fwd_idx)<<endl;
                    cerr<<"\tactual_N_idx\t"<<Uint64ToString(actual_N_idx)<<"\t"<<indexToSeq64(actual_N_idx)<<endl;
                    cerr<<"\tthis->fwd_idx\t"<<Uint64ToString(this->fwd_idx)<<"\t"<<indexToSeq64(this->fwd_idx)<<endl;
                    cerr<<"\tthis->mask\t"<<Uint64ToString(this->mask)<<"\t"<<indexToSeq64(this->mask)<<endl;
                    cerr<<"\tmasked_fwd_idx\t"<<Uint64ToString(masked_fwd_idx)<<"\t"<<indexToSeq64(masked_fwd_idx)<<endl;
                    cerr<<"\tactual_fwd_idx\t"<<Uint64ToString(actual_fwd_idx)<<"\t"<<indexToSeq64(actual_fwd_idx)<<endl;*/                    
                }
                

            }

            uint64_t masked_rev_idx=this->rev_idx & this->mask;
            if((masked_rev_idx&(~this->fwd_mask.AND_MASK))==this->fwd_mask.OR_MASK){
                //uint64_t actual_rev_idx=masked_rev_idx&(this->mask>>(2*(this->motifLength-this->kmer)));
                uint64_t actual_rev_idx=masked_rev_idx>>(2*(this->motifLength-this->kmer));


                //check N
                uint64_t masked_N_idx=this->N_idx & this->mask;
                uint64_t actual_N_idx=masked_N_idx&(this->mask>>(2*(this->motifLength-this->kmer)));

                if(actual_N_idx==0)
                {
                    this->bits.setBit(actual_rev_idx,1);

                    /*if(indexToSeq(actual_rev_idx)==checkSeq){
                        
                        cerr<<"SET REV\t"<<StringUtil::str(int(nucI+1))<<"\t"<<indexToSeq64(actual_rev_idx)<<endl;
                        //cerr<<"set Rev bit"<<endl;
                        cerr<<"\tthis->rev_idx\t"<<Uint64ToString(this->rev_idx)<<"\t"<<indexToSeq64(this->rev_idx)<<endl;
                        cerr<<"\tthis->mask\t"<<Uint64ToString(this->mask)<<"\t"<<indexToSeq64(this->mask)<<endl;
                        cerr<<"\tmasked_rev_idx\t"<<Uint64ToString(masked_rev_idx)<<"\t"<<indexToSeq64(masked_rev_idx)<<endl;
                        cerr<<"\tactual_rev_idx\t"<<Uint64ToString(actual_rev_idx)<<"\t"<<indexToSeq64(actual_rev_idx)<<endl;
                        cerr<<"\tget_bit at\t"<<Uint64ToString(actual_rev_idx)<<"\t"<<indexToSeq(actual_rev_idx)<<endl;
                    }*/
                }
                else
                {
                    /*cerr<<"has N REV\t"<<StringUtil::str(int(nucI+1))<<"\t"<<indexToSeq64(actual_rev_idx)<<endl;
                    //cerr<<"set Rev bit"<<endl;
                    cerr<<"\tactual_N_idx\t"<<Uint64ToString(actual_N_idx)<<"\t"<<indexToSeq64(actual_N_idx)<<endl;
                    cerr<<"\tthis->rev_idx\t"<<Uint64ToString(this->rev_idx)<<"\t"<<indexToSeq64(this->rev_idx)<<endl;
                    cerr<<"\tthis->mask\t"<<Uint64ToString(this->mask)<<"\t"<<indexToSeq64(this->mask)<<endl;
                    cerr<<"\tmasked_rev_idx\t"<<Uint64ToString(masked_rev_idx)<<"\t"<<indexToSeq64(masked_rev_idx)<<endl;
                    cerr<<"\tactual_rev_idx\t"<<Uint64ToString(actual_rev_idx)<<"\t"<<indexToSeq64(actual_rev_idx)<<endl;*/                    
                }
                

            }

            //cerr<<"\tFwd "<<Uint64ToString(this->fwd_idx & this->mask)<<" ~ "<<this->indexToSeq(this->fwd_idx & this->mask)<<endl;
            //cerr<<"\tRev "<<Uint64ToString(this->rev_idx & this->mask)<<" ~ "<<this->indexToSeq(this->rev_idx & this->mask)<<endl;
            this->idx_i--;
        }

    }

    static string getFwdMotif(int kmer,string PAM){
        string motif;
        for(int i=0;i<kmer;i++){
            motif+="N";
        }

        return motif+PAM;
    }

    static string getRevMotif(string rcPAM,int kmer){
        string motif=rcPAM;
        for(int i=0;i<kmer;i++){
            motif+="N";
        }

        return motif;
    }

    

};


class OffTargetEnumerator:public StepwiseSeqSpace
{
    class TaskMem{
        public:
            string misString;
            int level;
            int pos;
            TaskMem(string _misString,int _level, int _pos):misString(_misString),level(_level),pos(_pos){

            }
        
    };

    public:
        vector<SeqBitMask> seqBitMasks;
        int mthres;
        int *results;

        ~OffTargetEnumerator(){
            delete[] results;
        }

        void init_OffTargetEnumerator(int _kmer,int _mthres){
           results=new int[_mthres];
            for(int i=0;i<=_mthres;i++){
                results[i]=0;
            }

            queue<TaskMem> Q;
            string misString;
            string alphabet("ACGT");            

            for(int i=0;i<_kmer;i++){
                misString+="N";
            }

            for(int i=0;i<_kmer;i++){
                Q.push(TaskMem(misString,1,i));
            }

            int counter=0;
            
            while(!Q.empty()){
                TaskMem thisTask=Q.front();
                Q.pop();
                for(int k=0;k<4;k++){
                    
                    string newMisString(thisTask.misString);
                    newMisString[thisTask.pos]=alphabet[k];



                    counter++;
                    
                    SeqBitMask seqBitMask(newMisString);
                    this->seqBitMasks.push_back(seqBitMask);
                    //cerr<<counter<<"\t\t"<<newMisString<<"\t"<<seqBitMask.toString()<<endl; //do sth else;
                    //cerr<<counter<<"\t\t"<<seqBitMask.spacedString(newMisString)<<endl;
                    //cerr<<counter<<"\t\t"<<seqBitMask.ANDToString()<<endl;
                    //cerr<<counter<<"\t\t"<<seqBitMask.ORToString()<<endl;
                    
                    

                    if(thisTask.level<mthres){
                        for(int j=thisTask.pos+1;j<_kmer;j++){
                            TaskMem v(newMisString,thisTask.level+1,j);
                            Q.push(v);
                        }
                        
                    }
                }
            }
        }


        OffTargetEnumerator(int _kmer, string _inBitStringFileName, int _mthres, bool useGzip=false):StepwiseSeqSpace(_kmer,_inBitStringFileName,useGzip),mthres(_mthres)
        {
            init_OffTargetEnumerator(_kmer,_mthres);
        }

        OffTargetEnumerator(int _kmer,int _mthres):StepwiseSeqSpace(_kmer),mthres(_mthres)
        {
            init_OffTargetEnumerator(_kmer,_mthres);
        }




        inline void findOffTargetHits(const string &seq)
        {
            uint64_t seqIdx=this->seqToIndex(seq);
            findOffTargetHits(seqIdx);
        }

        inline void findOffTargetHits(uint64_t seqIdx, vector<uint64_t> *__mthres=NULL)
        {
            for(int i=0;i<=this->mthres;i++){
                results[i]=0;
            }

            

            if(this->bits.getBit(seqIdx)){
                results[0]=1;
            }

            for(vector<SeqBitMask>::iterator i=this->seqBitMasks.begin();i!=this->seqBitMasks.end();i++){

                bool doubleCounting=false;
                for(vector<uint64_t>::iterator it=i->positionMasks.begin();it!=i->positionMasks.end();it++){
                    uint64_t mask=(*it);
                    if((i->OR_MASK&mask) == (seqIdx&mask)){
                        doubleCounting=true;
                        break;
                    }
                    
                }

                if(doubleCounting){
                    continue;
                }


                uint64_t searchIdx=seqIdx;
                searchIdx&=i->AND_MASK;
                searchIdx|=i->OR_MASK;


                if(searchIdx!=seqIdx && this->bits.getBit(searchIdx)){
                    results[i->nmismatches]++;
                    if(__mthres){
                        if(results[i->nmismatches]>(*__mthres)[i->nmismatches]){
                            break;
                        }
                    }
                }
            }
        }

        string getResultString(){
            string ret;
            ret=StringUtil::str(results[0]);
            for(int i=1;i<=this->mthres;i++)
            {
                ret+="/"+StringUtil::str(results[i]);
            }
            return ret;
        }

        /*void enumerateAbsentSitesWithOffTargetProfiles(ostream &os)
        {
            uint64_t output_count=0;

            for(uint64_t i=0;i<bits.numBits;i++){

                if(!bits.getBit(i)){
                    os<<this->indexToSeq(i)<<"\t"<<this->findOffTargetHits(i)<<endl;
                    output_count++;
                }
            }

        }*/

        void enumerateAbsentSitesWithOffTargetProfiles(ostream &os)
        {
            uint64_t output_count=0;

            for(uint64_t i=0;i<bits.numBits;i++){

                if(!bits.getBit(i)){
                    this->findOffTargetHits(i);
                    os<<this->indexToSeq(i)<<"\t"<<this->getResultString()<<endl;
                    output_count++;
                }
            }

            cerr<<output_count<<" VaKmers outputted"<<endl;

        }

        bool enumerateAbsentSitesWithOffTargetProfiles(ostream &os,vector<uint64_t> __mthres,string seqPattern,uint64_t maxNumToOutput)
        {

            if(__mthres.size()!=this->mthres+1){
                return false;
            }

            if(seqPattern.length()!=this->kmer){
                return false;
            }

            SeqBitMask seqPatternMask(seqPattern);

            uint64_t output_count=0;

            //cerr<<"seqPattern"<<"\n"<<seqPatternMask.ANDToString()<<"\n"<<seqPatternMask.ORToString()<<endl;

            //return false;

            //cerr<<"noofbits "<<bits.numBits<<endl;

            for(uint64_t i=0;i<bits.numBits;i++){
                
                if((i&(~seqPatternMask.AND_MASK))!=seqPatternMask.OR_MASK){
                    //cerr<<"not match:"<<this->indexToSeq(i)<<" "<<seqPatternMask.bitStringToString(i&seqPatternMask.AND_MASK)<<endl;
                    //cerr<<"not match:"<<this->indexToSeq(i)<<" "<<seqPatternMask.bitStringToString(i)<<endl;
                    continue;        //does not match seqPattern
                    
                }

                if(!bits.getBit(i)){
                    


                    this->findOffTargetHits(i,&__mthres);
                    
                    //os<<"test "<<this->indexToSeq(i)<<"\t"<<this->getResultString()<<endl;

                    bool outputThis=true;
                    for(int j=0;j<__mthres.size();j++){
                        if(this->results[j]>__mthres[j]){
                            outputThis=false;
                            break;
                        }
                    }

                    if(outputThis){
                        os<<this->indexToSeq(i)<<"\t"<<this->getResultString()<<endl;
                        output_count++;
                    }

                    if(maxNumToOutput>0 && output_count>=maxNumToOutput){
                        break; //enough
                    }
                    
                }
            }

            cerr<<output_count<<" VaKmers outputted"<<endl;
            return true;
        }

        bool enumerateAbsentSitesWithOffTargetProfilesForPrefixedSeq(ostream &os,vector<uint64_t> __mthres,string seqPrefix,uint64_t maxNumToOutput)
        {

            if(__mthres.size()!=this->mthres+1){
                return false;
            }

            string lowBound=seqPrefix;
            string highBound=seqPrefix;

            for(int i=seqPrefix.length();i<this->kmer;i++){
                lowBound+="A";
                highBound+="T";
            }

            uint64_t lowBoundIdx=this->seqToIndex(lowBound);
            uint64_t highBoundIdx=this->seqToIndex(highBound);

            cerr<<"search ["<<lowBound<<"/"<<lowBoundIdx<<","<<highBound<<"/"<<highBoundIdx<<"]"<<endl;

            uint64_t output_count=0;

            //cerr<<"seqPattern"<<"\n"<<seqPatternMask.ANDToString()<<"\n"<<seqPatternMask.ORToString()<<endl;

            //return false;

            //cerr<<"noofbits "<<bits.numBits<<endl;

            for(uint64_t i=lowBoundIdx;i<=highBoundIdx;i++){
                
            

                if(!bits.getBit(i)){
                    


                    this->findOffTargetHits(i,&__mthres);
                    
                    //os<<"test "<<this->indexToSeq(i)<<"\t"<<this->getResultString()<<endl;

                    bool outputThis=true;
                    for(int j=0;j<__mthres.size();j++){
                        if(this->results[j]>__mthres[j]){
                            outputThis=false;
                            break;
                        }
                    }

                    if(outputThis){
                        os<<this->indexToSeq(i)<<"\t"<<this->getResultString()<<endl;
                        output_count++;
                    }

                    if(maxNumToOutput>0 && output_count>=maxNumToOutput){
                        break; //enough
                    }
                    
                }
            }

            cerr<<output_count<<" VaKmers outputted"<<endl;
            return true;
        }

        
        /*uint64_t simulate31220Ops(){
            uint64_t counter=0;
            for(uint64_t i=0;i<31220;i++){
                if(this->bits.getBit(i)){
                    counter++;
                }
            }

            return counter;
        }*/


};

#ifdef __SEQSPACE_COMPARE
int main(int argc,char **argv)
{
    cerr<<"SEQSPACE COMPARE"<<endl;
    
    if (argc<4)
    {
        cerr<<"Usage: "<<argv[0]<<" seqBaseFile1 seqBaseFile2 kmer"<<endl;
        exit(1);
    }
    
    string filename1=argv[1];
    string filename2=argv[2];
    int kmer=StringUtil::atoi(argv[3]);



    time_t prev_time, now_time;
    prev_time=time(NULL);
    SeqSpace seqspace1(kmer, filename1,StringUtil::toUpper(filename1.substr(filename1.length()-2,2))=="GZ" );
    now_time=time(NULL);
    cerr<<"time1 spent:"<<(now_time-prev_time)<<endl;
    prev_time=time(NULL);
    SeqSpace seqspace2(kmer, filename2, StringUtil::toUpper(filename2.substr(filename2.length()-2,2))=="GZ");
    now_time=time(NULL);
    cerr<<"time2 spent:"<<(now_time-prev_time)<<endl;
    //cerr<<"bistrings are "<<((e1.bits==e2.bits)?"identical":"different")<<endl;
    bool same=seqspace1.diff(seqspace2);
    cerr<<"bistrings are "<<(same?"identical":"different")<<endl;
    cerr<<"Done"<<endl;
    return 0;
}

#endif //__SEQSPACE_COMPARE


#ifdef __MODE1__

int main(int argc,char **argv)
{
    string extra_settings="";

    #ifdef __FAST_IO__
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    extra_settings=" FastIO option=True";
    #endif //__FAST_IO__

    cerr<<"VaKation (\e[4mVa\e[0mcant \e[4mK\e[0m-mers Identific\e[4mation\e[0m) vers 0.1 mode 1 "<<extra_settings<<endl;

    if(argc<3){
        
        return printUsageAndExit(argv[0]);
    }
    
    int kmer=StringUtil::atoi(argv[1]);
    
    time_t start_time=time(NULL);
   
    SeqSpace seqspace(kmer);
    
    
    
    for(int i=2;i<argc;i++){
        

       time_t prev_time=time(NULL);
        
        
        cerr<<"processing fasta file "<<argv[i]<<endl;
        FastaFile fastaFile(argv[i]);
        while (fastaFile.readEntry()){
            cerr<<"processing sequence "<<fastaFile.seqName<<endl;
            seqspace.enumRegisterSeq(fastaFile.seq);
        }

        time_t now_time=time(NULL);
        cerr<<"Finish "<<argv[i]<<" in "<<(now_time-prev_time)<<" second(s)"<<endl;

        
    }

    time_t now_time=time(NULL);
    cerr<<"Finish processing all fasta files in "<<(now_time-start_time)<<" second(s)"<<endl;   
    
    cerr<<"Start printing VaKmers"<<endl;
    time_t prev_time=time(NULL);
    seqspace.printAbsentKmers(cout);
    now_time=time(NULL);
    cerr<<"Finish printing VaKmers in "<<(now_time-prev_time)<<" second(s)"<<endl;
    cerr<<"Done in "<<(now_time-start_time)<<" second(s)"<<endl;
    
    
    return 0;
    

}

#endif //__MODE1__



#ifdef __MODE1_WRITERNGG

int main(int argc,char **argv)
{
    string extra_settings="";

    #ifdef __FAST_IO__
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    extra_settings+=" FastIO option=True";
    #endif //__FAST_IO__

    extra_settings+=" GZIP=True";

    cerr<<"VaKation (\e[4mVa\e[0mcant \e[4mK\e[0m-mers Identific\e[4mation\e[0m) vers 0.1 mode 1 Writer "<<extra_settings<<endl;

    if(argc<3){
        
        return printUsageAndExit_writer(argv[0]);
    }
    
    string outBitStringFileName=argv[1];

    int kmer=StringUtil::atoi(argv[2]);
    
    time_t start_time=time(NULL);
   
    SeqSpace seqspace(kmer);
    
    
    
    for(int i=3;i<argc;i++){
        

       time_t prev_time=time(NULL);
        
        
        cerr<<"processing fasta file "<<argv[i]<<endl;
        FastaFile fastaFile(argv[i]);
        while (fastaFile.readEntry()){
            cerr<<"processing sequence "<<fastaFile.seqName<<endl;
            seqspace.enumRegisterSeqNGG(fastaFile.seq);
        }

        time_t now_time=time(NULL);
        cerr<<"Finish "<<argv[i]<<" in "<<(now_time-prev_time)<<" second(s)"<<endl;

        
    }

    time_t now_time=time(NULL);
    cerr<<"Finish processing all fasta files in "<<(now_time-start_time)<<" second(s)"<<endl;   
    
    cerr<<"Start writing Bitstring"<<endl;
    time_t prev_time=time(NULL);

    if(StringUtil::toUpper(outBitStringFileName.substr(outBitStringFileName.length()-2,2))=="GZ")
        seqspace.bits.writeToFileGZ(outBitStringFileName);
    else
        seqspace.bits.writeToFile(outBitStringFileName);
    
    now_time=time(NULL);
    cerr<<"Finish writing Bitstring in "<<(now_time-prev_time)<<" second(s)"<<endl;
    cerr<<"Done in "<<(now_time-start_time)<<" second(s)"<<endl;
    
    
    return 0;
    

}

#endif //__MODE1_WRITERNGG




#ifdef __MODE2__

int main(int argc,char **argv)
{
    string extra_settings="";

    #ifdef __FAST_IO__
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    extra_settings=" FastIO option=True";
    #endif //__FAST_IO__

    cerr<<"VaKation (\e[4mVa\e[0mcant \e[4mK\e[0m-mers Identific\e[4mation\e[0m) vers 0.1 mode 2 "<<extra_settings<<endl;

    if(argc<3){
        
        return printUsageAndExit(argv[0]);
    }
    
    int kmer=StringUtil::atoi(argv[1]);
    
    time_t start_time=time(NULL);

    StepwiseSeqSpace seqspace(kmer);
    
    
    
    for(int i=2;i<argc;i++){

        time_t prev_time=time(NULL);

        cerr<<"processing fasta file "<<argv[i]<<endl;
        FastaFile fastaFile(argv[i]);

        while (fastaFile.readEntry()){
            cerr<<"processing sequence "<<fastaFile.seqName<<endl;
            seqspace.see('N'); //reinit 

            while(true)
            {
                char c=toupper(fastaFile.readNextNtFromLoadedSeq());
                
                
                if(c=='\0'){
                    break;
                }else{
                    //cerr<<c<<endl;
                    seqspace.see(c);
                    //cerr<<c<<endl;
                }
            }
        }

        time_t now_time=time(NULL);
        cerr<<"Finish "<<argv[i]<<" in "<<(now_time-prev_time)<<" second(s)"<<endl;

        
    }

    time_t now_time=time(NULL);
    cerr<<"Finish processing all fasta files in "<<(now_time-start_time)<<" second(s)"<<endl;

    cerr<<"Start printing VaKmers"<<endl;
    time_t prev_time=time(NULL);
    seqspace.printAbsentKmers(cout);
    now_time=time(NULL);
    cerr<<"Finish printing VaKmers in "<<(now_time-prev_time)<<" second(s)"<<endl;
    cerr<<"Done in "<<(now_time-start_time)<<" second(s)"<<endl;

    
    
    return 0;
    

}

#endif //__MODE2__

#ifdef __MODE3_OFFAWARE2

int main(int argc,char **argv)
{

    string extra_settings="";

    #ifdef __FAST_IO__
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    extra_settings=" FastIO option=True";
    #endif //__FAST_IO__

    cerr<<"VaKation (\e[4mVa\e[0mcant \e[4mK\e[0m-mers Identific\e[4mation\e[0m) vers 0.1 mode 3 offAware "<<extra_settings<<endl;
    
    if(argc<3){
        
        return printUsageAndExit_Mode3Offaware2(argv[0]);
    }
    
    int kmer=StringUtil::atoi(argv[1]);
    int mthres=StringUtil::atoi(argv[2]);


    vector<string> splits;
    vector<uint64_t> mthress; 
    StringUtil::split(argv[3],"/",splits);
    for(vector<string>::iterator i=splits.begin();i!=splits.end();i++){
        mthress.push_back(StringUtil::atoi(*i));
    }

    if(mthress.size()!=mthres+1){
        cerr<<"incorrect number of mthress. Abort"<<endl;
        return printUsageAndExit_Mode3Offaware2(argv[0]);
    }
    //enumerateAbsentSitesWithOffTargetProfiles(ostream &os,vector<uint64_t> __mthres,string seqPattern,uint64_t maxNumToOutput)
    string seqPrefix=argv[4];

    if(seqPrefix.length()>kmer){
        cerr<<"incorrect seqPrefix length>kmer. Abort"<<endl;
        return printUsageAndExit_Mode3Offaware2(argv[0]);
    }

    uint64_t maxNumToOutput=StringUtil::atoi(argv[5]);

    time_t start_time=time(NULL);

    OffTargetEnumerator seqspace(kmer,mthres);
    
    
    
    for(int i=6;i<argc;i++){
        

        time_t prev_time=time(NULL);
        
        cerr<<"Processing fasta file "<<argv[i]<<endl;
        FastaFile fastaFile(argv[i]);
        while (true){

            char c=toupper(fastaFile.readNextNt());
            if(c=='\0'){
                break;
            }else if(c=='>'){
                cerr<<"Processing sequence "<<fastaFile.seqName<<endl;
                seqspace.see('N');
            }else{
                //cerr<<c<<endl;
                seqspace.see(c);
                //cerr<<c<<endl;
            }
        }

        time_t now_time=time(NULL);
        cerr<<"Finish "<<argv[i]<<" in "<<(now_time-prev_time)<<" second(s)"<<endl;
        
    }

    time_t now_time=time(NULL);
    cerr<<"Finish processing all fasta files in "<<(now_time-start_time)<<" second(s)"<<endl;

    cerr<<"Start printing VaKmers"<<endl;
    time_t prev_time=time(NULL);
    seqspace.enumerateAbsentSitesWithOffTargetProfilesForPrefixedSeq(cout,mthress,seqPrefix,maxNumToOutput);
    now_time=time(NULL);
    cerr<<"Finish printing VaKmers in "<<(now_time-prev_time)<<" second(s)"<<endl;
    cerr<<"Done in "<<(now_time-start_time)<<" second(s)"<<endl;
    return 0;
    

}

#endif //__MODE3_OFFAWARE


#ifdef __MODE3_OFFAWARE

int main(int argc,char **argv)
{

    string extra_settings="";

    #ifdef __FAST_IO__
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    extra_settings=" FastIO option=True";
    #endif //__FAST_IO__

    //cerr<<"VaKation (\e[4mVa\e[0mcant \e[4mK\e[0m-mers Identific\e[4mation\e[0m) vers 0.1 mode 3 offAware "<<extra_settings<<endl;
    cerr<<"JACKIE.enumAbsentmers"<<extra_settings<<endl;

    if(argc<3){
        
        return printUsageAndExit_Mode3Offaware(argv[0]);
    }
    
    int kmer=StringUtil::atoi(argv[1]);
    int mthres=StringUtil::atoi(argv[2]);


    vector<string> splits;
    vector<uint64_t> mthress; 
    StringUtil::split(argv[3],"/",splits);
    for(vector<string>::iterator i=splits.begin();i!=splits.end();i++){
        mthress.push_back(StringUtil::atoi(*i));
    }

    if(mthress.size()!=mthres+1){
        cerr<<"incorrect number of mthress. Abort"<<endl;
        return printUsageAndExit_Mode3Offaware(argv[0]);
    }
    //enumerateAbsentSitesWithOffTargetProfiles(ostream &os,vector<uint64_t> __mthres,string seqPattern,uint64_t maxNumToOutput)
    string seqPattern=argv[4];

    if(seqPattern.length()!=kmer){
        cerr<<"incorrect seqPattern length!=kmer. Abort"<<endl;
        return printUsageAndExit_Mode3Offaware(argv[0]);
    }

    uint64_t maxNumToOutput=StringUtil::atoi(argv[5]);

    time_t start_time=time(NULL);

    OffTargetEnumerator seqspace(kmer,mthres);
    
    
    
    for(int i=6;i<argc;i++){
        

        time_t prev_time=time(NULL);
        
        cerr<<"Processing fasta file "<<argv[i]<<endl;
        FastaFile fastaFile(argv[i]);
        while (true){

            char c=toupper(fastaFile.readNextNt());
            if(c=='\0'){
                break;
            }else if(c=='>'){
                cerr<<"Processing sequence "<<fastaFile.seqName<<endl;
                seqspace.see('N');
            }else{
                //cerr<<c<<endl;
                seqspace.see(c);
                //cerr<<c<<endl;
            }
        }

        time_t now_time=time(NULL);
        cerr<<"Finish "<<argv[i]<<" in "<<(now_time-prev_time)<<" second(s)"<<endl;
        
    }

    time_t now_time=time(NULL);
    cerr<<"Finish processing all fasta files in "<<(now_time-start_time)<<" second(s)"<<endl;

    cerr<<"Start printing VaKmers"<<endl;
    time_t prev_time=time(NULL);
    seqspace.enumerateAbsentSitesWithOffTargetProfiles(cout,mthress,seqPattern,maxNumToOutput);
    now_time=time(NULL);
    cerr<<"Finish printing VaKmers in "<<(now_time-prev_time)<<" second(s)"<<endl;
    cerr<<"Done in "<<(now_time-start_time)<<" second(s)"<<endl;
    return 0;
    

}

#endif //__MODE3_OFFAWARE




#ifdef __MODE3__

int main(int argc,char **argv)
{

    string extra_settings="";

    #ifdef __FAST_IO__
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    extra_settings=" FastIO option=True";
    #endif //__FAST_IO__

    cerr<<"VaKation (\e[4mVa\e[0mcant \e[4mK\e[0m-mers Identific\e[4mation\e[0m) vers 0.1 mode 3 "<<extra_settings<<endl;
    
    if(argc<3){
        
        return printUsageAndExit(argv[0]);
    }
    
    int kmer=StringUtil::atoi(argv[1]);
    

    time_t start_time=time(NULL);

    StepwiseSeqSpace seqspace(kmer);
    
    
    
    for(int i=2;i<argc;i++){
        

        time_t prev_time=time(NULL);
        
        cerr<<"Processing fasta file "<<argv[i]<<endl;
        FastaFile fastaFile(argv[i]);
        while (true){

            char c=toupper(fastaFile.readNextNt());
            if(c=='\0'){
                break;
            }else if(c=='>'){
                cerr<<"Processing sequence "<<fastaFile.seqName<<endl;
                seqspace.see('N');
            }else{
                //cerr<<c<<endl;
                seqspace.see(c);
                //cerr<<c<<endl;
            }
        }

        time_t now_time=time(NULL);
        cerr<<"Finish "<<argv[i]<<" in "<<(now_time-prev_time)<<" second(s)"<<endl;
        
    }

    time_t now_time=time(NULL);
    cerr<<"Finish processing all fasta files in "<<(now_time-start_time)<<" second(s)"<<endl;

    cerr<<"Start printing VaKmers"<<endl;
    time_t prev_time=time(NULL);
    seqspace.printAbsentKmers(cout);
    now_time=time(NULL);
    cerr<<"Finish printing VaKmers in "<<(now_time-prev_time)<<" second(s)"<<endl;
    cerr<<"Done in "<<(now_time-start_time)<<" second(s)"<<endl;
    return 0;
    

}

#endif //__MODE3__

#ifdef __MODE3_WRITER


int main(int argc,char **argv)
{

    string extra_settings="";

    #ifdef __FAST_IO__
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    extra_settings=" FastIO option=True";
    #endif //__FAST_IO__

    extra_settings+=" GZIP=True";

    //cerr<<"VaKation (\e[4mVa\e[0mcant \e[4mK\e[0m-mers Identific\e[4mation\e[0m) vers 0.1 mode 3Write "<<extra_settings<<endl;
    cerr<<"JACKIE.encodeSeqSpace"<<extra_settings<<endl;

    if(argc<3){
        
        return printUsageAndExit_writer(argv[0]);
    }
    
    string outBitStringFileName=argv[1];

    int kmer=StringUtil::atoi(argv[2]);
    

    time_t start_time=time(NULL);

    StepwiseSeqSpace seqspace(kmer);
    
    
    
    for(int i=3;i<argc;i++){
        

        time_t prev_time=time(NULL);
        
        cerr<<"Processing fasta file "<<argv[i]<<endl;
        FastaFile fastaFile(argv[i]);
        while (true){

            char c=toupper(fastaFile.readNextNt());
            if(c=='\0'){
                break;
            }else if(c=='>'){
                cerr<<"Processing sequence "<<fastaFile.seqName<<endl;
                seqspace.see('N');
            }else{
                //cerr<<c<<endl;
                seqspace.see(c);
                //cerr<<c<<endl;
            }
        }

        time_t now_time=time(NULL);
        cerr<<"Finish "<<argv[i]<<" in "<<(now_time-prev_time)<<" second(s)"<<endl;
        
    }

    time_t now_time=time(NULL);
    cerr<<"Finish processing all fasta files in "<<(now_time-start_time)<<" second(s)"<<endl;

    cerr<<"Start printing VaKmers"<<endl;
    time_t prev_time=time(NULL);
    if(StringUtil::toUpper(outBitStringFileName.substr(outBitStringFileName.length()-2,2))=="GZ")
        seqspace.bits.writeToFileGZ(outBitStringFileName);
    else
        seqspace.bits.writeToFile(outBitStringFileName);
    now_time=time(NULL);
    cerr<<"Finish printing VaKmers in "<<(now_time-prev_time)<<" second(s)"<<endl;
    cerr<<"Done in "<<(now_time-start_time)<<" second(s)"<<endl;
    return 0;
    

}

#endif //__MODE3_WRITER


#ifdef __MODE3_WRITERNGG

int main(int argc,char **argv)
{

    string extra_settings="";

    #ifdef __FAST_IO__
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    extra_settings=" FastIO option=True";
    #endif //__FAST_IO__

    extra_settings+=" GZIP=True";

    //cerr<<"VaKation (\e[4mVa\e[0mcant \e[4mK\e[0m-mers Identific\e[4mation\e[0m) vers 0.1 mode 3Write NGG "<<extra_settings<<endl;
    cerr<<"JACKIE.encodeSeqSpaceNGG"<<extra_settings<<endl;


    if(argc<3){
        
        return printUsageAndExit_writer(argv[0]);
    }
    
    string outBitStringFileName=argv[1];

    int kmer=StringUtil::atoi(argv[2]);
    

    time_t start_time=time(NULL);

    StepwiseSeqSpaceWithMotifs seqspace(kmer,StepwiseSeqSpaceWithMotifs::getFwdMotif(kmer,"NGG"));
    
    
    
    for(int i=3;i<argc;i++){
        

        time_t prev_time=time(NULL);
        
        cerr<<"Processing fasta file "<<argv[i]<<endl;
        FastaFile fastaFile(argv[i]);
        while (true){

            char c=toupper(fastaFile.readNextNt());
            if(c=='\0'){
                break;
            }else if(c=='>'){
                cerr<<"Processing sequence "<<fastaFile.seqName<<endl;
                seqspace.see('>');
            }else{
                //cerr<<c<<endl;
                seqspace.see(c);
                //cerr<<c<<endl;
            }
        }

        time_t now_time=time(NULL);
        cerr<<"Finish "<<argv[i]<<" in "<<(now_time-prev_time)<<" second(s)"<<endl;
        
    }

    time_t now_time=time(NULL);
    cerr<<"Finish processing all fasta files in "<<(now_time-start_time)<<" second(s)"<<endl;

    cerr<<"Start printing VaKmers"<<endl;
    time_t prev_time=time(NULL);
    //seqspace.bits.writeToFile(outBitStringFileName);
    
    if(StringUtil::toUpper(outBitStringFileName.substr(outBitStringFileName.length()-2,2))=="GZ")
        seqspace.bits.writeToFileGZ(outBitStringFileName);
    else
        seqspace.bits.writeToFile(outBitStringFileName);

    now_time=time(NULL);
    cerr<<"Finish printing VaKmers in "<<(now_time-prev_time)<<" second(s)"<<endl;
    cerr<<"Done in "<<(now_time-start_time)<<" second(s)"<<endl;
    return 0;
    

}

#endif //__MODE3_WRITERNGG


#ifdef __MODE3_READER


int main(int argc,char **argv)
{

    string extra_settings="";

    #ifdef __FAST_IO__
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    extra_settings+=" FastIO option=True";
    #endif //__FAST_IO__

    #ifdef __USE_GZIP
    extra_settings+=" GZIP=True";
    #endif

    //cerr<<"VaKation (\e[4mVa\e[0mcant \e[4mK\e[0m-mers Identific\e[4mation\e[0m) vers 0.1 mode 3Read "<<extra_settings<<endl;
    cerr<<"JACKIE.countSeqNeighbors "<<extra_settings<<endl;


    if(argc<3){
        
        return printUsageAndExit_reader(argv[0]);
    }
    
    string inBitStringFileName=argv[1];

    int kmer=StringUtil::atoi(argv[2]);
    int mthres=StringUtil::atoi(argv[3]);

    string searchListFileName=argv[4];
  
    int col0=-1;
    string sep="";
    int splitCom0=-1;

    if(argc>=6){
        string argv_string(argv[5]);
        vector<string> v;
        StringUtil::split(argv_string,",",v);
        col0=StringUtil::atoi(v[0])-1;
        if(v.size()>1)
        {
            sep=v[1];
            splitCom0=StringUtil::atoi(v[2])-1;
        }

    }



    time_t start_time=time(NULL);
    time_t now_time=time(NULL);
    
    /*#ifdef __USE_GZIP
    bool useGzip=true;
    #else
    bool useGzip=false;
    #endif*/

    cerr<<"Start reading bitsring file"<<endl;
    time_t prev_time=time(NULL);

    bool useGzip=(StringUtil::toUpper(inBitStringFileName.substr(inBitStringFileName.length()-2,2))=="GZ");
   
    OffTargetEnumerator enumerator(kmer,inBitStringFileName,mthres,useGzip);

 

    now_time=time(NULL);
    cerr<<"Finish reading in "<<(now_time-prev_time)<<" second(s)"<<endl;
    prev_time=time(NULL);
    uint64_t counter=0;

    ifstream fin;
    fin.open(searchListFileName.c_str());

    uint64_t lino=0;
    while(!fin.eof()){
        lino++;
        if(lino%1000000==1){
            now_time=time(NULL);
            cerr<<"processed "<<(lino-1)<<" lines, time elapsed:"<<(now_time-prev_time)<<endl;
        }
        string lin;
        getline(fin,lin);
        if(lin!=""){

            if(col0>=0){
                vector<string> fields;
                StringUtil::split(lin,"\t",fields);
                string col_value=fields[col0];
                if(splitCom0>=0){
                    vector<string> splits;
                    StringUtil::split(col_value,sep,splits);
                    enumerator.findOffTargetHits(splits[splitCom0]);
                }else{
                    enumerator.findOffTargetHits(col_value);
                }
            }else{
                enumerator.findOffTargetHits(lin);
            }
            

            
            cout<<lin<<"\t"<<enumerator.getResultString()<<endl;
        }
    }

    now_time=time(NULL);
    cerr<<"Finish finding off-targets for "<<lino<< "lines in "<<(now_time-prev_time)<<" second(s)"<<endl;

    cerr<<"Done in "<<(now_time-start_time)<<" second(s)"<<endl;
    return 0;
    

}

#endif //__MODE3_READER

#ifdef __MODE3_SAMPLER

int main(int argc,char **argv)
{

    string extra_settings="";

    #ifdef __FAST_IO__
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    extra_settings=" FastIO option=True";
    #endif //__FAST_IO__

    cerr<<"VaKation (\e[4mVa\e[0mcant \e[4mK\e[0m-mers Identific\e[4mation\e[0m) vers 0.1 mode 3 Sampler "<<extra_settings<<endl;
    
    if(argc<3){
        
        return printUsageAndExit(string(argv[0])+" samplingCount");
    }
    
    int samplingCount=StringUtil::atoi(argv[1]);
    int kmer=StringUtil::atoi(argv[2]);
    

    time_t start_time=time(NULL);

    StepwiseSeqSpace seqspace(kmer);
    
    
    
    for(int i=3;i<argc;i++){
        

        time_t prev_time=time(NULL);
        
        cerr<<"Processing fasta file "<<argv[i]<<endl;
        FastaFile fastaFile(argv[i]);
        while (true){

            char c=toupper(fastaFile.readNextNt());
            if(c=='\0'){
                break;
            }else if(c=='>'){
                cerr<<"Processing sequence "<<fastaFile.seqName<<endl;
                seqspace.see('N');
            }else{
                //cerr<<c<<endl;
                seqspace.see(c);
                //cerr<<c<<endl;
            }
        }

        time_t now_time=time(NULL);
        cerr<<"Finish "<<argv[i]<<" in "<<(now_time-prev_time)<<" second(s)"<<endl;
        
    }

    time_t now_time=time(NULL);
    cerr<<"Finish processing all fasta files in "<<(now_time-start_time)<<" second(s)"<<endl;
    
    srand(time(NULL)); //seed random number generator

    cerr<<"Start printing VaKmers"<<endl;
    time_t prev_time=time(NULL);
    seqspace.sampleAbsentKmers(cout,samplingCount);
    now_time=time(NULL);
    cerr<<"Finish printing VaKmers in "<<(now_time-prev_time)<<" second(s)"<<endl;
    cerr<<"Done in "<<(now_time-start_time)<<" second(s)"<<endl;
    return 0;
    

}

#endif //__MODE3_SAMPLER


#ifdef __MODE4__

int main(int argc,char **argv)
{

    string extra_settings="";

    #ifdef __FAST_IO__
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    extra_settings=" FastIO option=True";
    #endif //__FAST_IO__

    cerr<<"VaKation (\e[4mVa\e[0mcant \e[4mK\e[0m-mers Identific\e[4mation\e[0m) vers 0.1 mode 4 "<<extra_settings<<endl;
    


    if(argc<3){
        
        return printUsageAndExit(argv[0]);
    }
    
    int kmer=StringUtil::atoi(argv[1]);
    

    time_t start_time=time(NULL);

    StepwiseSeqSpace seqspace(kmer);
    
    
    
    for(int i=2;i<argc;i++){
        

        time_t prev_time=time(NULL);
        
        cerr<<"Processing fasta file "<<argv[i]<<endl;
        FastaFile fastaFile(argv[i]);
        while (true){

            char c=toupper(fastaFile.readNextNt());
            if(c=='\0'){
                break;
            }else if(c=='>'){
                cerr<<"Processing sequence "<<fastaFile.seqName<<endl;
                seqspace.see('N');
            }else{
                //cerr<<c<<endl;
                seqspace.see(c);
                //cerr<<c<<endl;
            }
        }

        time_t now_time=time(NULL);
        cerr<<"Finish "<<argv[i]<<" in "<<(now_time-prev_time)<<" second(s)"<<endl;
        
    }

    time_t now_time=time(NULL);
    cerr<<"Finish processing all fasta files in "<<(now_time-start_time)<<" second(s)"<<endl;

    
    cerr<<"Start printing VaKmers"<<endl;
    time_t prev_time=time(NULL);
    seqspace.printAbsentKmers2(cout);
    now_time=time(NULL);
    cerr<<"Finish printing VaKmers in "<<(now_time-prev_time)<<" second(s)"<<endl;
    cerr<<"Done in "<<(now_time-start_time)<<" second(s)"<<endl;
    return 0;
    

}

#endif //__MODE4__

#ifdef __MODE4C__

int main(int argc,char **argv)
{

    string extra_settings="";

    #ifdef __FAST_IO__
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    extra_settings=" FastIO option=True";
    #endif //__FAST_IO__

    cerr<<"VaKation (\e[4mVa\e[0mcant \e[4mK\e[0m-mers Identific\e[4mation\e[0m) vers 0.1 mode 4 "<<extra_settings<<endl;
    


    if(argc<3){
        
        return printUsageAndExit(argv[0]);
    }
    
    int kmer=StringUtil::atoi(argv[1]);
    

    time_t start_time=time(NULL);

    StepwiseSeqSpace seqspace(kmer);
    
    
    
    for(int i=2;i<argc;i++){
        

        time_t prev_time=time(NULL);
        
        cerr<<"Processing fasta file "<<argv[i]<<endl;
        FastaFile fastaFile(argv[i]);
        while (true){

            char c=toupper(fastaFile.readNextNt());
            if(c=='\0'){
                break;
            }else if(c=='>'){
                cerr<<"Processing sequence "<<fastaFile.seqName<<endl;
                seqspace.see('N');
            }else{
                //cerr<<c<<endl;
                seqspace.see(c);
                //cerr<<c<<endl;
            }
        }

        time_t now_time=time(NULL);
        cerr<<"Finish "<<argv[i]<<" in "<<(now_time-prev_time)<<" second(s)"<<endl;
        
    }

    time_t now_time=time(NULL);
    cerr<<"Finish processing all fasta files in "<<(now_time-start_time)<<" second(s)"<<endl;

    
    cerr<<"Start counting VaKmers"<<endl;
    time_t prev_time=time(NULL);
    seqspace.countAbsentKmers2();
    now_time=time(NULL);
    cerr<<"Finish counting VaKmers in "<<(now_time-prev_time)<<" second(s)"<<endl;
    cerr<<"Done in "<<(now_time-start_time)<<" second(s)"<<endl;
    return 0;
    

}

#endif //__MODE4C__