#ifndef __BITSTRING_H
#define __BITSTRING_H

#include <iostream>
#include <stdint.h>
#include <fstream>

//#if __has_include(<zlib.h>)
#include <zlib.h>
#include <errno.h>
#include <string.h>
#define GZIP_BUFFER_LENGTH 4095
#define GZIP_WRITE_BUFFER_LENGTH 4095
//#endif

#include <sys/mman.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>


typedef unsigned char Byte;

using namespace std;


class BitString{
public:

    uint64_t numBytes;
    uint64_t numBits;
    struct stat statbuf;
    Byte* _data;
    
    #ifdef GZIP_BUFFER_LENGTH
    bool gzBitString(const string& filename,uint64_t _numBits)
    {
        numBits=_numBits;
        numBytes=numBits/8;
        if(numBits%8>0){
            numBytes++;
        }
        
        gzFile file;
        file=gzopen(filename.c_str(),"r");
        if(!file){
            cerr<<"gzopen of "<<filename<<" failed "<<strerror(errno)<<endl;
            return false;
        }else{
            cerr<<"gzopen of "<<filename<<" succeeded. Continue to load"<<endl;
        }

        _data=new Byte[numBytes];

        uint64_t byteRemaining=numBytes;

        for(uint64_t i=0;i<numBytes;i+=GZIP_BUFFER_LENGTH){
            int bytes_read=gzread(file,_data+i,(byteRemaining<GZIP_BUFFER_LENGTH)?byteRemaining:GZIP_BUFFER_LENGTH);
            byteRemaining-=bytes_read;
            if(!bytes_read){
                if(gzeof(file)){
                    break;
                }else{
                    const char*error_string;
                    int err;
                    error_string=gzerror(file,&err);
                    if(err){
                        cerr<<"Error: "<<error_string<<endl;
                        delete[] _data;
                        gzclose(file);
                        return false;
                    }
                }
            }
        }

        gzclose(file);
        return true;

    }
    #endif //GZIP_BUFFER_LENGTH

    bool operator == (const BitString& _right) const {
        //if(this->numBits!=_right.numBits)
        //    return false;
        
        
        if(this->numBytes!=_right.numBytes)
            return false;

        Byte* leftI=this->_data;
        Byte* rightI=_right._data;
        for(uint64_t i=0;i<this->numBytes;i++){
            if(*(leftI++)!=*(rightI++))
                return false;
        }

        return true;
    }

    BitString(const string& filename,uint64_t _numBits=0):_data(NULL){
        #ifndef GZIP_BUFFER_LENGTH
        _numBits=0;
        #endif 

        if(_numBits>0){
            gzBitString(filename,_numBits);
        }else{
            ifstream ifil(filename.c_str(),ios::binary|ios::in);
            ifil.seekg(0,ios::end);
            numBytes=ifil.tellg();
            numBits=numBytes*8;
            //ifil.close();
            _data=new Byte[numBytes];
            ifil.clear();
            ifil.seekg(0);
            ifil.read((char*)_data,numBytes);
            ifil.close();
        }
    }
    
    const BitString& operator-=(const BitString& _right){
        for(uint64_t i=0;i<numBytes;i++){
            _data[i]&=(~_right._data[i]);
        }
        
        return *this;
    }
    
    const BitString& operator&=(const BitString& _right){
        for(uint64_t i=0;i<numBytes;i++){
            _data[i]&=_right._data[i];
        }
        
        return *this;
    }
    
    const BitString& operator|=(const BitString& _right){
        for(uint64_t i=0;i<numBytes;i++){
            _data[i]|=_right._data[i];
        }
        
        return *this;
    }
    
    
    void init(uint64_t _numBits,Byte valuePerByte=0)
    {
        numBits=_numBits;
        numBytes=numBits/8;
        if(numBits%8>0){
            numBytes++;
        }
        
        cerr<<"trying map"<<endl;

        int fd = open("/Users/albert/mmap.bbb",O_RDWR);
        if(fd<0){
            cerr<<"cannot open ~/mmap.bbb"<<endl;
            exit(1);
        }

        
        int err = fstat(fd, &statbuf);
        if(err < 0){
            cerr<<"2 cannot open ~/mmap.bbb"<<endl;
             exit(2);
        }else{
            cerr<<"statbuf.st_size="<<statbuf.st_size<<endl;
        }

        //_data=new Byte[numBytes];
        _data=(Byte*)mmap(NULL,statbuf.st_size,PROT_READ|PROT_WRITE,MAP_SHARED,fd,0);
        if(_data==MAP_FAILED){
            cerr<<"mmap failed"<<endl;
        }else{
            cerr<<"mmap success"<<endl;
            
            /*int err=munmap(_data,statbuf.st_size);
            if(err!=0){
                cerr<<"unmap failed"<<endl;
            }*/
            
        }

        //exit(1);

        for(uint64_t i=0;i<numBytes;i++){
            if(i%10000000==1){
                cerr<<"set zero "<<i<<endl;
            }
            _data[i]=valuePerByte;
        }

        cerr<<"_data init success"<<endl;
    }
    
    BitString(uint64_t _numBits)
    {

        init(_numBits);
    }
    
    BitString():numBytes(0),numBits(0),_data(NULL){
    }
    
    
    ~BitString(){
        /*if(_data)
            delete[] _data;*/
        if(_data!=MAP_FAILED){
            int err=munmap(_data,statbuf.st_size);
            if(err!=0){
                cerr<<"unmap failed"<<endl;
            }
        }
        
    }
    inline void getByteBitIndex(uint64_t index,uint64_t& byteIndex,uint64_t& bitIndex){
        byteIndex=index/8;
        bitIndex=index%8;
    }
    void setBit(uint64_t index,bool value){
        uint64_t byteIndex;
        uint64_t bitIndex;
        getByteBitIndex(index,byteIndex,bitIndex);
        if(value){
            _data[byteIndex]|=(1<<bitIndex);
        }else{
            _data[byteIndex]&=(~(1<<bitIndex));
        }
    }
    
    bool getBit(uint64_t index){
        uint64_t byteIndex;
        uint64_t bitIndex;
        getByteBitIndex(index,byteIndex,bitIndex);
        return _data[byteIndex]&(1<<bitIndex);
    }


    
    void print(ostream& os,bool printReturn=false){
        for(uint64_t i=0;i<numBytes;i++){
            for (uint64_t j=0;j<8;j++){
                os<<bool(_data[i]&(1<<j));
            }
            
            os<<" ";
        }
        
        if(printReturn)
            os<<endl;
    }
    
    void writeToFile(const string& filename)
    {
        ofstream ofil(filename.c_str(),ios::binary|ios::out|ios::trunc);
        ofil.write((const char*)_data,numBytes);
        ofil.close();
    }

    bool writeToFileGZ(const string&filename)
    {
        gzFile file;
        file=gzopen(filename.c_str(),"w");
        if(!file){
            cerr<<"gzopen of "<<filename<<" failed "<<strerror(errno)<<endl;
            return false;
        }else{
            cerr<<"gzopen of "<<filename<<" succeeded. Proceed to writing"<<endl;
        }

        uint64_t numBytesRemaining=numBytes;

        for(uint64_t i=0;i<numBytes;i+=GZIP_WRITE_BUFFER_LENGTH){

            uint64_t byte_to_write=(numBytesRemaining<GZIP_WRITE_BUFFER_LENGTH)?numBytesRemaining:GZIP_WRITE_BUFFER_LENGTH;

            int byte_written=gzwrite(file,_data+i,byte_to_write);

            if(byte_written<byte_to_write){ //something is wrong!
                const char*error_string;
                int err;
                error_string=gzerror(file,&err);
                if(err){
                    cerr<<"Error: "<<error_string<<endl;
                    gzclose(file);
                    return false;
                    
                }else{
                    cerr<<"Unknown write error"<<endl;
                    gzclose(file);
                    return false;
                }
            }
        }        

        gzclose(file);
        return true;
    }
    
   	string getBitsAsString(uint64_t start0, uint64_t end1)
   	{
   		string s;
   		for(uint64_t i=start0;i<end1;i++){
   			s+=getBit(i)?"1":"0";	
   		}
   		return s;
   	}	
    void setBitsOnRange(uint64_t start0, uint64_t end1, bool value)
    {
        for(uint64_t i=start0;i<end1;i++){
            setBit(i,value);
        }
    }
    
};

class FileBitString{
private:
    
    string filename;
    uint64_t numBytes;
    uint64_t numBits;
    fstream *fil;
public:

    FileBitString(const string& _filename,uint64_t _numBits):numBits(_numBits),filename(_filename){
        numBytes=numBits/8;
        if(numBits%8>0){
            numBytes++;
        }
        
        //cerr<<numBytes<<"\t"<<numBits<<endl;
        
        fil=new fstream(_filename.c_str(),ios::in|ios::out|ios::binary);
        if(!fil->good())
        {
            //cerr<<"new file"<<endl;
            fil->close();
            delete fil;
            fil=new fstream(_filename.c_str(),ios::in|ios::out|ios::binary|ios::trunc);
        }
        fil->seekp(0,ios::end);
        
        uint64_t curFileSize=fil->tellp();
        while (fil->tellp()<numBytes) {
            fil->put(0);
            //cerr<<"put"<<fil->tellp()<<endl;
        }
        
        fil->clear();
        
        
    }

inline void getByteBitIndex(uint64_t index,uint64_t& byteIndex,uint64_t& bitIndex){
    byteIndex=index/8;
    bitIndex=index%8;
}
void setBit(uint64_t index,bool value){
    uint64_t byteIndex;
    uint64_t bitIndex;
    getByteBitIndex(index,byteIndex,bitIndex);
    
    fil->seekg(byteIndex);

    char _data=fil->get();
    
    if(value){
        _data|=(1<<bitIndex);
    }else{
        _data&=(~(1<<bitIndex));
    }
    fil->seekp(byteIndex);
    fil->put(_data);
}

bool getBit(uint64_t index){
    uint64_t byteIndex;
    uint64_t bitIndex;
    getByteBitIndex(index,byteIndex,bitIndex);
    
    fil->seekg(byteIndex);
    char _data=fil->get();
    return _data&(1<<bitIndex);
}
    
    void print(ostream& os,bool printReturn=false){
        fil->seekg(0);
        while(fil->good()){
            Byte _data=fil->get();
            if(fil->good()){
            for(int j=0;j<8;j++){
                os<<bool(_data&(1<<j));
                
            }
                os<<" ";
            }
        }
        fil->clear();
        if(printReturn)
            os<<endl;
    }


    ~FileBitString(){
        fil->close();
        delete fil;
    }

};

#endif