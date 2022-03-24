#include <fstream>
#include <iostream>

using namespace std;

int main(int argc, const char ** argv){

    cerr<<"size_t = "<<((size_t)-1)<<endl;

    if(argc<4){
        cerr<<argv[0]<<" filename sizeInGBytes byteValueToWrite"<<endl;
        return 1;
    }

    

    uint64_t multiplier=1024*1024*1024;

    const char*filename=argv[1];
    unsigned int sizeInGBytes=atoi(argv[2]);
    uint64_t sizeInBytes=sizeInGBytes*multiplier;
    unsigned char byteValueToWrite=atoi(argv[3]);
 

    ofstream fout;
    fout.open(filename,ios::binary|ios::out);
    cerr<<"starting to write "<<sizeInBytes<<" bytes"<<endl;

    for(uint64_t i=0;i<sizeInBytes;i++){
        if(i%multiplier==1){
            cerr<<"Writting Byte "<<i<<endl;
        }
        fout.write((char*)&byteValueToWrite,1);
    }
    fout.close();
    
    return 0;
}