/* bigBedToBedPlus - Convert from bigBed to ascii bed format.. */

/* Copyright (C) 2014 The Regents of the University of California 
 * See kent/LICENSE or http://genome.ucsc.edu/license/ for licensing information. */

/* 
*  2022/4/27: Albert Cheng, PhD, ASU added new options to bigBedToBed to make it bigBedToBed Plus, supporting a list of chrom start end ranges as well as filtering criteria for columns; also output is changed to STDOUT
*/

#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"
#include "localmem.h"
#include "udc.h"
#include "bigBed.h"
#include "asParse.h"
#include "obscure.h"


char *clRangeListFile = NULL;
int maxItems = 0;

int *gFilterMins=NULL;
int *gFilterMaxs=NULL;
int gFilterMaxColNum1=0;
char * gInBB=NULL;
void usage()
/* Explain usage and exit. */
{
errAbort(
  "bigBedToBedPlus v1 - Convert from bigBed to ascii bed format and print to STDOUT with ranges and filters provided in rangeListFile.\n"
  "usage:\n"
  "   bigBedToBedPlus input.bb\n"
  "options:\n"
  "   -maxItems=N - if set, restrict output to first N items\n"
  "   -udcDir=/dir/to/cache - place to put cache for remote bigBed/bigWigs\n"
  "   -filter=filename - if set restrict output to given list of chr start end in the specific file\n" 
  "   The range file contains each row tab delimited chr start end for chrom ranges or $col min max for min-max inclusive filtering for col\n"
  "   For example extract items within chr1:1-50000 or chr2:200-100000 with column 17 in range from 40 to 60 inclusive and column 14 in range from 0 to 3 inclusive:\n"
  "   run bigBedToBedPlus -filter=filter.txt input.bb\n"
  "   filter.txt:\n"
  "   chr1\t1\t50000\n"
  "   chr2\t200\t100000\n"
  "   $17\t40\t60\n"
  "   $14\t0\t3\n"
  "   if intput.bb (by adding a line preceded by colon,i.e., :input.bb) included in filter.txt, run bigBedToBedPlus -filter=filter.txt\n"
  "   filter.txt:\n"
  "   :input.bb\n"
  "   chr1\t1\t50000\n"
  "   chr1\t1\t50000\n"
  "   chr2\t200\t100000\n"
  "   $17\t40\t60\n"
  "   $14\t0\t3\n"

  );
}



struct FilterColumnCriteria{
    int colNum1;
    int min;
    int max;
    struct FilterColumnCriteria* next;
};

struct FilterColumnCriteria* gFirstFilter=NULL;

void addFilter(int _colNum1, int _min, int _max){
    struct FilterColumnCriteria*  newFilter=(struct FilterColumnCriteria*)malloc(sizeof(struct FilterColumnCriteria));
    newFilter->next=gFirstFilter;
    newFilter->colNum1=_colNum1;
    newFilter->min=_min;
    newFilter->max=_max;
    gFirstFilter=newFilter;
    if(_colNum1>gFilterMaxColNum1){
        gFilterMaxColNum1=_colNum1;
    }

}


void printFilter(struct FilterColumnCriteria* filter){
    printf("filter col %d for [%d,%d] inclusive\n",filter->colNum1,filter->min,filter->max);
}

void freeFilters(){
    struct FilterColumnCriteria *x=gFirstFilter;

    while(x){
        struct FilterColumnCriteria*xToFree=x;
        x=x->next;
        //printf("freeing :");printFilter(xToFree);
        free(xToFree);
    }

    gFirstFilter=NULL;
}

#define SMALLEST_INT -10000

void convertFiltersToArray(){
    if(!gFirstFilter)
        return;

    gFilterMins=(int*)malloc(gFilterMaxColNum1*sizeof(int)+sizeof(int));
    gFilterMaxs=(int*)malloc(gFilterMaxColNum1*sizeof(int)+sizeof(int));
    for(int i=0;i<=gFilterMaxColNum1;i++){
        gFilterMins[i]=SMALLEST_INT;
        gFilterMaxs[i]=SMALLEST_INT;
    }
    for(struct FilterColumnCriteria*x=gFirstFilter;x!=NULL;x=x->next){
        gFilterMins[x->colNum1]=x->min;
        gFilterMaxs[x->colNum1]=x->max;
    }

    freeFilters();
    
}

void freeFilterArrays(){
    if(gFilterMins)
        free(gFilterMins);
    
    if(gFilterMaxs)
        free(gFilterMaxs);
    
}





struct StartEnd{
    int start;
    int end;
    struct StartEnd* next;
};


struct RangeNode {
    struct RangeNode *next;
    char* chrom;
    struct StartEnd* firstStartEnd;
};

struct RangeNode *gHeadNode;

struct RangeNode* findChrom( struct RangeNode*headNode, const char*chrom){
    if(!headNode)
        return NULL;

    for(struct RangeNode*cur=headNode;cur!=NULL;cur=cur->next){
        if(!strcmp(chrom,cur->chrom)){
            return cur;
        }
    }

    return NULL;
};

void addStartEndToChrom(struct RangeNode*targetNode, int start, int end){
    struct StartEnd *newStartEnd=(struct StartEnd*)malloc(sizeof(struct StartEnd));
    newStartEnd->start=start;
    newStartEnd->end=end;
    struct StartEnd *oldFirst=targetNode->firstStartEnd;
    newStartEnd->next=oldFirst;
    targetNode->firstStartEnd=newStartEnd;
}

struct RangeNode* addChrom(struct RangeNode*headNode,const char* chrom){
    struct RangeNode *newRangeNode=(struct RangeNode*)malloc(sizeof(struct RangeNode));
    newRangeNode->chrom=(char*)malloc(strlen(chrom)+1);
    strcpy(newRangeNode->chrom,chrom);
    newRangeNode->firstStartEnd=NULL;
    newRangeNode->next=headNode;
    return newRangeNode;
}

void addRange(const char*chrom, int start, int end){
    struct RangeNode*foundChrom=findChrom(gHeadNode,chrom);
    if(!foundChrom){
        gHeadNode=addChrom(gHeadNode,chrom);
        foundChrom=gHeadNode;
    }

    addStartEndToChrom(foundChrom,start,end);
}

void printRangeList(){
    for(struct RangeNode*x=gHeadNode;x!=NULL;x=x->next){
        printf("%s:\n",x->chrom);
        for(struct StartEnd*y=x->firstStartEnd;y!=NULL;y=y->next){
            printf("\t%s\t%d\t%d\n",x->chrom,y->start,y->end);
        }
    }
}

void freeRangeList(){
    struct RangeNode *x=gHeadNode;
    if(!gHeadNode)
        return;

    do{
        struct RangeNode*xToFree=x;
        x=x->next;
        
        struct StartEnd *y=xToFree->firstStartEnd;
        do{
            struct StartEnd*yToFree=y;
            y=y->next;
            //printf("free %s:%d-%d\n",xToFree->chrom,yToFree->start,yToFree->end);
            free(yToFree);
        }while(y);

        //printf("free %s\n",xToFree->chrom);
        free(xToFree->chrom);
        free(xToFree);
    }while(x);
}



static struct optionSpec options[] = {
   {"filter", OPTION_STRING},
   {"maxItems", OPTION_INT},
   {"udcDir", OPTION_STRING},
   {"header", OPTION_BOOLEAN},
   {NULL, 0},
};


void writeHeader(struct bbiFile *bbi, FILE *f)
/* output a header from the autoSql in the file */
{
char *asText = bigBedAutoSqlText(bbi);
if (asText == NULL)
    errAbort("bigBed files does not contain an autoSql schema");
struct asObject *asObj = asParseText(asText);
char sep = '#';
for (struct asColumn *asCol = asObj->columnList; asCol != NULL; asCol = asCol->next)
    {
    fputc(sep, f);
    fputs(asCol->name, f);
    sep = '\t';
    }
fputc('\n', f);
}

void repairStringAfterStrtok(char *str,char delimiter,int length){
    for(int i=0;i<length;i++){
        if(*str=='\0'){
            *str=delimiter;
        }
        str++;
    }
}



void bigBedToBed(char *inFile)
/* bigBedToBed - Convert from bigBed to ascii bed format.. */
{
struct bbiFile *bbi = bigBedFileOpen(inFile);

/*if (header)
    writeHeader(bbi, f);*/

struct bbiChromInfo *chrom, *chromList = bbiChromList(bbi);
int itemCount = 0;
for (chrom = chromList; chrom != NULL; chrom = chrom->next)
    {
    /*if (clChrom != NULL && !sameString(clChrom, chrom->name))
        continue;*/
        //printf("traversing %s\n",chrom->name);

        if(!clRangeListFile){
            //no range defined, add psuedorange
            addRange(chrom->name,0,0);
        }

        struct RangeNode* chromRangeRequest=findChrom(gHeadNode,chrom->name);
        if(!chromRangeRequest) //this chromosome is not requested.
            continue;

        //printf("found\n");
        char *chromName = chrom->name;



        for(struct StartEnd*startend=chromRangeRequest->firstStartEnd;startend!=NULL;startend=startend->next){ 
            //printf("request %s:%d-%d\n",chromRangeRequest->chrom,startend->start,startend->end);
            int start = 0, end = chrom->size;
            if (startend->start > 0)
                start = startend->start;
            if (startend->end > 0)
                end = startend->end;
            int itemsLeft = 0;	// Zero actually means no limit.... 
            if (maxItems != 0)
                {
            itemsLeft = maxItems - itemCount;
            if (itemsLeft <= 0)
                break;
            }
            struct lm *lm = lmInit(0);
            struct bigBedInterval *interval, *intervalList = bigBedIntervalQuery(bbi, chromName, 
                start, end, itemsLeft, lm);
            for (interval = intervalList; interval != NULL; interval = interval->next)
            {
            
                
                char *rest = interval->rest;

                
                if(rest!=NULL){

                    if(gFilterMins){ //has some filters

                        int accepted=1; //assume all criteria OK

                        int restLength=strlen(rest);
                        char *token;

                        
                        token = strtok(rest, "\t");
                        int level=4;

                        while( token != NULL ) {
                            //printf( "\t%d:%s\n", level,token );
                            if(level>gFilterMaxColNum1){
                                break;
                            }
                            
                            if(gFilterMins[level]!=SMALLEST_INT && atoi(token)<gFilterMins[level]){
                                accepted=0; //failed min requirement
                                break;
                            }

                            if(gFilterMaxs[level]!=SMALLEST_INT && atoi(token)>gFilterMaxs[level]){
                                accepted=0;
                                break;
                            }


                            token = strtok(NULL, "\t");

                            level++;
                        }

                        if(accepted){
                            repairStringAfterStrtok(rest,'\t',restLength);
                            printf("%s\t%u\t%u", chromName, interval->start, interval->end);
                            printf("\t%s\n", rest);
                        }

                    }
                    else{
                        printf("%s\t%u\t%u", chromName, interval->start, interval->end);
                        printf("\t%s\n", rest);
                    }
                }
                else{
                    printf("%s\t%u\t%u", chromName, interval->start, interval->end);
                    printf( "\n");
                }



                    
            }
            lmCleanup(&lm);

        }
        
    }
bbiChromInfoFreeList(&chromList);
//carefulClose(&f);
bbiFileClose(&bbi);
}

void readChromRangeFile(const char*filename){

    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;

    fp = fopen(filename, "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);

    while ((read = getline(&line, &len, fp)) != -1) {
        //printf("Retrieved line of length %zu:\n", read);
        //printf("%s", line);
        char *token;
        
        char *chrom=NULL;
        int start=0;
        int end=0;
        int colNum1=0;
        int _min=0;
        int _max=0;

        int level=0;
        token = strtok(line, "\t");
        int readFilter=0;

        if(line[0]=='$'){
            //this is a filter
            readFilter=1;

        }

        if(line[0]==':'){
            gInBB=(char*)malloc(strlen(line));
            strcpy(gInBB,line+1);
            gInBB[strlen(line)-2]='\0';
            //printf("get file[%s]",gInBB);
            continue;
        }

        while( token != NULL ) {
            //printf( "\t%s\n", token );
            switch(level){
                case 0:
                    if(readFilter)
                        colNum1=atoi(token+1); //skip "$"
                    else
                        chrom=token;
                break;
                case 1:
                    if(readFilter)
                        _min=atoi(token);
                    else
                        start=atoi(token);
                break;
                case 2:
                    if(readFilter)
                        _max=atoi(token);
                    else
                        end=atoi(token);
                break;
            }

            token = strtok(NULL, "\t");
            level++;
        }     

        if(chrom){
            //printf("add range:%s:%d-%d\n",chrom,start,end);
            addRange(chrom,start,end);

        }else{
            //filter
            addFilter(colNum1,_min,_max);
        }

    }

    fclose(fp);
    if (line)
        free(line);

    fclose(fp);
}




int main(int argc, char *argv[])
/* Process command line. */
{




optionInit(&argc, argv, options);


clRangeListFile = optionVal("filter", clRangeListFile);
if(clRangeListFile){
    readChromRangeFile(clRangeListFile);
    //printRangeList();
    convertFiltersToArray();
}




maxItems = optionInt("maxItems", maxItems);
udcSetDefaultDir(optionVal("udcDir", udcDefaultDir()));

if(gInBB){
    bigBedToBed(gInBB);
}else{
    if (argc != 2)
        usage();
    bigBedToBed(argv[1]);
}

if (verboseLevel() > 1)
    printVmPeak();


freeRangeList();
//freeFilters();
freeFilterArrays();

if(gInBB)
    free(gInBB);

return 0;
}
