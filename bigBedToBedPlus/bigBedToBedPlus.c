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

int gWriteQueryToName=FALSE;
char *clRangeListFile = NULL;
char *clFilterFile=NULL;
int maxItems = 0;
boolean header=0;
boolean gPrintComment=0;
int *gFilterMins=NULL;
int *gFilterMaxs=NULL;
int gFilterMaxColNum1=0;
char * gInBB=NULL;

void replaceCharInString(char *target,char from,char to){
    while((*target)!='\0'){
        
        if((*target)==from){
            //fprintf(stderr,"%c->%c",from,to);
            (*target)=to;
        }
        target++;
    }
}

void usage()
/* Explain usage and exit. */
{
errAbort(
  "bigBedToBedPlus v1.2 - Convert from bigBed to ascii bed format and print to STDOUT with ranges and filters provided in rangeListFile.\n"
  "usage:\n"
  "   bigBedToBedPlus input.bb [inlineFiltersOrSelections(/separated) ... ]\n"
  "options:\n"
  "   -addQueryRangeToName - if set, add query name to name field of output\n"
  "   -maxItems=N - if set, restrict output to first N items\n"
  "   -udcDir=/dir/to/cache - place to put cache for remote bigBed/bigWigs\n"
  "   -ranges=filename - if set restrict output to given list of chr start end in the specific file and other filter or sorting selection critiera\n" 
  "   -filters=filename - same effect as ranges\n"
  "   The range file contains each row tab delimited chr start end for chrom ranges or $col min max for min-max inclusive filtering for col\n"
  "   For example extract items within chr1:1-50000 or chr2:200-100000 with column 17 in range from 40 to 60 inclusive and column 14 in range from 0 to 3 inclusive:\n"
  "   run bigBedToBedPlus -filters=filter.txt input.bb\n"
  "   filter.txt:\n"
  "   chr1\t1\t50000\n"
  "   chr2\t200\t100000\n"
  "   $17\t40\t60\n"
  "   $14\t0\t3\n"
  "   include input.bb in filter.txt by adding a line preceded by colon,i.e., :input.bb, run bigBedToBedPlus -filter=filter.txt\n"
  "   filter.txt:\n"
  "   :input.bb\n"
  "   chr1\t1\t50000\n"
  "   chr1\t1\t50000\n"
  "   chr2\t200\t100000\n"
  "   $17\t40\t60\n"
  "   $14\t0\t3\n"
  "   include !MIN min, !MAX max and !BEST numItems to select the best numItems items sorting with ascending(min) or descending(max) on the specificed columns in the order as they appear in the filter.txt.\n"
  "   For example, 4 best items sorting first ascending with col 14, then descending with col 17\n"
  "   filter.txt:\n"
  "   :input.bb\n"
  "   chr1\t1\t50000\n"
  "   chr1\t1\t50000\n"
  "   chr2\t200\t100000\n"
  "   $17\t40\t60\n"
  "   $14\t0\t3\n"
  "   !BEST\t4\n"
  "   !MIN\t14\n"
  "   !MAX\t17\n"  
  );
}



typedef struct LineData {
    int *colValues;
    int numColValues;
    char *line;
    int lineBufferLength;
} LineData;

LineData* createLineData(int _numColValues,const char*_line){
    LineData* newLineData=(LineData*)malloc(sizeof(LineData));
    newLineData->numColValues=_numColValues;
    newLineData->colValues=(int*)malloc(sizeof(int)*_numColValues);
    newLineData->line=(char*)malloc(strlen(_line)+1);
    newLineData->lineBufferLength=strlen(_line)+1;
    strcpy(newLineData->line,_line);
    return newLineData;
}

LineData* duplicateLineData(const LineData*src){
    LineData* newLineData=createLineData(src->numColValues,src->line);
    for(int i=0;i<src->numColValues;i++){
        newLineData->colValues[i]=src->colValues[i];
    }
    return newLineData;
}

void freeLineData(LineData *target){
    free(target->colValues);
    free(target->line);
    free(target);
}

int getColValueFromLineData(LineData* src,int _col1){
    return src->colValues[_col1-1];
}

void setLineDataColValue(LineData* target,int _col1, int _setvalue){
    target->colValues[_col1-1]=_setvalue;
}

void setLineDataLineString(LineData* target,const char*_str){
    int incomingLength=strlen(_str);
    if(incomingLength+1>target->lineBufferLength){
        free(target->line);
        target->line=(char*)malloc(incomingLength+1);
    }
    strcpy(target->line,_str);
}


struct FilterColumnCriteria{
    char colString[100];
    int colNum1;
    int min;
    int max;
    struct FilterColumnCriteria* next;
};

struct FilterColumnCriteria* gFirstFilter=NULL;





void addFilter(char *_colString,int _colNum1, int _min, int _max){
    struct FilterColumnCriteria*  newFilter=(struct FilterColumnCriteria*)malloc(sizeof(struct FilterColumnCriteria));
    strcpy(newFilter->colString,_colString);
    newFilter->next=gFirstFilter;
    newFilter->colNum1=_colNum1;
    newFilter->min=_min;
    newFilter->max=_max;
    gFirstFilter=newFilter;

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



void freeFilterArrays(){
    if(gFilterMins)
        free(gFilterMins);
    
    if(gFilterMaxs)
        free(gFilterMaxs);
    
}

int gBestNum=0;

#define SORT_MAX 1
#define SORT_MIN 2

typedef struct SortCriterion {
    char colString[100];
    int type;
    int colNum1;
    struct SortCriterion* next;

} SortCriterion;

SortCriterion *gLastSortCriterion=NULL;
SortCriterion *gSortCriteriaArray=NULL;
int gNumSortCriteria=0;
int gSortCrtieriaMaxColNum1=0;
int *gSortWatchList=NULL;

void addSortCriterion(int _type,char* _colString,int _colNum1 ){
    SortCriterion *newSortCriterion=(SortCriterion*)malloc(sizeof(SortCriterion));
    strcpy(newSortCriterion->colString,_colString);
    newSortCriterion->type=_type;
    newSortCriterion->colNum1=_colNum1;
    /*if(_colNum1>gSortCrtieriaMaxColNum1){
        gSortCrtieriaMaxColNum1=_colNum1;
    }*/
    newSortCriterion->next=gLastSortCriterion;
    gLastSortCriterion=newSortCriterion;
    gNumSortCriteria++;
}






void printSortCriteria(){
    if(gSortCriteriaArray){
        for(int i=0;i<gNumSortCriteria;i++){
            printf("Sort Criterion %d: %s on Column %d\n",i+1,(gSortCriteriaArray[i].type==SORT_MAX)?"MAX":"MIN",gSortCriteriaArray[i].colNum1);
        }
    }
}

void freeSortCriteria(){
    if(gSortCriteriaArray){
        free(gSortCriteriaArray);
    }
    if(gSortWatchList){
        free(gSortWatchList);
    }
}


int better(LineData* thisData,LineData* thatData){  //is thisData better than thatData according to the sort Criteria?
    for(int i=0;i<gNumSortCriteria;i++){
        int thisCol1=gSortCriteriaArray[i].colNum1;
        int thisValue=getColValueFromLineData(thisData,thisCol1);
        int thatValue=getColValueFromLineData(thatData,thisCol1);

        switch(gSortCriteriaArray[i].type){
            case SORT_MIN:
            if(thisValue<thatValue)
                return 1;
            else if(thisValue>thatValue)
                return 0;
            break;
            case SORT_MAX:
            if(thisValue>thatValue)
                return 1;
            else if(thisValue<thatValue)
                return 0;
            break;
        }
    }

    return 0;
}


void makeSureLineDataStringHaveAtLeast(LineData*_data,int _numchars){
    if(_data->lineBufferLength<_numchars){
        _data->lineBufferLength=_numchars;
        if(_data->line){
            free(_data->line);
        }
        _data->line=(char*)malloc(_numchars);
    }
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

struct RangeNode *gHeadNode=NULL;

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
    //printf("add range %s:%d-%d\n",chrom,start,end);
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


LineData** gBests=NULL;
int gNumFilled=0;


void updateBests(LineData* thisData){

    //printf("insert %d %d ",getColValueFromLineData(thisData,14),getColValueFromLineData(thisData,17));

    int betterThan=-1;
    for(int i=0;i<gNumFilled;i++){
        if(better(thisData,gBests[i])){
            //printf("bt [%d]:%d %d ,",i,getColValueFromLineData(gBests[i],14),getColValueFromLineData(gBests[i],17));
            //thisData better than bests[i]
            betterThan=i;
        }else{
            break; //not better than this, no need to go forward
        }
    }
    
    /*if(betterThan>=0)
        printf("better than [%d]:%d %d\n",betterThan,getColValueFromLineData(gBests[betterThan],14),getColValueFromLineData(gBests[betterThan],17));
    else
        printf("better than nothing\n");*/

    if(gNumFilled<gBestNum){
        gNumFilled++;
        for(int i=gNumFilled-1;i>betterThan+1;i--){
            gBests[i]=gBests[i-1]; //move up
        }
        gBests[betterThan+1]=duplicateLineData(thisData);
    }
    else{
        //numFilled==gBestNum
        if(betterThan>=0){
            //thisData better than the least best
            freeLineData(gBests[0]);
            for(int i=0;i<betterThan;i++){
                gBests[i]=gBests[i+1];
            }
            gBests[betterThan]=duplicateLineData(thisData);
        }

    }
    
    /*for(int i=gNumFilled-1;i>=0;i--){
        printf("%d %d\n",getColValueFromLineData(gBests[i],14),getColValueFromLineData(gBests[i],17));
        
    }
    printf("------\n");*/

}



static struct optionSpec options[] = {
   {"filters", OPTION_STRING},
   {"ranges", OPTION_STRING},
   {"maxItems", OPTION_INT},
   {"udcDir", OPTION_STRING},
   {"addQueryRangeToName", OPTION_BOOLEAN},
   {"header", OPTION_BOOLEAN},
   {"printComments", OPTION_BOOLEAN},
   {NULL, 0},
};

char *gAsText;
struct asObject *gAsObj;
int gNeedAsObject=0;

void writeHeader(struct bbiFile *bbi)
/* output a header from the autoSql in the file */
{

char* sep=(char*)malloc(2);
strcpy(sep,"#");

    for (struct asColumn *asCol = gAsObj->columnList; asCol != NULL; asCol = asCol->next)
    {
        printf("%s%s",sep,asCol->name);
        strcpy(sep,"\t");
    }
    printf("\n");

free(sep);

}

int colStringToColNum1(char *colString){
    int col1=1;
    for (struct asColumn *asCol = gAsObj->columnList; asCol != NULL; asCol = asCol->next)
    {
        
        if(!strcmp(asCol->name,colString)){
            if(gPrintComment) printf("#%s->%d\n",colString,col1);
            return col1;
        }
        col1++;
    }
    return -1;

}


void convertFiltersToArray(){
    if(!gFirstFilter)
        return;

    for(struct FilterColumnCriteria*x=gFirstFilter;x!=NULL;x=x->next){
        if(x->colNum1==SMALLEST_INT){
            x->colNum1=colStringToColNum1(x->colString);
            if(x->colNum1<1){
                
                errAbort("failed to find column with name %s",x->colString);
                
            }
        }

        if(x->colNum1>gFilterMaxColNum1){
            gFilterMaxColNum1=x->colNum1;
        }
    }    



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


void transferSortCriteriaToArray(){
    if(!gLastSortCriterion){
        return;
    }
    gSortWatchList=(int*)malloc(gSortCrtieriaMaxColNum1*sizeof(int)+sizeof(int));
    for(int i=0;i<=gSortCrtieriaMaxColNum1;i++){
        gSortWatchList[i]=0;
    }
    gSortCriteriaArray=(SortCriterion*)malloc(gNumSortCriteria*sizeof(SortCriterion));
    SortCriterion*x=gLastSortCriterion;
    for(int i=gNumSortCriteria-1;i>=0;i--){
        SortCriterion*xToFree=x;
        x=x->next;

        //xToFree convert colString -> colNum1 TODO!!
        if(xToFree->colNum1==SMALLEST_INT){
            xToFree->colNum1=colStringToColNum1(xToFree->colString);
            if(xToFree->colNum1<1){
                

                errAbort("failed to find column with name %s",xToFree->colString);
                
            }
        }

        if(xToFree->colNum1>gSortCrtieriaMaxColNum1){
            gSortCrtieriaMaxColNum1=xToFree->colNum1;
        }
        gSortCriteriaArray[i].type=xToFree->type;
        gSortCriteriaArray[i].colNum1=xToFree->colNum1;
        gSortWatchList[xToFree->colNum1]=1;
        free(xToFree);
    }
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

if (header || gNeedAsObject){

    gAsText= bigBedAutoSqlText(bbi);
    if (gAsText == NULL)
        errAbort("bigBed files does not contain an autoSql schema");

    gAsObj = asParseText(gAsText);


}

//remember header before this.
convertFiltersToArray();
transferSortCriteriaToArray();


//create buffer for bests

LineData *thisData=NULL;

if(gBestNum>0){
    gBests=(LineData**)malloc(sizeof(LineData*)*gBestNum);
    thisData=createLineData(gSortCrtieriaMaxColNum1,"NNN");
}

char queryRange[100];
queryRange[0]='\0';

if(header)
    writeHeader(bbi);

struct bbiChromInfo *chrom, *chromList = bbiChromList(bbi);
int itemCount = 0;



for (chrom = chromList; chrom != NULL; chrom = chrom->next)
    {
    /*if (clChrom != NULL && !sameString(clChrom, chrom->name))
        continue;*/
        //printf("traversing %s\n",chrom->name);
        

        struct RangeNode* chromRangeRequest;

        if(gHeadNode) // ranges present
        {
            //printf("has head node %s",gHeadNode->chrom);
            chromRangeRequest=findChrom(gHeadNode,chrom->name);
            if(!chromRangeRequest) //this chromosome is not requested.
                continue;
        }
        else{ //no ranges, setup a fake one
            //printf("set up fake rarnge for %s\n",chrom->name);
            chromRangeRequest=addChrom(NULL,chrom->name);
            addStartEndToChrom(chromRangeRequest,0,0);
        }

        //printf("found\n");
        char *chromName = chrom->name;



        for(struct StartEnd*startend=chromRangeRequest->firstStartEnd;startend!=NULL;startend=startend->next){ 
            //printf("request %s:%d-%d\n",chromRangeRequest->chrom,startend->start,startend->end);
            if(gPrintComment) printf("#range=%s:%d-%d\n",chromRangeRequest->chrom,startend->start,startend->end);

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


            gNumFilled=0;

            if(gWriteQueryToName){
                sprintf(queryRange,"%s:%d-%d/",chromName,start,end);
            }

            for (interval = intervalList; interval != NULL; interval = interval->next)
            {
            
                
                char *rest = interval->rest;

                
                if(rest!=NULL){

                    if(gFilterMins || gBests ){ //has some filters or select bests

                        int accepted=1; //assume all criteria OK

                        int restLength=strlen(rest);
                        char *token;

                        
                        token = strtok(rest, "\t");
                        int level=4;

                        while( token != NULL ) {
                            //printf( "\t%d:%s\n", level,token );
                            /*if(level>gFilterMaxColNum1){
                                break;
                            }*/
                            
                            if(gFilterMins && level<=gFilterMaxColNum1 && gFilterMins[level]!=SMALLEST_INT && atoi(token)<gFilterMins[level]){
                                accepted=0; //failed min requirement
                                break;
                            }

                            if(gFilterMaxs && level<=gFilterMaxColNum1 && gFilterMaxs[level]!=SMALLEST_INT && atoi(token)>gFilterMaxs[level]){
                                accepted=0;
                                break;
                            }

                            if(thisData && level<=gSortCrtieriaMaxColNum1 && gSortWatchList[level]){
                                setLineDataColValue(thisData,level,atoi(token));
                            }

                            token = strtok(NULL, "\t");

                            level++;
                        }

                        if(accepted){
                            repairStringAfterStrtok(rest,'\t',restLength);

                            if(gBestNum>0){
                                makeSureLineDataStringHaveAtLeast(thisData,strlen(rest)+1000);
                                sprintf(thisData->line,"%s\t%u\t%u\t%s%s\n", chromName, interval->start, interval->end, queryRange, rest);

                                updateBests(thisData);

                            }else{
                                printf("%s\t%u\t%u", chromName, interval->start, interval->end);
                                printf("\t%s%s\n", queryRange,rest);
                            }
                        }

                    }
                    else{
                        printf("%s\t%u\t%u", chromName, interval->start, interval->end);
                        printf("\t%s%s\n", queryRange, rest);
                    }
                }
                else{
                    printf("%s\t%u\t%u", chromName, interval->start, interval->end);
                    printf( "\n");
                }



                    
            }
            lmCleanup(&lm);


            for(int i=gNumFilled-1;i>=0;i--){
                printf("%s",gBests[i]->line);
                freeLineData(gBests[i]);
            }

        }

        if(!gHeadNode){
            free(chromRangeRequest->firstStartEnd);
            free(chromRangeRequest->chrom);
            free(chromRangeRequest);
        }
        
    }
bbiChromInfoFreeList(&chromList);
//carefulClose(&f);
bbiFileClose(&bbi);

if(thisData){
    freeLineData(thisData);
}

}

#define READ_CHROM_RANGE 0
#define READ_FILTER 1
#define READ_MIN 2
#define READ_MAX 3
#define READ_NUMBEST 4



void parseFilterLine(char *line){
        char colStringBuff[100];
        colStringBuff[0]='\0';

        char *token;
        
        char *chrom=NULL;
        int start=0;
        int end=0;
        int colNum1=0;
        int _min=0;
        int _max=0;

        int level=0;
        
        int readType=READ_CHROM_RANGE;

        //replaceCharInString(line,'\r','\0');
        replaceCharInString(line,'\n',' ');
        if(strlen(line)<2)
            return;

        if(line[0]=='#'){
            //comment
            if(gPrintComment && strlen(line)>1 && line[1]=='#'){
                //transfer comment into bed
                printf("%s\n",line);
            }
            return;
        }

        if(line[0]=='$'){
            //this is a filter
            readType=READ_FILTER;
            if(gPrintComment) printf("#%s\n",line);

        }

        if(line[0]==':'){
            gInBB=(char*)malloc(strlen(line));
            strcpy(gInBB,line+1);
            gInBB[strlen(line)-2]='\0';
            //printf("get file[%s]",gInBB);
            return;
        }

        if(!strncmp(line,"!MIN",4)){
            readType=READ_MIN;
            if(gPrintComment) printf("#%s\n",line);
        }
        if(!strncmp(line,"!MAX",4)){
            readType=READ_MAX;
            if(gPrintComment) printf("#%s\n",line);
        }
        if(!strncmp(line,"!BEST",5)){
            readType=READ_NUMBEST;
            if(gPrintComment) printf("#%s\n",line);
        }

        token = strtok(line, "\t");

        while( token != NULL ) {
            //printf( "\t%s\n", token );
            switch(level){
                case 0:
                    if(readType==READ_CHROM_RANGE){
                        chrom=token;
                    }
                    else if(readType==READ_FILTER){
                        char*colString=token+1; //skip "$"
                        if(colString[0]<'0' || colString[0]>'9') //not a number
                        {
                            gNeedAsObject=1;
                            colNum1=SMALLEST_INT;
                            strcpy(colStringBuff,colString);
                        }
                        else{
                            colNum1=atoi(colString);
                            colStringBuff[0]='\0';
                        }
                    }
                        
                break;
                case 1:
                    if(readType==READ_CHROM_RANGE){
                        start=atoi(token);
                    }else if(readType==READ_FILTER){
                        _min=atoi(token);
                    }else if(readType==READ_MIN){
                        if(token[0]<'0' || token[0]>'9') //not a number
                        {
                            gNeedAsObject=1;
                            addSortCriterion(SORT_MIN,token,SMALLEST_INT);
                        }else{
                            addSortCriterion(SORT_MIN,"",atoi(token));
                        }
                    }else if(readType==READ_MAX){
                        if(token[0]<'0' || token[0]>'9') //not a number
                        {
                            gNeedAsObject=1;
                            addSortCriterion(SORT_MAX,token,SMALLEST_INT);
                        }else{                        
                            addSortCriterion(SORT_MAX,"",atoi(token));
                        }
                    }else if(readType==READ_NUMBEST){
                        gBestNum=atoi(token);
                    }
                    
                        
                break;
                case 2:
                    if(readType==READ_CHROM_RANGE){
                        end=atoi(token);
                    }else if(readType==READ_FILTER){
                        _max=atoi(token);
                    }
                    
                        
                break;
            }

            token = strtok(NULL, "\t");
            level++;
        }     

        if(readType==READ_CHROM_RANGE){
            //printf("add range:%s\n",line);
            addRange(chrom,start,end);

        }else if(readType==READ_FILTER){
            //filter
            addFilter(colStringBuff,colNum1,_min,_max);
        }
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
        parseFilterLine(line);


    }

    fclose(fp);
    if (line)
        free(line);


}

void clean_up(){
    freeSortCriteria();
    freeRangeList();
    freeFilterArrays();

    if(gBests){
        free(gBests);
    }

    if(gInBB)
        free(gInBB);
    
}



int main(int argc, char *argv[])
/* Process command line. */
{


gHeadNode=NULL;

optionInit(&argc, argv, options);


clRangeListFile = optionVal("ranges", clRangeListFile);
clFilterFile = optionVal("filters", clFilterFile);
maxItems = optionInt("maxItems", maxItems);
udcSetDefaultDir(optionVal("udcDir", udcDefaultDir()));
gWriteQueryToName = optionExists("addQueryRangeToName");
header = optionExists("header");
gPrintComment = optionExists("printComments");

if(clRangeListFile || clFilterFile){

    if(clRangeListFile)
        readChromRangeFile(clRangeListFile);
    if(clFilterFile)
        readChromRangeFile(clFilterFile);
    //printRangeList();
    
}

if(gBestNum>0){
    if(gNumSortCriteria==0){
        printf("Error: !BEST specified but no sort criteria !MIN or !MAX specified\n");
        clean_up();
        usage();
        return 1;
    }
}

if(gNumSortCriteria>0){
    if(gBestNum==0){
        printf("Error: sort criteria !MIN or !MAX specified but !BEST not specified\n");
        clean_up();
        usage();
        return 1;
    }
}



//transferSortCriteriaToArray();
//printSortCriteria();

//printf("NUMBEST=%d\n",gBestNum);

if(argc>2){
    for(int i=2;i<argc;i++){
        replaceCharInString(argv[i],'/','\t');
        parseFilterLine(argv[i]);
    }
}

if(gInBB){
    if(gPrintComment) printf("#inputbb=%s\n",gInBB);
    bigBedToBed(gInBB);
    
}else{
    if (argc < 2)
        usage();
    if(gPrintComment) printf("#inputbb=%s\n",argv[1]);
    bigBedToBed(argv[1]);
}

if (verboseLevel() > 1)
    printVmPeak();

clean_up();

return 0;
}
