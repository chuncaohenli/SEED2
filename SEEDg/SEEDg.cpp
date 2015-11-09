//**********************************************************************************
//* Title: SEEDg
//* Platform: 32-Bit/64-Bit Windows/Linux/Mac
//* Author: Gao Xiaoxu
//* Affliation: Bao Lab, BJTU
//* Date: 2015.11
//* Copy Right: For Purpose of Study Only
//**********************************************************************************


#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <cstdlib>
#include <cstring>
#include <stdlib.h>

#include <string>
#include <algorithm>
#include <map>
#include "k-means.h"

#define LENGTH_DIFFERENCE 10
#define BOWTIE_DIFFERENCE 5

using namespace std;

int QV = 0;
int reversed = 0;
int paired = 0;

int seedsWeight = 16 * 1024;

//GXX
bool is_firstSetSeed = true;
//char core_seed[100][100];
ofstream ofile;
ofstream ofile_test("output_sid_seq.txt");


#define OFFSET 33
#define RANGE 94 
//GXX
int **seeds;
int seedsCount = 10;
int seeds_length = 30;
int seeds_weight = 12;

static int fastSeeds[4][52] =
{
     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
};

static int shortSeeds[10][15] =
{
     1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0,
     1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0,
     1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,
     0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0,
     0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1,
     0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1  
};

class Hash
{
     ifstream in;
     int lowerSizeInBit, lowerSizeInChar, upperSizeInBit, upperSizeInChar, num, mismatchAllowed;
     char * seqHead;
     unsigned int ** indexHead;
//   unsigned int * offsetCount;
public:
     unsigned int * offsetCount;
     Hash(char [], int, int, int);
     void build();
     void seqInsert(char [], int, char *, unsigned int);
     void indexInsert(char [], unsigned int);
     char change(char);
     char changeBack(char);
     unsigned int calOffset(int, char []);
     int searchByIndex(int, unsigned int, int, char []);
     int searchBySeq(unsigned int, char []);
     void deleteByIndex(int, unsigned int, int);
     void deleteBySeq(unsigned int);
     int calSeqID(int, unsigned int, int);
     void QVInsert(char [], int, char *, unsigned int);
     int searchByIndex(int, unsigned int, int, char [], char [], unsigned int &);
     int searchBySeq(unsigned int, char [], char [], unsigned int &);
     void seqInsert(char [], int, char *, unsigned int, unsigned int);
     int searchByIndex(int, unsigned int, int, char [], unsigned int &);
     int searchBySeq(unsigned int, char [], unsigned int &);
     void adjust();
     void tmpDeleteByIndex(int, unsigned int, int);
     void recoverByIndex(int, unsigned int, int);
     unsigned int getOffsetCount(unsigned int, int);

     ~Hash();
};

class FastqGenerator
{
     ifstream in;
     char addiInput[100];
     ifstream addiIn;
     char outputq[100];
     ofstream out;
     int num;
     char * seq;
public:
     FastqGenerator(char [], char [], int);
     FastqGenerator(char [], char [], int, int);
     void record();
     void generateFastq();
};

class Cluster
{
     ofstream out;
//   ofstream dis;
     int lowerSizeInChar;
     int upperSizeInChar;
     int lowerSizeInBit;
     int upperSizeInBit;
     int num;
     int mismatchAllowed;
     int shiftAllowed;
     int CLID;
//   int numInCL;
//   Hash * h;
     int lowerQV;
     int upperQV;
     ofstream addiOut;
     char addiOutput[100];
     unsigned long seqNum;
     unsigned long adjustNum;
     unsigned int ** mappingTable;
     unsigned int * mappingNum;
     char midInput[100];
public:
     Hash * h;
     Cluster(char [], char [], int, int, int, int, int, unsigned int **, unsigned int *);
     void cluster();
     int compare(char [], char [], int, int, int &);
     void clusterWithMismatches(char []);
     void clusterWithShifts(char []);
     char max(int, int, int, int);
     void clusterByConsensus();
     void calConsensus(char [], unsigned int, int &);
     void calConsensus(char [], char [], unsigned int, int &);
     void preprocess();
     Cluster(char [], char [], int, int, int, int, int, int, int, unsigned int **, unsigned int *);
     int compare(char [], char [], char [], char [], int, int, int &);
     void clusterWithMismatches(char [], char []);
     void clusterWithShifts(char [], char []);
     char reverseChange(char);

};

class FileAnalyzer
{
public:
     void inputAnalyze(char [], int &, int &, int &, int &);
     void outputAnalyze(int, int);
     void PECombine(char [], int, char [], int, char *, int &, int &, int &);
};

class Sorter
{
     ifstream in;
     ofstream midOut;
     char midOutput[100];
     int num;
     unsigned int ** mappingTable;
     unsigned int * mappingNum;
     int lowerSizeInChar;
     int realNum;
     typedef struct 
     {
          int ID;
          int realID;
     } Order;
public:
     Sorter(char [], int, int);
     void sort();
     void suffixSort(int, int, int, char [], Order []);
     int getRealNum();
     unsigned int ** getMappingTable();
     unsigned int * getMappingNum();
};

Hash::Hash(char input[], int num, int lowerSizeInChar, int upperSizeInChar)
{
     long int i;

     if(upperSizeInChar % 4)
          this->upperSizeInBit = upperSizeInChar / 4 + 1;
     else
          this->upperSizeInBit = upperSizeInChar / 4;
     if(lowerSizeInChar % 4)
          this->lowerSizeInBit = lowerSizeInChar / 4 + 1;
     else
          this->lowerSizeInBit = lowerSizeInChar / 4;
     this->lowerSizeInChar = lowerSizeInChar;
     this->upperSizeInChar = upperSizeInChar;
     this->num = num;
     if(QV)
     {
          seqHead = new char[(long int)num * (upperSizeInBit + 5 + upperSizeInChar)];
          for(i = 0; i < (long int)num * (upperSizeInBit + 5 + upperSizeInChar); i ++)
               seqHead[i] = 0;
     }
     else
     {
          seqHead = new char[(long int)num * (upperSizeInBit + 5)];
          for(i = 0; i < (long int)num * (upperSizeInBit + 5); i ++)
               seqHead[i] = 0;
     }
//   indexHead = new unsigned int * [1024 * 1024 * 16 * 10];
//   indexHead = new unsigned int * [1024 * 1024 * 64 * 4];
     indexHead = new unsigned int * [1024 * seedsWeight * seedsCount];
     offsetCount = new unsigned int [1024 * seedsWeight * seedsCount];
     for(i = 0; i < (long int)1024 * seedsWeight * seedsCount; i ++)
          offsetCount[i] = 0;
     in.open(input);
}

Hash::~Hash(){
     delete[] indexHead;
     delete[] offsetCount;
     delete[] seqHead;
}

void Hash::build()
{
     char buf[1001];
     int i, count = 0, tag = 1;
     unsigned int seqOffset = 0, seqID;

//   int seqNext;
//   char base[4];
//   ofstream out;

     int totalLength, j;

     if(in.is_open())
     {
          seqID = 0;
cont:
          while(seqID < num * 4)
          {
               in.getline(buf, 1001);
               if(seqID % 4 == 1)
               {
                    for(i = 0; i < in.gcount() - 1; i ++)
                    {
                         if(buf[i] == 'N')
                         {
                              tag = 0;
                              seqID ++;
                              goto cont;
                         }//to deal with N
                         else
                         {
                              tag = 1;
                              buf[i] = change(buf[i]);
                         }
                    }
                    seqInsert(buf, in.gcount() - 1, seqHead, seqOffset, seqID/4);
                    indexInsert(buf, seqOffset);
                    seqOffset = seqOffset + (upperSizeInBit + 5);
               }
               if(QV)
               {
                    if(seqID % 4 == 3 && tag == 1)
                    {
                         if(buf[in.gcount()] == '\n')
                              QVInsert(buf, in.gcount() - 1, seqHead, seqOffset);
                         else
                              QVInsert(buf, in.gcount(), seqHead, seqOffset);//for the last line in Windows OS
                         seqOffset = seqOffset + upperSizeInChar;
                    }
                    for(i = 0; i < 1001; i ++)
                         buf[i] = 0;
               }
               seqID ++;
          }
     }
     else
     {
          cout << "CANNOT OPEN INPUT FILE!" << endl;
          exit(-1);
     }

//verification
//   cout << endl;
//#ifdef QV
//#ifdef REALID
//   for(seqOffset = 0; seqOffset < num * (upperSizeInBit + 5 + upperSizeInChar); seqOffset = seqOffset + upperSizeInBit + 5 + upperSizeInChar)
//#else
//   for(seqOffset = 0; seqOffset < num * (upperSizeInBit + 1 + upperSizeInChar); seqOffset = seqOffset + upperSizeInBit + 1 + upperSizeInChar)
//#endif
//#else
//#ifdef REALID
//   for(seqOffset = 0; seqOffset < num * (upperSizeInBit + 5); seqOffset = seqOffset + upperSizeInBit + 5)
//#else
//   for(seqOffset = 0; seqOffset < num * (upperSizeInBit + 1); seqOffset = seqOffset + upperSizeInBit + 1)
//#endif
//#endif
//   {
//        cout << (int)*(seqHead + seqOffset) << ": " << endl;
//        for(seqNext = 1; seqNext < 5; seqNext ++)
//        {
//             cout << (int)*(seqHead + seqOffset + seqNext) << " ";
//        }
//        cout << endl;
//#ifdef REALID
//        for(seqNext = 5; seqNext < upperSizeInBit + 5; seqNext ++)
//#else
//        for(seqNext = 1; seqNext < upperSizeInBit + 1; seqNext ++)
//#endif
//        {
//             base[0] = (*(seqHead + seqOffset + seqNext) >> 6) & 0x03;
//             base[1] = (*(seqHead + seqOffset + seqNext) >> 4) & 0x03;
//             base[2] = (*(seqHead + seqOffset + seqNext) >> 2) & 0x03;
//             base[3] = (*(seqHead + seqOffset + seqNext)) & 0x03;
//             cout << (unsigned int)base[0] << " " << (unsigned int)base[1] << " " << (unsigned int)base[2] << " " << (unsigned int)base[3] << " ";
//        }
//#ifdef QV
//        cout << endl;
//#ifdef REALID
//        for(; seqNext < upperSizeInBit + 1 + upperSizeInChar; seqNext ++)
//#else
//        for(; seqNext < upperSizeInBit + 5 + upperSizeInChar; seqNext ++)
//#endif
//             cout << *(seqHead + seqOffset + seqNext) << " ";
//#endif
//        cout << endl;
//   }
//verification
/*
        totalLength = 0;
        count = 0;
        for(i = 0; i < seedsCount; i ++)
        {
                for(j = 0; j < 1024 * seedsWeight; j ++)
                        if(c.h->offsetCount[j * seedsCount + i] > 1000)
                        {
                                totalLength = totalLength + c.h->offsetCount[j * seedsCount + i];
                                count ++;
                        }
        }
        cout << "#buckets longer than 1000: " << count << endl;
     if(count != 0)
             cout << "average length of the buckets: " << totalLength / count << endl;
*/
}

void Hash::adjust()
{
     int i, j, k ,p;

     for(i = 0; i < seedsCount; i ++)
          for(j = 0; j < 1024 * seedsWeight; j ++)
               if(offsetCount[j * seedsCount + i] != 0)
               {
                    p = offsetCount[j * seedsCount + i] - 1;
                    for(k = 0; k < offsetCount[j * seedsCount + i] && p > k; k ++)
                    {
                         if(seqHead[indexHead[j * seedsCount + i][k]] == 0)
                         {
                              while(seqHead[indexHead[j * seedsCount + i][p]] == 0 && p > k)
                              {
                                   p --;
                                   offsetCount[j * seedsCount + i] --;
                              }
                              if(p > k)
                              {
                                   indexHead[j * seedsCount + i][k] = indexHead[j * seedsCount + i][p];
                                   p --;
                                   offsetCount[j * seedsCount + i] --;
                              }
                         }
                    }
                    if(seqHead[indexHead[j * seedsCount + i][k]] == 0)
                         offsetCount[j * seedsCount + i] --;
               }
}

void Hash::seqInsert(char buf[], int realSize, char * seqHead, unsigned int seqOffset, unsigned int seqID)
{
     int i, seqNext;
     char bitBuf = 0x00;

     *(seqHead + seqOffset) = 0x01;// the first byte was intended to record size of the seq, but is now used for existence of the seq. If it is changed for the recording purpose in the future, more bytes would be required.
     *(seqHead + seqOffset + 1) = (char) ((seqID & 0xff000000) >> 24);
     *(seqHead + seqOffset + 2) = (char) ((seqID & 0x00ff0000) >> 16);
     *(seqHead + seqOffset + 3) = (char) ((seqID & 0x0000ff00) >> 8);
     *(seqHead + seqOffset + 4) = (char) (seqID & 0x000000ff);
     for(i = 0, seqNext = 5; i < realSize; i ++)
     {
          if((i + 1) % 4 == 0)
          {
               bitBuf = (bitBuf | buf[i]);
               *(seqHead + seqOffset + seqNext) = bitBuf;
               seqNext ++;
               bitBuf = 0x00;
          }
          else if(i == realSize - 1)
          {
               bitBuf = (bitBuf | buf[i]) << (4 - realSize % 4) * 2;
               *(seqHead + seqOffset + seqNext) = bitBuf;
               seqNext ++;
               bitBuf = 0x00;
          }
          else
               bitBuf = (bitBuf | buf[i]) << 2;
     }
}

void Hash::QVInsert(char buf[], int realSize, char * seqHead, unsigned int seqOffset)
{
     int i, seqNext = 0;

     for(i = 0; i < realSize; i ++)
     {
          *(seqHead + seqOffset + seqNext) = buf[i];
          seqNext ++;
     }
}

void Hash::indexInsert(char buf[], unsigned int seqOffset)
{
     int indexNext;
     unsigned int indexOffset = 0;

     for(indexNext = 0; indexNext < seedsCount; indexNext ++)
     {
          indexOffset = calOffset(indexNext, buf);
          if(offsetCount[indexOffset + indexNext] == 0)
          {
               indexHead[indexOffset + indexNext] = (unsigned int *) malloc((++ offsetCount[indexOffset + indexNext]) * sizeof(unsigned int));
               if(indexHead[indexOffset + indexNext] == NULL)
               {
                    cout << "CANNOT ALLOCATE MEMORY!" << endl;
                    exit(-1);
               }
          }
          else
          {
               indexHead[indexOffset + indexNext] = (unsigned int *) realloc(indexHead[indexOffset + indexNext], (++ offsetCount[indexOffset + indexNext]) * sizeof(unsigned int));
               if(indexHead[indexOffset + indexNext] == NULL)
               {
                    cout << "CANNOT ALLOCATE MEMORY!" << endl;
                    exit(-1);
               }
          }
          indexHead[indexOffset + indexNext][offsetCount[indexOffset + indexNext] - 1] = seqOffset;
     }
}

char Hash::change(char base)
{
     switch(base)
     {
          case 'A': return 0x00;
          case 'C': return 0x01;
          case 'G': return 0x02;
          case 'T': return 0x03;
          default: cout << "INPUT ERROR!" << endl; exit(-1);
     }
}

char Hash::changeBack(char base)
{
     switch(base)
     {
          case 0x00: return 'A';
          case 0x01: return 'C';
          case 0x02: return 'G';
          case 0x03: return 'T';
          default: cout << "MEMORY ERROR!" << endl; exit(-1);
     }
}

unsigned int Hash::calOffset(int indexNext, char buf[])
{
     unsigned int indexOffset = 0;
     int i, j = 0;

     if(seedsWeight == 1024 * 16)
     {
          for(i = 0; i < 30; i ++)
               if(seeds[indexNext][i] == 1)
                    indexOffset = indexOffset + buf[i + lowerSizeInChar - 33] * (unsigned int)pow(4, j ++);
     }
     else if(seedsWeight == 1024 * 64)
     {
          for(i = 0; i < 52; i ++)
               if(fastSeeds[indexNext][i] == 1)
                    indexOffset = indexOffset + buf[i + lowerSizeInChar - 55] * (unsigned int)pow(4, j ++);
     }
     else
     {
          for(i = 0; i < 15; i ++)
                        if(shortSeeds[indexNext][i] == 1)
                                indexOffset = indexOffset + buf[i + lowerSizeInChar - 18] * (unsigned int)pow(4, j ++);
     }
     return indexOffset * seedsCount;
}

int Hash::searchByIndex(int indexNext, unsigned int indexOffset, int no, char buf[], unsigned int & seqID)
{
     int seqNext, p = 0;
     int i;

     seqID = 0;
     for(i = 1; i < 5; i ++)
     {
          seqID = seqID | (((unsigned int) seqHead[indexHead[indexOffset + indexNext][no] + i]) & 0x000000ff);
          if(i < 4)
               seqID = seqID << 8;
     }
     for(seqNext = 5; seqNext < upperSizeInBit + 5; seqNext ++)
     {
          buf[p ++] = (seqHead[indexHead[indexOffset + indexNext][no] + seqNext] >> 6) & 0x03;
          buf[p ++] = (seqHead[indexHead[indexOffset + indexNext][no] + seqNext] >> 4) & 0x03;
          buf[p ++] = (seqHead[indexHead[indexOffset + indexNext][no] + seqNext] >> 2) & 0x03;
          buf[p ++] = (seqHead[indexHead[indexOffset + indexNext][no] + seqNext]) & 0x03;
     }
     return (int)seqHead[indexHead[indexOffset + indexNext][no]];
}

int Hash::searchBySeq(unsigned int seqOffset, char buf[], unsigned int & seqID)
{
     int seqNext, p = 0;
     int i;

     seqID = 0;
     for(i = 1; i < 5; i ++)
     {
          seqID = seqID | (((unsigned int) seqHead[seqOffset + i]) & 0x000000ff);
          if(i < 4)
               seqID = seqID << 8;
     }
     for(seqNext = 5; seqNext < upperSizeInBit + 5; seqNext ++)
     {
          buf[p ++] = (seqHead[seqOffset + seqNext] >> 6) & 0x03;
          buf[p ++] = (seqHead[seqOffset + seqNext] >> 4) & 0x03;
          buf[p ++] = (seqHead[seqOffset + seqNext] >> 2) & 0x03;
          buf[p ++] = (seqHead[seqOffset + seqNext]) & 0x03;
     }
     return (int)seqHead[seqOffset];
}

int Hash::searchByIndex(int indexNext, unsigned int indexOffset, int no, char buf[], char QVBuf[], unsigned int & seqID)
{
     int seqNext, p = 0;
     int i;

     seqID = 0;
     for(i = 1; i < 5; i ++)
     {
          seqID = seqID | (((unsigned int) seqHead[indexHead[indexOffset + indexNext][no] + i]) & 0x000000ff);
          if(i < 4)
               seqID = seqID << 8;
     }
     for(seqNext = 5; seqNext < upperSizeInBit + 5; seqNext ++)
     {
          buf[p ++] = (seqHead[indexHead[indexOffset + indexNext][no] + seqNext] >> 6) & 0x03;
          buf[p ++] = (seqHead[indexHead[indexOffset + indexNext][no] + seqNext] >> 4) & 0x03;
          buf[p ++] = (seqHead[indexHead[indexOffset + indexNext][no] + seqNext] >> 2) & 0x03;
          buf[p ++] = (seqHead[indexHead[indexOffset + indexNext][no] + seqNext]) & 0x03;
     }
     p = 0;
     for(; seqNext < upperSizeInBit + 5 + upperSizeInChar; seqNext ++)
          QVBuf[p ++] = seqHead[indexHead[indexOffset + indexNext][no] + seqNext];
     return (int)seqHead[indexHead[indexOffset + indexNext][no]];
}

int Hash::searchBySeq(unsigned int seqOffset, char buf[], char QVBuf[], unsigned int & seqID)
{
     int seqNext, p = 0;
     int i;

     seqID = 0;
     for(i = 1; i < 5; i ++)
     {
          seqID = seqID | (((unsigned int) *(seqHead + seqOffset + i)) & 0x000000ff);
          if(i < 4)
               seqID = seqID << 8;
     }
     for(seqNext = 5; seqNext < upperSizeInBit + 5; seqNext ++)
     {
          buf[p ++] = (seqHead[seqOffset + seqNext] >> 6) & 0x03;
          buf[p ++] = (seqHead[seqOffset + seqNext] >> 4) & 0x03;
          buf[p ++] = (seqHead[seqOffset + seqNext] >> 2) & 0x03;
          buf[p ++] = (seqHead[seqOffset + seqNext]) & 0x03;
     }
     p = 0;
     for(; seqNext < upperSizeInBit + 5 + upperSizeInChar; seqNext ++)
          QVBuf[p ++] = seqHead[seqOffset + seqNext];
     return (int)seqHead[seqOffset];
}

void Hash::deleteByIndex(int indexNext, unsigned int indexOffset, int no)
{
     seqHead[indexHead[indexOffset + indexNext][no]] = 0;
}

void Hash::tmpDeleteByIndex(int indexNext, unsigned int indexOffset, int no)
{
     seqHead[indexHead[indexOffset + indexNext][no]] = seqHead[indexHead[indexOffset + indexNext][no]] | 0x80;
}

void Hash::recoverByIndex(int indexNext, unsigned int indexOffset, int no)
{
     seqHead[indexHead[indexOffset + indexNext][no]] = seqHead[indexHead[indexOffset + indexNext][no]] & 0x7f;
}

void Hash::deleteBySeq(unsigned int seqOffset)
{
     seqHead[seqOffset] = 0;
}

int Hash::calSeqID(int indexNext, unsigned int indexOffset, int no)
{
     if(QV)
          return indexHead[indexOffset + indexNext][no]/(upperSizeInBit + 1 + upperSizeInChar);
     else
          return indexHead[indexOffset + indexNext][no]/(upperSizeInBit + 1);
}

unsigned int Hash::getOffsetCount(unsigned int indexOffset, int indexNext)
{
     return offsetCount[indexOffset + indexNext];
}

Cluster::Cluster(char input[], char output[], int num, int lowerSizeInChar, int upperSizeInChar, int mismatchAllowed, int shiftAllowed, int lowerQV, int upperQV, unsigned int ** mappingTable, unsigned int * mappingNum)
{
     if(lowerSizeInChar % 4)
          this->lowerSizeInBit = lowerSizeInChar / 4 + 1;
     else
          this->lowerSizeInBit = lowerSizeInChar / 4;
     if(upperSizeInChar % 4)
          this->upperSizeInBit = upperSizeInChar / 4 + 1;
     else
          this->upperSizeInBit = upperSizeInChar / 4;
     this->lowerSizeInChar = lowerSizeInChar;
     this->upperSizeInChar = upperSizeInChar;
     this->num = num;
     this->mismatchAllowed = mismatchAllowed;
     this->shiftAllowed = shiftAllowed;
     this->CLID = 0;
     this->lowerQV = lowerQV;
     this->upperQV = upperQV;
     this->seqNum = 0;
     this->adjustNum = num / 20;
     this->mappingTable = mappingTable;
     this->mappingNum = mappingNum;
     strcpy(midInput, input);
     strcat(midInput, ".mid.fastq");
     out.open(output);
     strcpy(addiOutput, output);
     strcat(addiOutput, ".fasta");
     addiOut.open(addiOutput);
//   dis.open("distribution.txt");
     out << "CLID   SeqID" << endl;
//   dis << "CLID   No" << endl;
     h = new Hash(midInput, num, lowerSizeInChar, upperSizeInChar);
     h->build();
}

Cluster::Cluster(char input[], char output[], int num, int lowerSizeInChar, int upperSizeInChar, int mismatchAllowed, int shiftAllowed, unsigned int ** mappingTable, unsigned int * mappingNum)
{
     if(lowerSizeInChar % 4)
          this->lowerSizeInBit = lowerSizeInChar / 4 + 1;
     else
          this->lowerSizeInBit = lowerSizeInChar / 4;
     if(upperSizeInChar % 4)
          this->upperSizeInBit = upperSizeInChar / 4 + 1;
     else
          this->upperSizeInBit = upperSizeInChar / 4;
     this->lowerSizeInChar = lowerSizeInChar;
     this->upperSizeInChar = upperSizeInChar;
     this->num = num;
     this->mismatchAllowed = mismatchAllowed;
     this->shiftAllowed = shiftAllowed;
     this->CLID = 0;
     this->seqNum = 0;
     this->adjustNum = num / 20;
     this->mappingTable = mappingTable;
     this->mappingNum = mappingNum;
     strcpy(midInput, input);
     strcat(midInput, ".mid.fastq");
     out.open(output);
     strcpy(addiOutput, output);
     strcat(addiOutput, ".fasta");
     addiOut.open(addiOutput);
//   dis.open("distribution.txt");
     out << "CLID   SeqID" << endl;
//   dis << "CLID   No" << endl;
     h = new Hash(midInput, num, lowerSizeInChar, upperSizeInChar);
     h->build();
}

char Cluster::max(int a, int b, int c, int d)
{
     if(b >= a && b >= c && b >= d) return 0x01;
     if(c >= a && c >= b && c >= d) return 0x02;
     if(d >= a && d >= b && d >= c) return 0x03;
     if(a >= b && a >= c && a >= d) return 0x00;
     cout << "UNKNOWN ERROR!" << endl; exit(-1);
}

void Cluster::calConsensus(char sBuf[], unsigned int sSeqID, int & tagReverse)
{
     long int A[1000] = {0}, C[1000] = {0}, G[1000] = {0}, T[1000] = {0};
     int indexNext, no, i, realSize, similarity = 1000, diff, j;
     char sBufBak[1000], tBuf[1000], buf[1000];
     unsigned int indexOffset;
     unsigned int seqID, centerSeqID;

     for(indexNext = 0; indexNext < seedsCount; indexNext ++)
     {
          indexOffset = h->calOffset(indexNext, sBuf);
          for(no = 0; no < h->getOffsetCount(indexOffset, indexNext); no ++)
          {
               realSize = h->searchByIndex(indexNext, indexOffset, no, tBuf, seqID);
               if(realSize > 0 && compare(sBuf, tBuf, 0, lowerSizeInChar, tagReverse) <= mismatchAllowed * 2)
               {
                    if(reversed && tagReverse)
                         for(i = 0, j = lowerSizeInChar - 1; i < lowerSizeInChar; i ++, j --)
                              switch(tBuf[j])
                              {
                                   case 0x00: T[i] = T[i] + mappingNum[seqID]; break;
                                   case 0x01: G[i] = G[i] + mappingNum[seqID]; break;
                                   case 0x02: C[i] = C[i] + mappingNum[seqID]; break;
                                   case 0x03: A[i] = A[i] + mappingNum[seqID]; break;
                                   default: cout << "MEMORY ERROR!" << endl; exit(-1);
                              }    
                    else
                         for(i = 0; i < lowerSizeInChar; i ++)
                              switch(tBuf[i])
                              {
                                   case 0x00: A[i] = A[i] + mappingNum[seqID]; break;
                                   case 0x01: C[i] = C[i] + mappingNum[seqID]; break;
                                   case 0x02: G[i] = G[i] + mappingNum[seqID]; break;
                                   case 0x03: T[i] = T[i] + mappingNum[seqID]; break;
                                   default: cout << "MEMORY ERROR!" << endl; exit(-1);
                              }
                    h->tmpDeleteByIndex(indexNext, indexOffset, no);
               }
          }
     }
     for(i = 0; i < lowerSizeInChar; i ++)
     {
          sBufBak[i] = sBuf[i];
          sBuf[i] = max(A[i], C[i], G[i], T[i]);
     }

     //GXX 5-6  
     if (ofile.is_open() && is_firstSetSeed){
          
          ofile << "@1" << endl;
          for(i = 0; i < lowerSizeInChar; i ++)
          {

               // 转码 GXX   
               if ((int)sBuf[i] == 0)
               {
                    //cout << 'A';
                    ofile << 'A';

               } 
               else if((int)sBuf[i] == 1)
               {
                    ofile << 'C';
               } else if ((int)sBuf[i] == 2){
                    ofile << 'G';
               } else{
                    ofile << 'T';
               }
               //cout <<  (int)sBuf[i];
          }
          ofile << endl;
          ofile << "+" << endl <<endl;
     } else {
          //cout << "Unable to open file"; 
          return ;
     }    
//find the most similar seq to the consensus and write to the output
     for(indexNext = 0; indexNext < seedsCount; indexNext ++)
     {
          indexOffset = h->calOffset(indexNext, sBufBak);
          for(no = 0; no < h->getOffsetCount(indexOffset, indexNext); no ++)
          {
               realSize = h->searchByIndex(indexNext, indexOffset, no, tBuf, seqID);
               if(realSize < 0)
               {
                    h->recoverByIndex(indexNext, indexOffset, no);
                    diff = compare(sBuf, tBuf, 0, lowerSizeInChar, tagReverse);
                    if(diff < similarity)
                    {
                         similarity = diff;
                         for(i = 0; i < lowerSizeInChar; i ++)
                              buf[i] = tBuf[i];
                         centerSeqID = seqID;
                    }
               }
          }
     }
//if the virtual center and the center are too far
     if(compare(sBuf, buf, 0, lowerSizeInChar, tagReverse) <= mismatchAllowed)// bug exists here: the same center sequence is found for many clusters since the center sequence is not deleted
     {
          if(paired == 0)
          {
               addiOut << ">" << mappingTable[centerSeqID][0] << endl;
               for(i = 0; i < lowerSizeInChar; i ++)
               {
                    out << h->changeBack(buf[i]);
                    addiOut << h->changeBack(buf[i]);
               }
               out << endl;
               addiOut << endl;
          }
          else
          {
               addiOut << ">" << mappingTable[centerSeqID][0] << endl;
               for(i = 0; i < paired; i ++)//paired == lower
               {
                    out << h->changeBack(buf[i]);
                    addiOut << h->changeBack(buf[i]);
               }
               out << endl;
               addiOut << endl;
               addiOut << ">" << mappingTable[centerSeqID][0] << endl;
               for(; i < lowerSizeInChar; i ++)
               {
                    out << h->changeBack(buf[i]);
                    addiOut << h->changeBack(buf[i]);
               }
               out << endl;
               addiOut << endl;
          }
     }
     else
     {
          if(paired == 0)
          {
               addiOut << ">" << mappingTable[sSeqID][0] << endl;
               for(i = 0; i < lowerSizeInChar; i ++)
               {
                    out << h->changeBack(sBufBak[i]); //out << h->changeBack(sBuf[i]);
                    addiOut << h->changeBack(sBufBak[i]); //addiOut << h->changeBack(sBuf[i]);
               }
               out << endl;
               addiOut << endl;
          }
          else
          {
               addiOut << ">" << mappingTable[sSeqID][0] << endl;
               for(i = 0; i < paired; i ++)//paired == lower
               {
                    out << h->changeBack(sBufBak[i]);
                    addiOut << h->changeBack(sBufBak[i]);
               }
               out << endl;
               addiOut << endl;
               addiOut << ">" << mappingTable[sSeqID][0] << endl;
               for(; i < lowerSizeInChar; i ++)
               {
                    out << h->changeBack(sBufBak[i]);
                    addiOut << h->changeBack(sBufBak[i]);
               }
               out << endl;
               addiOut << endl;
          }
//still need to decide if the source sequence is reverse complementary to the virtual center
          compare(sBuf, sBufBak, 0, lowerSizeInChar, tagReverse);
     }
}

void Cluster::calConsensus(char sBuf[], char sQVBuf[], unsigned int sSeqID, int & tagReverse)
{
        long int A[1000] = {0}, C[1000] = {0}, G[1000] = {0}, T[1000] = {0}, QV[1000] = {0};
        int indexNext, no, i, realSize, similarity = 1000, diff, j;
        char sBufBak[1000], tBuf[1000], tQVBuf[1000], buf[1000], QVBuf[1000];
        unsigned int indexOffset;
        unsigned int seqID, centerSeqID;

        for(indexNext = 0; indexNext < seedsCount; indexNext ++)
        {
                indexOffset = h->calOffset(indexNext, sBuf);
                for(no = 0; no < h->getOffsetCount(indexOffset, indexNext); no ++)
                {
                        realSize = h->searchByIndex(indexNext, indexOffset, no, tBuf, tQVBuf, seqID);
                        if(realSize > 0 && compare(sBuf, sQVBuf, tBuf, tQVBuf, 0, lowerSizeInChar, tagReverse) <= mismatchAllowed * 2)
                        {
                                if(reversed && tagReverse)
                                        for(i = 0, j = lowerSizeInChar - 1; i < lowerSizeInChar; i ++, j --)
                         {
                                                switch(tBuf[j])
                                                {
                                                        case 0x00: T[i] = T[i] + mappingNum[seqID]; break;
                                                        case 0x01: G[i] = G[i] + mappingNum[seqID]; break;
                                                        case 0x02: C[i] = C[i] + mappingNum[seqID]; break;
                                                        case 0x03: A[i] = A[i] + mappingNum[seqID]; break;
                                                        default: cout << "MEMORY ERROR!" << endl; exit(-1);
                                                }
                              QV[i] = QV[i] + mappingNum[seqID] * tQVBuf[j];
                         }
                                else
                                        for(i = 0; i < lowerSizeInChar; i ++)
                         {
                                                switch(tBuf[i])
                                                {
                                                        case 0x00: A[i] = A[i] + mappingNum[seqID]; break;
                                                        case 0x01: C[i] = C[i] + mappingNum[seqID]; break;
                                                        case 0x02: G[i] = G[i] + mappingNum[seqID]; break;
                                                        case 0x03: T[i] = T[i] + mappingNum[seqID]; break;
                                                        default: cout << "MEMORY ERROR!" << endl; exit(-1);
                                                }
                              QV[i] = QV[i] + mappingNum[seqID] * tQVBuf[i];
                         }
                                h->tmpDeleteByIndex(indexNext, indexOffset, no);
                        }
                }
        }
        for(i = 0; i < lowerSizeInChar; i ++)
        {
                sBufBak[i] = sBuf[i];
                sBuf[i] = max(A[i], C[i], G[i], T[i]);
          sQVBuf[i] = QV[i] / (A[i] + C[i] + G[i] + T[i]);
        }
//find the most similar seq to the consensus and write to the output
        for(indexNext = 0; indexNext < seedsCount; indexNext ++)
        {
                indexOffset = h->calOffset(indexNext, sBufBak);
                for(no = 0; no < h->getOffsetCount(indexOffset, indexNext); no ++)
                {
                        realSize = h->searchByIndex(indexNext, indexOffset, no, tBuf, tQVBuf, seqID);
                        if(realSize < 0)
                        {
                                h->recoverByIndex(indexNext, indexOffset, no);
                                diff = compare(sBuf, sQVBuf, tBuf, tQVBuf, 0, lowerSizeInChar, tagReverse);
                                if(diff < similarity)
                                {
                                        similarity = diff;
                                        for(i = 0; i < lowerSizeInChar; i ++)
                         {
                                                buf[i] = tBuf[i];
                              QVBuf[i] = tQVBuf[i];
                         }
                                        centerSeqID = seqID;
                                }
                        }
                }
        }
//if the virtual center and the center are too far
        if(compare(sBuf, sQVBuf, buf, QVBuf, 0, lowerSizeInChar, tagReverse) <= mismatchAllowed)
        {
                if(paired == 0)
                {
                        addiOut << ">" << mappingTable[centerSeqID][0] << endl;
                        for(i = 0; i < lowerSizeInChar; i ++)
                        {
                                out << h->changeBack(buf[i]);
                                addiOut << h->changeBack(buf[i]);
                        }
                        out << endl;
                        addiOut << endl;
                }
                else
                {
                        addiOut << ">" << mappingTable[centerSeqID][0] << endl;
                        for(i = 0; i < paired; i ++)//paired == lower
                        {
                                out << h->changeBack(buf[i]);
                                addiOut << h->changeBack(buf[i]);
                        }
                        out << endl;
                        addiOut << endl;
                        addiOut << ">" << mappingTable[centerSeqID][0] << endl;
                        for(; i < lowerSizeInChar; i ++)
                        {
                                out << h->changeBack(buf[i]);
                                addiOut << h->changeBack(buf[i]);
                        }
                        out << endl;
                        addiOut << endl;
                }
        }
        else
        {
                if(paired == 0)
                {
                        addiOut << ">" << mappingTable[sSeqID][0] << endl;
                        for(i = 0; i < lowerSizeInChar; i ++)
                        {
                                out << h->changeBack(sBufBak[i]); //out << h->changeBack(sBuf[i]);
                                addiOut << h->changeBack(sBufBak[i]); //addiOut << h->changeBack(sBuf[i]);
                        }
                        out << endl;
                        addiOut << endl;
                }
                else
                {
                        addiOut << ">" << mappingTable[sSeqID][0] << endl;
                        for(i = 0; i < paired; i ++)//paired == lower
                        {
                                out << h->changeBack(sBufBak[i]);
                                addiOut << h->changeBack(sBufBak[i]);
                        }
                        out << endl;
                        addiOut << endl;
                        addiOut << ">" << mappingTable[sSeqID][0] << endl;
                        for(; i < lowerSizeInChar; i ++)
                        {
                                out << h->changeBack(sBufBak[i]);
                                addiOut << h->changeBack(sBufBak[i]);
                        }
                        out << endl;
                        addiOut << endl;
                }
//still need to decide if the source sequence is reverse complementary to the virtual center
                compare(sBuf, sBufBak, 0, lowerSizeInChar, tagReverse);
        }
}

void Cluster::cluster()
{
     char sBuf[1000];
     int realSize, tagReverse, i;
     unsigned int seqOffset;
     char sQVBuf[1000];
     unsigned int seqID;
//   long int big = 0, small = 0;

     if(QV)
     {
          for(seqOffset = 0; seqOffset < num * (upperSizeInBit + 5 + upperSizeInChar); seqOffset = seqOffset + (upperSizeInBit + 5 + upperSizeInChar))
          {
               realSize = h->searchBySeq(seqOffset, sBuf, sQVBuf, seqID);
               if(realSize == 0) continue;
               else realSize = 0;
               calConsensus(sBuf, sQVBuf, seqID, tagReverse);
               if(out.is_open())
                    if(reversed)
                         for(i = 0; i < mappingNum[seqID]; i ++)
                              out << CLID << "    " << mappingTable[seqID][i] << "   " << tagReverse << endl;
                    else
                         for(i = 0; i < mappingNum[seqID]; i ++)
                              out << CLID << "    " << mappingTable[seqID][i] << endl;
               else
               {
                    cout << "CANNOT OPEN OUTPUT FILE!" << endl;
                    exit(-1);
               }
               h->deleteBySeq(seqOffset);
               seqNum ++;
               if(seqNum > adjustNum)
               {
                    h->adjust();
                    adjustNum = adjustNum + num / 20;
               }
//forcefully write this seq to avoid lost of it
//             numInCL = 1;

               clusterWithMismatches(sBuf, sQVBuf);
               clusterWithShifts(sBuf, sQVBuf);

//             if(numInCL > 100)
//             {
//                  dis << CLID << " " << numInCL << endl;
//                  big ++;
//             }
//             else if(numInCL == 1)
//             {
//                  small ++;
//             }
               CLID ++;
          }
//        cout << "#clusters of more than 100 seqs is " << big <<     endl;
//        cout << "#singleton clusters is " << small << endl;
     }
     else
     {
          for(seqOffset = 0; seqOffset < num * (upperSizeInBit + 5); seqOffset = seqOffset + (upperSizeInBit + 5))
          {
               realSize = h->searchBySeq(seqOffset, sBuf, seqID);
               if(realSize == 0) continue;
               else realSize = 0;
               calConsensus(sBuf, seqID, tagReverse);
               if(out.is_open())
                    if(reversed)
                         for(i = 0; i < mappingNum[seqID]; i ++)
                              out << CLID << "    " << mappingTable[seqID][i] << "   " << tagReverse << endl;
                    else
                         for(i = 0; i < mappingNum[seqID]; i ++)
                              out << CLID << "    " << mappingTable[seqID][i] << endl;
               else
               {
                    cout << "CANNOT OPEN OUTPUT FILE!" << endl;
                    exit(-1);
               }
               h->deleteBySeq(seqOffset);
               seqNum ++;
               if(seqNum > adjustNum)
               {
                    h->adjust();
                    adjustNum = adjustNum + num / 20;
               }
               clusterWithMismatches(sBuf);
               clusterWithShifts(sBuf);
               CLID ++;
          }
     }
     delete h;
     addiOut.close();
}

void Cluster::clusterWithMismatches(char sBuf[])
{
     char tBuf[1000];
     int indexNext, no, realSize, tagReverse, i;
     unsigned int indexOffset;
     unsigned int seqID;

     for(indexNext = 0; indexNext < seedsCount; indexNext ++)
     {
          indexOffset = h->calOffset(indexNext, sBuf);
          for(no = 0; no < h->getOffsetCount(indexOffset, indexNext); no ++)
          {
               realSize = h->searchByIndex(indexNext, indexOffset, no, tBuf, seqID);
               if(realSize && compare(sBuf, tBuf, 0, lowerSizeInChar, tagReverse) <= mismatchAllowed)
               {
                    h->deleteByIndex(indexNext, indexOffset, no);
                    if(out.is_open())
                         if(reversed)
                              for(i = 0; i < mappingNum[seqID]; i ++)
                                   out << CLID << "    " << mappingTable[seqID][i] << "   " << tagReverse << endl;
                         else
                              for(i = 0; i < mappingNum[seqID]; i ++)
                                   out << CLID << "    " << mappingTable[seqID][i] << endl;
                    else
                    {
                         cout << "CANNOT OPEN OUTPUT FILE!" << endl;
                         exit(-1);
                    }
                    seqNum ++;
//                  numInCL ++;
               }
          }
     }
}

void Cluster::clusterWithShifts(char sBuf[])
{
     char tBuf[1000], buf[1000];
     int no, i, lShift, rShift, realSize, indexNext, tagReverse;
     unsigned int indexOffset;
     unsigned int seqID;

     for(lShift = 1; lShift <= shiftAllowed; lShift ++)
     {
          for(i = 0; i < lowerSizeInChar - lShift; i ++)
               buf[i] = sBuf[i + lShift];
          for(indexNext = 0; indexNext < seedsCount; indexNext ++)
          {
               indexOffset = h->calOffset(indexNext, buf);
               for(no = 0; no < h->getOffsetCount(indexOffset, indexNext); no ++)
               {
                    realSize = (int)h->searchByIndex(indexNext, indexOffset, no, tBuf, seqID);
                    if(realSize && compare(buf, tBuf, 0, lowerSizeInChar - lShift, tagReverse) <= mismatchAllowed)//from 0 to 3 - lShift
                    {
                         h->deleteByIndex(indexNext, indexOffset, no);
                         if(out.is_open())
                              if(reversed)
                                   for(i = 0; i < mappingNum[seqID]; i ++)
                                        out << CLID << "    " << mappingTable[seqID][i] << "   " << tagReverse << endl;
                              else
                                   for(i = 0; i < mappingNum[seqID]; i ++)
                                        out << CLID << "    " << mappingTable[seqID][i] << endl;
                         else
                         {
                              cout << "CANNOT OPEN OUTPUT FILE!" << endl;
                              exit(-1);
                         }
//                       numInCL ++;
                         seqNum ++;
                    }
               }
          }
     }

     for(rShift = 1; rShift <= shiftAllowed; rShift ++)
     {
          for(i = rShift; i < lowerSizeInChar; i ++)
               buf[i] = sBuf[i - rShift];
          for(indexNext = 0; indexNext < seedsCount; indexNext ++)
          {
               indexOffset = h->calOffset(indexNext, buf);
               for(no = 0; no < h->getOffsetCount(indexOffset, indexNext); no ++)
               {
                    realSize = (int)h->searchByIndex(indexNext, indexOffset, no, tBuf, seqID);
                    if(realSize && compare(buf, tBuf, rShift, lowerSizeInChar, tagReverse) <= mismatchAllowed)//from rShift to 3
                    {
                         h->deleteByIndex(indexNext, indexOffset, no);
                         if(out.is_open())
                              if(reversed)
                                   for(i = 0; i < mappingNum[seqID]; i ++)
                                        out << CLID << "    " << mappingTable[seqID][i] << "   " << tagReverse << endl;
                              else
                                   for(i = 0; i < mappingNum[seqID]; i ++)
                                        out << CLID << "    " << mappingTable[seqID][i] << endl;
                         else
                         {
                              cout << "CANNOT OPEN OUTPUT FILE!" << endl;
                              exit(-1);
                         }
//                       numInCL ++;
                         seqNum ++;
                    }
               }
          }
     }
}

int Cluster::compare(char sBuf[], char tBuf[], int start, int end, int & tagReverse)
{
     int i, j, count = 0, reverseCount = 0;

     for(i = start; i < end; i ++)
          if(sBuf[i] != tBuf[i])
               count ++;

     if(reversed)
     {
          for(i = start, j = end - 1; i < end; i ++, j --)
               if(sBuf[i] != reverseChange(tBuf[j]))
                    reverseCount ++;
          if(count <= reverseCount) 
          {
               tagReverse = 0;
               return count;
          }
          else
          {
               tagReverse = 1;
               return reverseCount;
          }
     }
     else
          return count;
}

char Cluster::reverseChange(char base)
{
     switch(base)
     {
          case 0x00: return 0x03;
          case 0x01: return 0x02;
          case 0x02: return 0x01;
          case 0x03: return 0x00;
          default: cout << "MEMORY ERROR!!" << endl; exit(-1); 
     }
}

void Cluster::clusterWithMismatches(char sBuf[], char sQVBuf[])
{
     char tBuf[1000], tQVBuf[1000];
     int indexNext, no, realSize, tagReverse, i;
     unsigned int indexOffset;
     unsigned int seqID;

     for(indexNext = 0; indexNext < seedsCount; indexNext ++)
     {
          indexOffset = h->calOffset(indexNext, sBuf);
          for(no = 0; no < h->getOffsetCount(indexOffset, indexNext); no ++)
          {
               realSize = h->searchByIndex(indexNext, indexOffset, no, tBuf, tQVBuf, seqID);
               if(realSize && compare(sBuf, sQVBuf, tBuf, tQVBuf, 0, lowerSizeInChar, tagReverse) <= mismatchAllowed)
               {
                    h->deleteByIndex(indexNext, indexOffset, no);
                    if(out.is_open())
                         if(reversed)
                              for(i = 0; i < mappingNum[seqID]; i ++)
                                   out << CLID << "    " << mappingTable[seqID][i] << "   " << tagReverse << endl;
                         else
                              for(i = 0; i < mappingNum[seqID]; i ++)
                                   out << CLID << "    " << mappingTable[seqID][i] << endl;
                    else
                    {
                         cout << "CANNOT OPEN OUTPUT FILE!" << endl;
                         exit(-1);
                    }
//                  numInCL ++;
                    seqNum ++;
               }
          }
     }
}

void Cluster::clusterWithShifts(char sBuf[], char sQVBuf[])
{
     char tBuf[1000], buf[1000], tQVBuf[1000], QVBuf[1000];
     int no, i, lShift, rShift, realSize, indexNext, tagReverse;
     unsigned int indexOffset;
     unsigned int seqID;

     for(lShift = 1; lShift <= shiftAllowed; lShift ++)
     {
          for(i = 0; i < lowerSizeInChar - lShift; i ++)
          {
               buf[i] = sBuf[i + lShift];
               QVBuf[i] = sQVBuf[i + lShift];
          }
          for(indexNext = 0; indexNext < seedsCount; indexNext ++)
          {
               indexOffset = h->calOffset(indexNext, buf);
               for(no = 0; no < h->getOffsetCount(indexOffset, indexNext); no ++)
               {
                    realSize = (int)h->searchByIndex(indexNext, indexOffset, no, tBuf, tQVBuf, seqID);
                    if(realSize && compare(buf, QVBuf, tBuf, tQVBuf, 0, lowerSizeInChar - lShift, tagReverse) <= mismatchAllowed)//from 0 to 3 - lShift
                    {
                         h->deleteByIndex(indexNext, indexOffset, no);
                         if(out.is_open())
                              if(reversed)
                                   for(i = 0; i < mappingNum[seqID]; i ++)
                                        out << CLID << "    " << mappingTable[seqID][i] << "   " << tagReverse << endl;
                              else
                                   for(i = 0; i < mappingNum[seqID]; i ++)
                                        out << CLID << "    " << mappingTable[seqID][i] << endl;
                         else
                         {
                              cout << "CANNOT OPEN OUTPUT FILE!" << endl;
                              exit(-1);
                         }
//                       numInCL ++;
                         seqNum ++;
                    }
               }
          }
     }

     for(rShift = 1; rShift <= shiftAllowed; rShift ++)
     {
          for(i = rShift; i < lowerSizeInChar; i ++)
          {
               buf[i] = sBuf[i - rShift];
               QVBuf[i] = sQVBuf[i - rShift];
          }
          for(indexNext = 0; indexNext < seedsCount; indexNext ++)
          {
               indexOffset = h->calOffset(indexNext, buf);
               for(no = 0; no < h->getOffsetCount(indexOffset, indexNext); no ++)
               {
                    realSize = (int)h->searchByIndex(indexNext, indexOffset, no, tBuf, tQVBuf, seqID);
                    if(realSize && compare(buf, QVBuf, tBuf, tQVBuf, rShift, lowerSizeInChar, tagReverse) <= mismatchAllowed)//from rShift to 3
                    {
                         h->deleteByIndex(indexNext, indexOffset, no);
                         if(out.is_open())
                              if(reversed)
                                   for(i = 0; i < mappingNum[seqID]; i ++)
                                        out << CLID << "    " << mappingTable[seqID][i] << "   " << tagReverse << endl;
                              else
                                   for(i = 0; i < mappingNum[seqID]; i ++)
                                        out << CLID << "    " << mappingTable[seqID][i] << endl;
                         else
                         {
                              cout << "CANNOT OPEN OUTPUT FILE!" << endl;
                              exit(-1);
                         }
//                       numInCL ++;
                         seqNum ++;
                    }
               }
          }
     }
}

int Cluster::compare(char sBuf[], char sQVBuf[], char tBuf[], char tQVBuf[], int start, int end, int & tagReverse)
{
     int i, j, count = 0, qv = 0, reverseCount = 0, reverseQv = 0;

     for(i = start; i < end; i ++)
          if(sBuf[i] != tBuf[i] && !((sQVBuf[i] - OFFSET) + (tQVBuf[i] - OFFSET) < lowerQV))
          {
               count ++;
               qv = qv + (sQVBuf[i] - OFFSET) + (tQVBuf[i] - OFFSET);
          }

     if(reversed)
     {
          for(i = start, j = end - 1; i < end; i ++, j --)
               if(sBuf[i] != reverseChange(tBuf[j]) && !((sQVBuf[i] - OFFSET) + (tQVBuf[j] - OFFSET) < lowerQV))
               {
                    reverseCount ++;
                    reverseQv = reverseQv + (sQVBuf[i] - OFFSET) + (tQVBuf[j] - OFFSET);
               }
          if(count <= reverseCount)
          {
               tagReverse = 0;
               return qv > upperQV ? 10 : count;
          }
          else
          {
               tagReverse = 1;
               return reverseQv > upperQV ? 10 : reverseCount;
          }
     }
     else
     {
          return qv > upperQV ? 10 : count;
     }
}

void FileAnalyzer::outputAnalyze(int num, int correctCluster)
{
     ifstream in;
     char buf[20] = {0}, CLID[10] = {0}, seqID[10] = {0};
     int i = 0, j = 0, clusterID = -1, correctLower, correctUpper, subCluster = 0, multiCluster = 0;

     in.open("output.txt");
     if(in.is_open())
     {
          if(in.good())
          {
               in.getline(buf, 20);
               for(i = 0; i < 20; i ++)
                    buf[i] = 0;
               i = 0;
          }
          while(in.good())
          {
               in.getline(buf, 20);
               if(buf[0] == 0)
                    break;
               while(buf[i] != ' ')
                    CLID[i ++] = buf[i];
               i ++;
               while(buf[i] != 0)
                    seqID[j ++] = buf[i ++];

               if(atoi(CLID) > clusterID)
               {
                    clusterID = atoi(CLID);
                    correctLower = atoi(seqID) / (num / correctCluster) * (num / correctCluster);
                    correctUpper = correctLower + num / correctCluster - 1;
                    if(multiCluster == 0)
                         subCluster ++;
                    multiCluster = 0;
               }
               else
               {
                    if(atoi(seqID) < correctLower || atoi(seqID) > correctUpper)
                         multiCluster = 1;
               }
               for(i = 0; i < 10; i ++)
                    CLID[i] = seqID[i] = 0;
               for(i = 0; i < 20; i ++)
                    buf[i] = 0;
               i = j = 0;
          }
     }
     else
     {
          cout << "CANNOT OPEN INPUT FILE!" << endl;
          exit(-1);
     }
     cout << "(4) output analysis finished" << endl;
     cout << " - " << clusterID + 1 << " clusters in total" << endl;
     cout << " - " << subCluster << " sub clusters of degree " << (double)subCluster / (clusterID + 1) << endl;
}

void FileAnalyzer::inputAnalyze(char input[], int & num, int & tNum, int & lower, int & upper)
{
     ifstream in;
     int seqID = 0, i;
     char buf[1001];

     num = upper = 0;
     lower = 1000;
     in.open(input);

     if(in.is_open())
     {
//cont:
          while(in.good())
          {
               in.getline(buf, 1001);
               if(seqID % 4 == 1)
               {
//                  for(i = 0; i < in.gcount() - 1; i ++)
//                       if(buf[i] == 'N') 
//                       {
//                            seqID ++;
//                            goto cont;
//                       }
                    if(in.gcount() - 1 < lower) lower = in.gcount() - 1;
                    if(in.gcount() - 1 > upper) upper = in.gcount() - 1;
//                  num ++;
               }
               seqID ++;
          }
     }
     else
     {
          cout << "CANNOT OPEN INPUT FILE!" << endl;
          exit(-1);
     }
//   tNum = seqID;

//Filter reads based on first "lower" bases rather than all bases to avoid underestimate the number of valid reads
     in.clear(); in.seekg(0); seqID = 0; num = 0;
     if(in.is_open())
     {
conti:
          while(in.good())
          {
               in.getline(buf, 1001);
               if(seqID % 4 == 1)
               {
                    for(i = 0; i < lower; i ++)
                         if(buf[i] == 'N')
                         {
                              seqID ++;
                              goto conti;
                         }
                    num ++;
               }
               seqID ++;
          }
     }
     else
     {
          cout << "CANNOT OPEN INPUT FILE!" << endl;
          exit(-1);
     }
     tNum = seqID;
}

void FileAnalyzer::PECombine(char input1[], int lower1, char input2[], int lower2, 
char * input, int & num, int & lower, int & upper)
{
     ifstream in1, in2;
     ofstream out;
     char buf1[1000], buf2[1000];
     int seqID = 0, i, NBase;

     in1.open(input1);
     in2.open(input2);
     out.open("combined.fastq");

     strcpy(input, "combined.fastq");
     num = 0;
     lower = 1000;
     upper = 0;

     if(in1.is_open() && in2.is_open())
     {
          while(in1.good() && in2.good())
          {
               in1.getline(buf1, 1000);
               in2.getline(buf2, 1000);
               if(buf1[0] == 0 || buf2[0] == 0) break;

               if(seqID % 4 == 0)
                    out << "@" << seqID / 4 << endl;
               else if(seqID % 4 == 2)
                    out << "+" << seqID / 4 << endl;
               else
               {
                    NBase = 0;
                    for(i = 0; i < lower1; i ++)
                    {
                         out << buf1[i];
                         if(buf1[i] == 'N') NBase = 1;
                    }
                    for(i = 0; i < lower2; i ++)
                    {
                         out << buf2[i];
                         if(buf2[i] == 'N') NBase = 1;
                    }
                    out << endl;
                    if(NBase == 0) num ++;
                    if(in1.gcount() - 1 + in2.gcount() - 1 < lower) lower = in1.gcount() - 1 + in2.gcount() - 1;
                    if(in1.gcount() - 1 + in2.gcount() - 1 > upper) upper = in1.gcount() - 1 + in2.gcount() - 1;
               }
               seqID ++;
          }
     }
     else
     {
          cout << "CANNOT OPEN INPUT FILE!" << endl;
          exit(-1);
     }
}

FastqGenerator::FastqGenerator(char input[], char output[], int num)
{
        long i;

        in.open(input);
        strcpy(addiInput, output);
        strcat(addiInput, ".fasta");
        addiIn.open(addiInput);
        strcpy(outputq, output);
        strcat(outputq, ".fastq");
        out.open(outputq);
        this->num = num;
        seq = new char [num];
        for(i = 0; i < num; i ++)
                seq[i] = 0;
}

FastqGenerator::FastqGenerator(char input[], char output[], int num, int pair)
{
     long i;

     in.open(input);
     strcpy(addiInput, output);
     strcat(addiInput, ".fasta");
     addiIn.open(addiInput);
     strcpy(outputq, output);
     strcat(outputq, "."); 
     if(pair == 1) strcat(outputq, "1"); else strcat(outputq, "2");
     strcat(outputq, ".fastq");
     out.open(outputq);
     this->num = num;
     seq = new char [num];
     for(i = 0; i < num; i ++)
          seq[i] = 0;
}

void FastqGenerator::record()
{
     long i;
     char buf[1001];

     if(addiIn.is_open())
     {
          while(addiIn.good())
          {
               addiIn.getline(buf, 1001);
               if(buf[0] == 0)
                    break;
               for(i = 1; i < addiIn.gcount() - 1; i ++)
                    buf[i - 1] = buf[i];
               buf[i - 1] = '\0';
               seq[atoi(buf)] = 1;
               addiIn.getline(buf, 1001);
          }
     }
     else
     {
          cout << "CANNOT OPEN INPUT FILE!" << endl;
          exit(-1);
     }

//   for(i = 0; i < num; i ++)
//        if(seq[i])
//             cout << i << endl;
}

void FastqGenerator::generateFastq()
{
     unsigned long seqID, i;
     char buf[1001];

     record();
     if(in.is_open())
     {
          for(seqID = 0; seqID < num * 4; seqID ++)
          {
               in.getline(buf, 1001);
               if(seq[seqID / 4])
               {
                    for(i = 0; i < in.gcount() - 1; i ++)
                         out << buf[i];
                    out << endl;
               }
          }
     }
     else
     {
          cout << "CANNOT OPEN INPUT FILE!" << endl;
          exit(-1);
     }
}

Sorter::Sorter(char input[], int num, int lower)
{
     in.open(input);
     strcpy(midOutput, input);
     strcat(midOutput, ".mid.fastq");
     midOut.open(midOutput);
     this->num = num;
     mappingTable = new unsigned int * [num];
     mappingNum = new unsigned int [num];
     this->lowerSizeInChar = lower;
     this->realNum = 0;
}

int Sorter::getRealNum()
{
     return realNum;
}

unsigned int ** Sorter::getMappingTable()
{
     return mappingTable;
}

unsigned int * Sorter::getMappingNum()
{
     return mappingNum;
}

void Sorter::sort()
{
     char * seqs;
     Order * order;
     int i, j, seqID = 0, tag = 1, localNum = 0;
     char buf[1001];

     seqs = new char [(long int)num * lowerSizeInChar * 2];//2 means both bases and QVs
     order = new Order [num];

     if(in.is_open())
     {
cont:
          while(in.good())
          {
               in.getline(buf, 1001);
               if(seqID % 4 == 1)
                    for(i = 0; i < lowerSizeInChar; i ++)
                    {
                         if(buf[i] == 'N')
                         {
                              seqID ++;
                              tag = 0;
                              goto cont;
                         }
                         else
                              tag = 1;
                         seqs[(long int)localNum * lowerSizeInChar * 2 + i] = buf[i];
                    }
               if(seqID % 4 == 3 && tag == 1)
               {
                    for(i = lowerSizeInChar; i < lowerSizeInChar * 2; i ++)
                         seqs[(long int)localNum * lowerSizeInChar * 2 + i] = buf[i - lowerSizeInChar];
                    order[localNum ++].realID = seqID / 4;
               }
               seqID ++;
               if(localNum == num) break;
//Must finish to avoid crash. Otherwise:
//if there are reads with N in the end, they are not counted to initialize the seqs array but are put in seqs array
          }
     }
     else
     {
          cout << "CANNOT OPEN INPUT FILE!" << endl;
          exit(-1);
     }
     for(i = 0; i < localNum; i ++)
          order[i].ID = i;
//verification
//   cout << "-------------------------------------" << endl;
//   for(i = 0; i < localNum; i ++)
//        cout << order[i].ID << "|" << order[i].realID << " ";
//   cout << endl;
//   cout << "-------------------------------------" << endl;
//   for(i = 0; i < localNum; i ++)
//   {
//        for(j = 0; j < lowerSizeInChar * 2; j ++)
//             cout << seqs[i * lowerSizeInChar * 2 + j];
//        cout << endl;
//   }
//   cout << "-------------------------------------" << endl;
//verification

     suffixSort(0, localNum - 1, 0, seqs, order);
//verification
//   cout << "-------------------------------------" << endl;
//   for(i = 0; i < localNum; i ++)
//        cout << order[i].ID << "|" << order[i].realID << " ";
//   cout << endl;
//   cout << "-------------------------------------" << endl;
//   for(i = 0; i < realNum; i ++)
//   {
//        cout << mappingNum[i] << ": ";
//        for(j = 0; j < mappingNum[i]; j ++)
//             cout << mappingTable[i][j] << " ";
//        cout << endl;
//   }
//   cout << "-------------------------------------" << endl;
//verification
     delete seqs;
     delete order;
}

void Sorter::suffixSort(int start, int end, int depth, char seqs[], Order order[])
{
     int i, seqID = 0, tag;
     Order * buf;
     int s[4], e[4];
     int startBuf;
     int j;

     if(start == -1 && end == -1)
          return;

     if(start == end || depth == lowerSizeInChar - 1)
     {
          mappingTable[realNum] = new unsigned int [end - start + 1];
          for(i = start; i <= end; i ++)
               mappingTable[realNum][i - start] = order[i].realID;
          mappingNum[realNum] = end - start + 1;

          if(midOut.is_open())
          {
               midOut << "@" << realNum << endl;
               for(i = 0; i < lowerSizeInChar; i ++)
                    midOut << seqs[(long int)order[start].ID * lowerSizeInChar * 2 + i];
               midOut << endl;
               midOut << "+" << realNum << endl;
               for(i = lowerSizeInChar; i < lowerSizeInChar * 2; i ++)
                    midOut << seqs[(long int)order[start].ID * lowerSizeInChar * 2 + i];
               midOut << endl;
          }
          else
          {
               cout << "CANNOT OPEN OUTPUT FILE!" << endl;
               exit(-1);
          }

          realNum ++;
          return;
     }

     buf = new Order [end - start + 1];

     startBuf = start;
     tag = 0;
     for(i = start; i <= end; i ++)
          if(seqs[(long int)order[i].ID * lowerSizeInChar * 2 + depth] == 'A')
          {
               tag = 1;
               buf[seqID].ID = order[i].ID;
               buf[seqID ++].realID = order[i].realID;
          }
     if(tag == 1)
     {
          s[0] = startBuf;
          e[0] = start + seqID - 1;
          startBuf = start + seqID;
     }
     else
          s[0] = e[0] = -1;

     tag = 0;
     for(i = start; i <= end; i ++)
          if(seqs[(long int)order[i].ID * lowerSizeInChar * 2 + depth] == 'C')
          {
               tag = 1;
               buf[seqID].ID = order[i].ID;
               buf[seqID ++].realID = order[i].realID;
          }
     if(tag == 1)
     {
          s[1] = startBuf;
          e[1] = start + seqID - 1;
          startBuf = start + seqID;
     }
     else
          s[1] = e[1] = -1;

     tag = 0;
     for(i = start; i <= end; i ++)
          if(seqs[(long int)order[i].ID * lowerSizeInChar * 2 + depth] == 'G')
          {
               tag = 1;
               buf[seqID].ID = order[i].ID;
               buf[seqID ++].realID = order[i].realID;
          }
     if(tag == 1)
     {
          s[2] = startBuf;
          e[2] = start + seqID - 1;
          startBuf = start + seqID;
     }
     else
          s[2] = e[2] = -1;

     tag = 0;
     for(i = start; i <= end; i ++)
          if(seqs[(long int)order[i].ID * lowerSizeInChar * 2 + depth] == 'T')
          {
               tag = 1;
               buf[seqID].ID = order[i].ID;
               buf[seqID ++].realID = order[i].realID;
          }
     if(tag == 1)
     {
          s[3] = startBuf;
          e[3] = start + seqID - 1;
          startBuf = start + seqID;
     }
     else
          s[3] = e[3] = -1;

     for(i = start, seqID = 0; i <= end; i ++, seqID ++)
     {
          order[i].ID = buf[seqID].ID;
          order[i].realID = buf[seqID].realID;
     }

//   cout << s[0] << ", " << s[1] << ", " << s[2] << ", " << s[3] << endl;
//   cout << e[0] << ", " << e[1] << ", " << e[2] << ", " << e[3] << endl;

     suffixSort(s[0], e[0], depth + 1, seqs, order);
     suffixSort(s[1], e[1], depth + 1, seqs, order);
     suffixSort(s[2], e[2], depth + 1, seqs, order);
     suffixSort(s[3], e[3], depth + 1, seqs, order);
}

char change(char base)
{
     switch(base)
     {
          case 0x00: return 'A';
          case 0x01: return 'C';
          case 0x02: return 'G';
          case 0x03: return 'T';
          default: cout << "UNKNOWN ERROR!" << endl; exit(-1);
     }
}

bool within(int p, int pos[], int size)
{
     int i;
     for(i = 0; i < size; i ++)
          if(p == pos[i]) return true;
     return false;
}

void introduceMismatches(char sBuf[], char buf[], int mismatch, int size)
{
     int i, j;
     int pos[1000], mBuf[1000] = {0};

     for(i = 0; i < mismatch; i ++)
     {
          do
               pos[i] = rand() % size;
          while(mBuf[pos[i]] == 1);
          mBuf[pos[i]] = 1;
     }
     for(i = 0, j = 0; i < size; i ++)
     {
          if(within(i, pos, mismatch) && j < mismatch)
          {
               do
                    buf[i] = change((char)(rand() % 4));
               while(sBuf[i] == buf[i]);
               j ++;
          }
          else
               buf[i] = sBuf[i];
     }
}

void introduceShifts(char sBuf[], char buf[], int shift, int size)
{
     int i;

     if(shift < 0)
     {
          for(i = 0; i < size - abs(shift); i ++)
               buf[i] = sBuf[i + abs(shift)];
          for(i = size - abs(shift); i < size; i ++)
               buf[i] = change((char)(rand() % 4));
     }
     else
     {
          for(i = 0; i < shift; i ++)
               buf[i] = change((char)(rand() % 4));
          for(i = shift; i < size; i ++)
               buf[i] = sBuf[i - shift];
     }
}

#ifdef WITHSIMILARITY
void generateClusteredSeq(int num, int lower, int upper, int mismatchAllowed, int shiftAllowed, int correctCluster, int distance)
#else
void generateClusteredSeq(int num, int lower, int upper, int mismatchAllowed, int shiftAllowed, int correctCluster)
#endif
{
     ofstream out;
     int size, mismatch, shift, i, j, k;
     char s[1000], sBuf[1000], mismatchBuf[1000], shiftBuf[1000], buf[1000];
     int mBuf[1000] = {0};

#ifdef WITHSIMILARITY
     if(distance < 2) 
     {
          cout << "incorrect distance" << endl;
          return;
     }
     for(i = 0; i < 1000; i ++)
          s[i] = change((char)(rand() % 4));
#endif
     out.open("input.txt");

     for(k = 0; k < correctCluster; k ++)
     {
#ifdef WITHSIMILARITY
          generateCenter(s, sBuf, distance, mBuf);
#else
          for(i = 0; i < 1000; i ++)
               sBuf[i] = change((char)(rand() % 4));
#endif
          for(i = 0; i < num / correctCluster; i ++)
          {
               size = lower + rand() % (upper - lower + 1);
               mismatch = rand() % (mismatchAllowed + 1);
               introduceMismatches(sBuf, mismatchBuf, mismatch, size);
               shift = rand() % (shiftAllowed * 2 + 1) - shiftAllowed;
               introduceShifts(mismatchBuf, shiftBuf, shift, size);
//             if(rand() % 2 == 1)
//                  introduceMismatches(sBuf, buf, mismatch, size);
//             else
//                  introduceShifts(sBuf, buf, shift, size);
               out << "@title " << k * num / correctCluster + i << " size = " << size << " mismatches = " << mismatch << " shifts = " << shift << endl;
               out.write(shiftBuf, size);
//             out.write(buf, size);
               out << endl;
               out << "+title " << k * num / correctCluster + i << " size = " << size << " mismatches = " << mismatch << " shifts = " << shift << endl;
               for(j = 0; j < size; j ++)
                    shiftBuf[j] = (char)(rand() % RANGE) + OFFSET;
               out.write(shiftBuf, size);
               out << endl;
          }
     }
}

void itoa(char buf[], unsigned int v)
{
     if(v / 1000)
     {
          buf[0] = v / 1000 + 48;
          buf[1] = (v % 1000) / 100 + 48;
          buf[2] = (v % 100) / 10 + 48;
          buf[3] = v % 10 + 48;
          buf[4] = '\0';
     }
     else if(v / 100)
     {
          buf[0] = v / 100 + 48;
          buf[1] = (v % 100) / 10 + 48;
          buf[2] = v % 10 + 48;
          buf[3] = '\0';
     }
     else if(v / 10)
     {
          buf[0] = v / 10 + 48;
          buf[1] = v % 10 + 48;
          buf[2] = '\0';
     }
     else
     {
          buf[0] = v + 48;
          buf[1] = '\0';
     }
}

void print()
{

     cout << "Usage:" << endl;
     cout << "\t./seed --input <FILE> --output <FILE> [options]*" << endl;
     cout << "Options:" << endl;
     cout << "\t--M1 <int>\tmax # mismatches in first seed alignment (0 - 3, default 3)" << endl;
     cout << "\t--S1 <int>\tmax # shifts in first seed alignment (0 - 6, default 3)" << endl;
     cout << "\t--M2 <int>\tmax # mismatches in first seed alignment (0 - 20, default 10)" << endl;
     cout << "\t--S2 <int>\tmax # shifts in first seed alignment (0 - 6, default 3)" << endl;
     cout << "\t--QV1    \tthreshold for the base call quality values (QV)" << endl;
     cout << "\t--QV2    \tanother QV threshold" << endl;
     cout << "\t--default\trun in default mode" << endl;
     cout << "\t--fast   \trun in fast mode" << endl;
     cout << "\t--short  \trun in short mode" << endl;
     cout << "\t--reverse\t" << endl;
     cout << "\t--input2 <FILE>\tspecifies the paired sequences" << endl;
     cout << "\t--KM1  \tuse 1st mapping strategy for K-means part (default choice)" << endl;
     cout << "\t--KM2  \tuse 2nd mapping strategy for K-means part" << endl;

     //cout << "\t--L <int>\tuse 1st mapping algorithm" << endl;
     //cout << "\t--N <int>\tuse 2nd mapping algorithm" << endl;
     cout << "Without options, seed will run in default mode (see more in the manual)" << endl;

}

int old_seed(int argc, char * argv[])
{
     time_t start, end;
     int num, num1, num2, tNum, tNum1, tNum2, lower, lower1, lower2, upper, upper1, upper2, mismatch = 3, shift = 3, lowerQV = 0, upperQV = 6 * 93, i, tagMismatch = 0, tagShift = 0, tagInput = 0, tagInput2 = 0, tagOutput = 0, tagFast = 0, tagShort = 0, tagReverse = 0;
     int tagQV1 = 0, tagQV2 = 0;
     char buf[5], input[100], input1[100], input2[100], output[100], midOutput[100];
     ifstream in, in2;
     int io = 0;
     int totalLength, count, j;

     for(i = 1; i < argc; i ++)
          if(strcmp(argv[i], "--input") == 0)
          {
               if(tagInput == 1 || i == argc - 1)
               {
                    print();
                    return 0;
               }
               in.open(argv[++ i]);
               if(!in.is_open())
               {
                    cout << "CANNOT OPEN INPUT FILE!" << endl;
                    print();
                    return 0;
               }
               in.close();
               strcpy(input, argv[i]);
               tagInput = 1;
          }
                else if(strcmp(argv[i], "--input2") == 0)
                {
                        if(tagInput2 == 1 || i == argc - 1)
                        {
                                print();
                                return 0;
                        }
                        in2.open(argv[++ i]);
                        if(!in2.is_open())
                        {
                                cout << "CANNOT OPEN PAIRED INPUT FILE!" << endl;
                                print();
                                return 0;
                        }
                        in2.close();
                        strcpy(input2, argv[i]);
                        tagInput2 = 1;
               paired = 1;
                }
          else if(strcmp(argv[i], "--output") == 0)
          {
               if(tagOutput == 1 || i == argc - 1)
               {
                    print();
                    return 0;
               }
               strcpy(output, argv[++ i]);
               tagOutput = 1;
          }
          else if(strcmp(argv[i], "--mismatch") == 0)
          {
               if(tagMismatch == 1 || i == argc - 1)
               {
                    print();
                    return 0;
               }
               mismatch = atoi(argv[++ i]);
               itoa(buf, mismatch);
               if(strcmp(argv[i], buf) != 0)
               {
                    print();
                    return 0;
               }
               tagMismatch = 1;
          }
          else if(strcmp(argv[i], "--shift") == 0)
          {
               if(tagShift == 1 || i == argc - 1)
               {
                    print();
                    return 0;
               }
               shift = atoi(argv[++ i]);
               itoa(buf, shift);
               if(strcmp(argv[i], buf) != 0)
               {
                    print();
                    return 0;
               }
               tagShift = 1;
          }
          else if(strcmp(argv[i], "--QV1") == 0)
          {
               if(tagQV1 == 1 || i == argc - 1)
               {
                    print();
                    return 0;
               }
               lowerQV = atoi(argv[++ i]);
               itoa(buf, lowerQV);
               if(strcmp(argv[i], buf) != 0)
               {
                    print();
                    return 0;
               }
               tagQV1 = 1;
               QV = 1;
          }
          else if(strcmp(argv[i], "--QV2") == 0)
          {
               if(tagQV2 == 1 || i == argc - 1)
               {
                    print();
                    return 0;
               }
               upperQV = atoi(argv[++ i]);
               itoa(buf, upperQV);
               if(strcmp(argv[i], buf) != 0)
               {
                    print();
                    return 0;
               }
               tagQV2 = 1;
               QV = 1;
          }
          else if(strcmp(argv[i], "--fast") == 0)
          {
               if(tagFast == 1 || tagShort == 1)
               {
                    print();
                    return 0;
               }
               tagFast = 1;
               seedsCount = 4;
               seedsWeight = 64 * 1024;
          }
          else if(strcmp(argv[i], "--short") == 0)
                {
                        if(tagFast == 1 || tagShort == 1)
                        {
                                print();
                                return 0;
                        }
                        tagShort = 1;
                        seedsWeight = 4;
                }
          else if(strcmp(argv[i], "--reverse") == 0)
          {
               if(tagReverse == 1)
               {
                    print();
                    return 0;
               }
               tagReverse = 1;
               reversed = 1;
          }
          else
          {
               print();
               return 0;
          }

     if(tagInput == 0 || tagOutput == 0 || mismatch < 0 || shift < 0 || lowerQV < 0 || lowerQV > 2 * 93 || upperQV < 0 || upperQV > 6 * 93)
     {
          print();
          return 0;
     }

     if(QV)
          cout << "#mismatch = " << mismatch << "; #shift = " << shift << "; QV1 = " << lowerQV << "; QV2 = " << upperQV << endl;
     else
          cout << "#mismatch = " << mismatch << "; #shift = " << shift << endl;

//   generateClusteredSeq(1000, 95, 100, 0, 0, 100);
//   return 0;

     FileAnalyzer fa;
     start = time(NULL);
     if(paired == 0)
     {
          fa.inputAnalyze(input, num, tNum, lower, upper);

          cout << "(1) input analysis finished" << endl;
          cout << " - " << num << " valid seqs with lengths between " << lower << " and " << upper << endl;

          if(num == 0)
          {
               cout << "INSUFFICIENT VALID READS!" << endl;
               return 0;
          }
          if(upper - lower > LENGTH_DIFFERENCE)
          {
               cout << "INVALID READ LENGTH DIFFERENCE (ABOVE 5)!" << endl;
               return 0;
          }
          if(lower < 36 && seedsWeight == 1024 * 16)
          {
               cout << "INVALID READ LENGTH (BELOW 36) IN ORDINARY MODE!" << endl;
               return 0;
          }
          if(lower < 58 && seedsWeight == 1024 * 64)
          {
               cout << "INVALID READ LENGTH (BELOW 58) IN FAST MODE!" << endl;
               return 0;
          }
          if(lower < 21 && seedsWeight == 4)
          {
               cout << "INVALID READ LENGTH (BELOW 21) IN SHORT MODE!" << endl;
               return 0;
          }
          if(upper > 1000)
          {
               cout << "INVALID READ LENGTH (ABOVE 1000)!" << endl;
               return 0;
          }
     }
     else
     {
          strcpy(input1, input);
          fa.inputAnalyze(input1, num1, tNum1, lower1, upper1);
          paired = lower1;//keep lower1 in paired to separate read pairs
          fa.inputAnalyze(input2, num2, tNum2, lower2, upper2);
          fa.PECombine(input1, lower1, input2, lower2, input, num, lower, upper);
          //combine both pairs and trim to keep reads in the same pair same length

          cout << "(1) input analysis finished" << endl;
//        cout << " - " << num1 << " valid seqs with lengths between " << lower1 << " and " << upper1 << " in left pair" << endl;
//        cout << " - " << num2 << " valid seqs with lengths between " << lower2 << " and " << upper2 << " in right pair" << endl;
          cout << " - " << num << " valid seqs with combined lengths between " << lower << " and " << upper << endl;

          if(num1 == 0)
          {
               cout << "INSUFFICIENT VALID READS IN LEFT PAIR!" << endl;
               return 0;
          }
          if(num2 == 0)
          {
               cout << "INSUFFICIENT VALID READS IN RIGHT PAIR!" << endl;
               return 0;
          }
          if(tNum1 != tNum2)
          {
               cout << "DIFFERENT NUMBER OF READS IN LEFT AND RIGHT PAIRS!" << endl;
               return 0;
          }
          if(upper1 - lower1 > 5)
          {
               cout << "INVALID READ LENGTH DIFFERENCE (ABOVE 5) IN LEFT PAIR!" << endl;
               return 0;
          }
          if(upper2 - lower2 > 5)
          {
               cout << "INVALID READ LENGTH DIFFERENCE (ABOVE 5) IN RIGHT PAIR!" << endl;
               return 0;
          }
          if(lower < 36 && seedsWeight == 1024 * 16)
          {
               cout << "INVALID COMBINED READ LENGTH (BELOW 36) IN ORDINARY MODE!" << endl;
               return 0;
          }
          if(lower < 58 && seedsWeight == 1024 * 64)
          {
               cout << "INVALID COMBINED READ LENGTH (BELOW 58) IN FAST MODE!" << endl;
               return 0;
          }
          if(lower < 21 && seedsWeight == 4)
          {
               cout << "INVALID COMBINED READ LENGTH (BELOW 21) IN SHORT MODE!" << endl;
               return 0;
          }
          if(upper > 1000)
          {
               cout << "INVALID COMBINED READ LENGTH (ABOVE 1000)!" << endl;
               return 0;
          }

          if(shift > 0)
          {
               cout << "In current implementation, #shift must be 0 for paired-end clustering. Please wait for SEED2 to solve this issue." << endl;
               return 0;
          }
          if(reversed == 1 && lower1 != lower2)
          {
               cout << "In current implementation, if reverse complementary is considered in clustering, lower bounds of both pairs should be the same. Please wait for SEED2 to solve this issue." << endl;
               return 0;
          }

          upper = lower;
     }

//   produce realNum, mappingTable and mappingNum here, and the intermediate file is produced/opened by protocol
     Sorter s(input, num, lower);
     s.sort();// if seqs x and y are the same but with differnt QVs, then y's QV will be represented by x's QV and not be considered in clustering
     cout << "(2) sorting finished" << endl;

     Cluster c(input, output, s.getRealNum(), lower, upper, mismatch, shift, lowerQV, upperQV, s.getMappingTable(), s.getMappingNum());
     cout << "(3) init finished" << endl;

     c.cluster();
     end = time(NULL);
     cout << "(4) clustering finished" << endl;

     if(paired == 0)
     {
          FastqGenerator f(input, output, tNum);
          f.generateFastq();
     }
     else
     {
          FastqGenerator f1(input1, output, tNum1, 1);
          f1.generateFastq();
          FastqGenerator f2(input2, output, tNum2, 2);
          f2.generateFastq();
     }
     cout << "(5) fastq file generated" << endl;

     cout << " - " << end - start << " seconds" << endl;

}

int f_gcd(int l, int w){
     
     while(l % w != 0){
          l = l-w;

          if(l<w){
               int k=l;
               l=w;
               w=k;
          }
     }

     return w;

}

int factorial(int num){
     int result=1;

     if (num == 0)
     {
          return 1;
     }

     while(num != 0){
          result *= num;
          num--;
     }

     return result;
}

int combination(int n,int m){

     if (m==0 || m==n)
     {
          return 1;
     } else if (m==1 || m==n-1)
     {
          return n;
     } 

     if (n==2)
     {
          return factorial(n) / (factorial(m) * factorial(n-m));
     } 
     else
     {
          return combination(n-1,m-1) + combination(n-1,m);
     }
}

void set_seeds(int l,int w){

     seeds_length = l;
     seeds_weight = w;

     int gcd = f_gcd(l,w);
     int num_block_seed = l / gcd;
     int num_block_one = w / gcd;
     int num_block_zero = (l-w) / gcd;

     int num_seeds = combination(num_block_seed,num_block_one);
     //int num_seeds = factorial(num_block_seed) / (factorial(num_block_one) * factorial(num_block_zero));

     seedsCount = num_seeds;

     int length_block = w / num_block_one;
     


     seeds = new int*[num_seeds];
     for(int i=0; i<num_seeds; i++){
          seeds[i] = new int[length_block * num_block_seed]; // 6 * 5
     }

     int *origin_seed = new int[num_block_seed];

     int index = 0;
     for(; index < num_block_one; index++){
          origin_seed[index] = 1;
     }
     // 如果数组值为0的话，next_permutation函数会忽略之，比如 '11000' 只会有一个排列输出
     for(; index < num_block_seed; index++){
          origin_seed[index] = 2;
     }


     index = 0;
     do{
          for(int i=0; i<num_block_seed; i++){

               if(origin_seed[i] == 2){
                    for(int j=0; j<length_block; j++){
                         seeds[index][i*length_block+j] = 0;
                    }
               } else{
                    for(int j=0; j<length_block; j++){
                         seeds[index][i*length_block+j] = 1;
                    }
               }

          }
          index++;
     }while(next_permutation(origin_seed,origin_seed+num_block_seed));

     // 显示 seed 数组结果
     /*
     for (int i=0; i<num_seeds; i++)
     {
          for (int j=0; j<length_block * num_block_seed; j++)
          {
               cout << seeds[i][j] << " ";
          }
          cout << endl;
     }*/


}

multimap<int,string> K_Means_2(multimap<int,string> clid_seq,int expectation_int){

     map<int, string>::iterator iter;
     int data_size = clid_seq[0]->second.length()  *  clid_seq.size(); 
     /*
     for(iter = clid_seq.begin(); iter !=clid_seq.end(); iter++){
          string seq = iter->second;
          if (seq.length() > data_size){
               data_size = seq.length();
          }
     }*/

     double *data = new double[data_size];

     for (int i=0; i<data_size; i++){
          data[i] = 0.0;
     }


     int counter_data = 0;
     // A-1 C-2 G-3 T-4
     for(iter = clid_seq.begin(); iter !=clid_seq.end(); iter++){
          string seq = iter->second;
          for (int i=0; i<seq.size();i++)
          {
               if (seq[i] == 'A')
               {
                    data[counter_data++] = 1;
               } else if(seq[i] == 'C'){
                    data[counter_data++] = 2;
               } else if(seq[i] == 'G'){
                    data[counter_data++] = 3;
               } else if(seq[i] == 'T'){
                    data[counter_data++] = 4;
               }
          }

     }

     const int size = clid_seq.size(); //Number of samples
     const int dim = clid_seq[0]->second.length();   //Dimension of feature
     const int cluster_num = clid_seq.size() / expectation_int + 1; //Cluster number

     KMeans* kmeans = new KMeans(dim,cluster_num);
     int* labels = new int[size];
     kmeans->SetInitMode(KMeans::InitRandom);
     kmeans->Cluster(data,size,labels);

     ofstream ofile;
     ofile.open("k-means.txt",ios::app);

     multimap<int,string> clid_seq_kmeans;
     int i=0;
     ofile << "CLID\t" << clid_seq.begin()->first << endl;
     for(iter = clid_seq.begin(); iter !=clid_seq.end(); iter++){
          string seq = iter->second;
          ofile << seq << endl;
          ofile << labels[i] << "\t["<< data[i*dim+0] << ','<< data[i*dim+1]<< ','<< data[i*dim+2]<< ','<< data[i*dim+3] << "]" <<endl;

          clid_seq_kmeans.insert(pair<int,string>(labels[i],seq));

          i++;
     }

     clid_seq.clear();
     delete []labels;
     delete kmeans;

     return clid_seq_kmeans;
}

multimap<int,string> K_Means_1(multimap<int,string> clid_seq,int expectation_int){

     map<int, string>::iterator iter;

     int data_size = clid_seq.size()*4; // seq size * 4(acgt)
     double *data = new double[data_size];
     int A_num=0,G_num=0,C_num=0,T_num = 0;
     double A_ratio,G_ratio,C_ratio,T_ratio;

     int counter_data = 0;

     for(iter = clid_seq.begin(); iter !=clid_seq.end(); iter++){
          string seq = iter->second;
          for (int i=0; i<seq.size();i++)
          {
               if (seq[i] == 'A')
               {
                    A_num++;
               } else if(seq[i] == 'C'){
                    C_num++;
               } else if(seq[i] == 'G'){
                    G_num++;
               } else if(seq[i] == 'T'){
                    T_num++;
               }
          }
          A_ratio = (double) A_num / seq.size();
          C_ratio = (double) C_num / seq.size();
          G_ratio = (double) G_num / seq.size();
          T_ratio = (double) T_num / seq.size();

          data[counter_data++] = A_ratio;
          data[counter_data++] = C_ratio;
          data[counter_data++] = G_ratio;
          data[counter_data++] = T_ratio;

          A_num = 0;
          C_num = 0;
          G_num = 0;
          T_num = 0;

     }

     const int size = clid_seq.size(); //Number of samples
     const int dim = 4;   //Dimension of feature
     const int cluster_num = clid_seq.size() / expectation_int + 1; //Cluster number

     KMeans* kmeans = new KMeans(dim,cluster_num);
     int* labels = new int[size];
     kmeans->SetInitMode(KMeans::InitRandom);
     kmeans->Cluster(data,size,labels);

     ofstream ofile;
     ofile.open("k-means.txt",ios::app);

     multimap<int,string> clid_seq_kmeans;
     int i=0;
     ofile << "CLID\t" << clid_seq.begin()->first << endl;
     for(iter = clid_seq.begin(); iter !=clid_seq.end(); iter++){
          string seq = iter->second;
          ofile << seq << endl;
          ofile << labels[i] << "\t["<< data[i*dim+0] << ','<< data[i*dim+1]<< ','<< data[i*dim+2]<< ','<< data[i*dim+3] << "]" <<endl;

          clid_seq_kmeans.insert(pair<int,string>(labels[i],seq));

          i++;
     }

     clid_seq.clear();
     delete []labels;
     delete kmeans;

     return clid_seq_kmeans;
}





/************************************************************************/
/* argv[]
   0 ...\\Seed.exe
   1 --input
   2 input.fastq
   3 --output
   4 output.txt
   5 --shift
   6 3
   */
/************************************************************************/
int main(int argc, char * argv[])
{
     string input,output;

     int tag_input=0,tag_output=0,tag_mismatch1=0,tag_shift1=0,tag_mismatch2=0,tag_shift2=0;
     int tag_L = 0,tag_N = 0,tag_fast=0,tag_short=0,tag_default=0;
     int tag_KM1 = 0,tag_KM2 = 0,KM_flag = 1;
     int argv_mismatch1=3,argv_shift1=3,argv_mismatch2=10,argv_shift2=3,argv_L=30,argv_N=12;

     for(int i=1; i<argc; i++){
          if(strcmp(argv[i], "--input") == 0){
               if(tag_input > 0 || i == argc - 1)
               {
                    print();
                    return 0;
               }
               ifstream in;
               in.open(argv[i+1]);
               if(!in.is_open())
               {
                    cout << "CANNOT OPEN INPUT FILE!" << endl;
                    print();
                    return 0;
               }
               in.close();
               input = argv[i+1];
               tag_input = i;
               i++;

          } 
          else if(strcmp(argv[i], "--output") == 0){
               if(tag_output > 0 || i == argc - 1)
               {
                    print();
                    return 0;
               }

               output = argv[i+1];
               tag_output = i;
               i++;

          } 
          else if(strcmp(argv[i], "--M1") == 0){
               if(tag_mismatch1 > 0 || i == argc - 1)
               {
                    print();
                    return 0;
               }
               argv_mismatch1 = atoi(argv[i+1]);
               tag_mismatch1 = i;
               i++;

          } 
          else if(strcmp(argv[i], "--S1") == 0){
               if(tag_shift1 > 0|| i == argc - 1)
               {
                    print();
                    return 0;
               }
               argv_shift1 = atoi(argv[i+1]);
               tag_shift1 = i;
               i++;

          } 
          else if(strcmp(argv[i], "--M2") == 0){
               if(tag_mismatch2 > 0 || i == argc - 1)
               {
                    print();
                    return 0;
               }

               argv_mismatch2 = atoi(argv[i+1]);
               tag_mismatch2 = i;
               i++;

          }
          else if(strcmp(argv[i], "--S2") == 0){
               if(tag_shift2 > 0 || i == argc - 1)
               {
                    print();
                    return 0;
               }

               argv_shift2 = atoi(argv[i+1]);
               tag_shift2 = i;
          }
          else if(strcmp(argv[i], "--fast") == 0){
               if(tag_fast == 1 || i == argc - 1)
               {
                    print();
                    return 0;
               }

               seeds_length = 52;
               seeds_weight = 13;

               tag_fast = i;
          }
          else if(strcmp(argv[i], "--short") == 0){
               if(tag_short > 0 || i == argc - 1)
               {
                    print();
                    return 0;
               }

               seeds_length = 15;
               seeds_weight = 6;

               tag_short = i;
          }
          else if(strcmp(argv[i], "--default") == 0){
               if(tag_default > 0 || i == argc - 1)
               {
                    print();
                    return 0;
               }

               seeds_length = 30;
               seeds_weight = 12;

               tag_default = i;
          }
          else if(strcmp(argv[i], "--KM1") == 0){
               if(tag_KM1 > 0 || i == argc - 1)
               {
                    print();
                    return 0;
               }

               tag_KM1 = 1;
               KM_flag = 1;
          }
          else if(strcmp(argv[i], "--KM2") == 0){
               if(tag_KM2 > 0 || i == argc - 1)
               {
                    print();
                    return 0;
               }

               tag_KM2 = 1;
               KM_flag = 2;
          }
     

     }

     if(tag_input == 0 || tag_output==0){
          print();
          return 0;
     }



     if(tag_KM2 == 1 && tag_KM1 ==1){
          cout << "KM1 and KM2 cannot be set at the same time" << endl;
          print();
          return 0;
     } 

     if((tag_fast != 0 && tag_short != 0) || (tag_default != 0 && tag_short != 0) || (tag_fast != 0 && tag_default != 0) ){
          cout << "default mode, fast mode and short mode cannot be set at the same time" << endl;
          print();
          return 0;
     }

     if(tag_L == 1 && tag_N ==1){
          if(tag_fast > 0){
               strcpy(argv[tag_fast],"");
          }

          if(tag_short > 0){
               strcpy(argv[tag_short],"");
          }

          if(tag_default > 0){
               strcpy(argv[tag_default],"");
          }

          seeds_length = argv_L;
          seeds_weight = argv_N;
          //GXX
          //seedsWeight = 

     }

     int new_argc = argc + 4;
     char ** new_argv = new char*[new_argc];

     for (int i=0; i<new_argc; i++)
     {
          new_argv[i] = new char[20];
          
     }
     //the argv[0] char too long
     for (int i=1; i<new_argc; i++)
     {
          strcpy(new_argv[i],"");
     }


     for (int i=1; i<argc; i++)
     {
          strcpy(new_argv[i],argv[i]);
     }

     if (tag_mismatch1 == 0)
     {
          strcpy(new_argv[new_argc-4],"--mismatch");
          itoa(argv_mismatch1, new_argv[new_argc-3], 10);
     } 

     if (tag_shift1 == 0)
     {
          strcpy(new_argv[new_argc-2],"--shift");
          itoa(argv_shift1, new_argv[new_argc-1], 10);                          
     }

     //return 0;

     string input_firsttime_kmeans = input;

     ofstream statis_data;
     statis_data.open("statictic_data.txt",ios::out|ios::trunc);
     ofstream input_from_output;
     input_from_output.open("input_from_output.fastq",ios::out|ios::trunc);

     ofile.open ("input2.fastq", ios::out|ios::trunc);
     int num_cluster=0;
     int num_per_cluster=0;

     multimap<int,int> cluster_seed_output;   // 将output.txt 读入   left: Cluster 编号  right：原有sequence id    
     multimap<int,int> cluster_seed_output3,cluster_seed_output2;
     
     cout << "Running Seed For the First Time" << endl;
     set_seeds(seeds_length,seeds_weight);


     for (int i=0; i<argc; i++)
     {
          cout << argv[i] << endl;
     }

     for (int i=0; i<new_argc; i++)
     {
          cout << new_argv[i] << endl;
     }

     old_seed(new_argc,new_argv);
     // 避免 input2.fastq 写两遍
     ofile.open ("input3.fastq", ios::out|ios::trunc);
     // statictic for fisrt seed
     statis_data << "First Seed Result" << endl;



     string line;
     ifstream ifile (output);
     if (ifile.is_open())
     {
          getline (ifile,line); // 读掉首行
          while ( getline (ifile,line) )
          {

               if (line[0]>= 'A' && line[0] <='T')
               {
                    num_cluster++;
                    //ofile_input1_info << "cluster " << num_cluster;
                    input_from_output << "@1" << endl << line << endl << "+" << endl << endl;
               }
               else if(line[0]>='0' && line[0] <= '9'){

                    int index_space = line.find_first_of('\t');
                    string left_str = line.substr(0,index_space);
                    string right_str = line.substr(index_space+1,line.length()-index_space);
                    int left = stoi(left_str);
                    int right = stoi(right_str);
                    cluster_seed_output.insert(pair<int,int>(left,right));
               }
          }
          ifile.close();
     }else {
          cout << "Unable to open file"; 
     }
     multimap<int,int> cluster_seed_size_output;  // <size,cluster id>

     for (int i=0; i<cluster_seed_output.size();i++)
     {
          int num = cluster_seed_output.count(i);

          if (num > 0)
          {
               cluster_seed_size_output.insert(pair<int,int>(num,i));
          }
     }

     statis_data << "MISMATCH\t" << argv[8] << "\tSHIFT\t" << argv[6] << endl;

     for (int i=1; i<=cluster_seed_size_output.size();i++)
     {
          int num = cluster_seed_size_output.count(i);

          if (num > 0)
          {
               statis_data << "size\t" << i << "\tnum\t" << num << endl;
          }
     }

     // continue for seed2 or stop at seed
     
     char continue_flag = ''
     while (1){
          cout << "The old SEED finished, will you continue for SEED2 ? (Y/N)" << endl;
          cin >> continue_flag;
          if (continue_flag == 'N' || continue_flag == 'n'){
               exit();
          } else if (continue_flag == 'Y' || continue_flag == 'y'){
               break;
          }
     }

     statis_data << "Second Seed Result"<< endl;

     /*
     argv[2] = "input_from_output.fastq";    //输入每个cluster的core
     argv[4] = "output2.txt";  
     argv[6] = "3";   // SHIFT
     argv[8] = "5";   // MISMATCH*/

     strcpy(new_argv[tag_input],"input_from_output.fastq");
     strcpy(new_argv[tag_output],"output2.txt");
     itoa(argv_mismatch2, new_argv[new_argc-3], 10);
     itoa(argv_shift2, new_argv[new_argc-1], 10);

     
     // output.txt --> input2.fastq

     cout << endl << "Running Seed For the Second Time" << endl;
     old_seed(new_argc,new_argv);         //对第一遍seed后每个cluster的core进行seed

     /************************************************************************/
     /*  根据output2 再次合并                                                                    */
     /************************************************************************/
     time_t start,end;
     cout << endl << "Integrate the results of first seed and second seed" << endl;
     start = time(NULL);
     ifile.open("output2.txt",std::ifstream::in);
     ofile.close();
     ofile.open ("output3.txt", ios::out|ios::trunc);

     map<int, int>::iterator iter;

     int size_eachCluster = 0;
     int counter_output3 = 0;
     int cluster_id = 0;
     multimap<int,int> cluster_seed_size;
     if (ifile.is_open())
     {
          while ( getline (ifile,line) )
          {
               if (line[0]>= 'A' && line[0] <='T')
               {
                    if(counter_output3 != 0 && counter_output3 != 1){

                         ofile << "size " << size_eachCluster << endl;
                         cluster_seed_size.insert(pair<int,int>(size_eachCluster,cluster_id));
                         size_eachCluster = 0;
                    }

                    ofile << line << endl;
               }

               else if(line[0]>='0' && line[0] <= '9'){

                    int index_space = line.find_first_of('\t');
                    string left_str = line.substr(0,index_space);
                    string right_str = line.substr(index_space+1,line.length()-index_space);
                    int left = stoi(left_str);
                    int right = stoi(right_str);
                    cluster_seed_output2.insert(pair<int,int>(right,left));   // <new cluster id, old cluster id>

                    multimap<int, int>::iterator iter_lower,iter_upper,iter_output;
                    iter_lower = cluster_seed_output.lower_bound(right);
                    iter_upper = cluster_seed_output.upper_bound(right);

                    for(iter_output = iter_lower; iter_output != iter_upper; iter_output++){
                         cluster_seed_output3.insert(pair<int,int>(left,iter_output->second));
                         ofile << left << "\t" << iter_output->second << endl;
                         size_eachCluster++;
                    }

                    cluster_id = left;
               }
               counter_output3++;
          }
          ofile << "size " << size_eachCluster << endl;
          cluster_seed_size.insert(pair<int,int>(size_eachCluster,cluster_id));
          ifile.close();
     }else {
          cout << "Unable to open file"; 
     }
     //TEST
     //cluster_seed_output.clear();

     ofile.close();
     end = time(NULL);
     
     cout << endl << " " << end - start << "seconds" << endl << endl;



     cout << "MISMATCH " << argv[8] << " SHIFT " << argv[6] << endl;
     statis_data << "MISMATCH\t" << argv[8] << "\tSHIFT\t" << argv[6] << endl;

     double expectation=0;

     for (int i=1; i<=cluster_seed_size.size();i++)
     {
          int num = cluster_seed_size.count(i);
          
          if(num != 0){
               //cout << "size\t" << i << "\tnum\t" << cluster_seed_size.count(i) << "\tratio\t"<< (double) cluster_seed_size.count(i)/cluster_seed_size.size()<< endl;
               statis_data << "size\t" << i << "\tnum\t" << cluster_seed_size.count(i) << "\tratio\t"<< (double) cluster_seed_size.count(i)/cluster_seed_size.size()<< endl;
          }    
          // 排除 size 为1 的 cluster
          if(i!=0){
               expectation += (double) cluster_seed_size.count(i)/cluster_seed_size.size() * i;
          }

     }

     //expectation /= 1-(double) cluster_seed_size.count(1)/cluster_seed_size.size();

     int expectation_int = (int) (expectation+0.5);  //cluster 大小期望值


     cout << "expectation\t" << expectation_int << endl;

     /************************************************************************/
     /* k-means part    
        1. 
        2.*/
     /************************************************************************/

     // 1. prepare data
     

     ifile.open(input_firsttime_kmeans,std::ifstream::in);
     ofile.open ("output_after_twice_seed.txt", ios::out|ios::trunc);
     


     multimap<int,string> clid_seq;   // <cluster id, sequence>
     int max_clid_clid_seq = 0;
     int max_clid_clid_seq_snd = 1;
     multimap<int,string> sid_seq;   // read data from input.fastq into memory
     multimap<string,int> seq_sid;


     if (ifile.is_open())
     {
          int input_counter = 0;
          while(getline(ifile,line)){

               if (line.find("length") != line.npos && line[0] == '@')
               {
                    getline(ifile,line);
                    sid_seq.insert(pair<int,string>(input_counter,line));
                    seq_sid.insert(pair<string,int>(line,input_counter));
                    ofile_test << input_counter << "\t" << line << endl;
                    input_counter++;
               }
               
               
          }

          int couter_iter=0;
          for(iter = cluster_seed_output3.begin(); iter !=cluster_seed_output3.end(); iter++)
          {
               string seq = sid_seq.find(iter->second)->second;
               clid_seq.insert(pair<int,string>(iter->first,seq));
               couter_iter++;
               ofile << iter->first << " " << seq << endl;
               //cout << iter->first << ' ' << seq << endl;

               if(max_clid_clid_seq < iter->first){
                    max_clid_clid_seq = iter->first;
               }
          }
     

          ifile.close();
          cluster_seed_output3.clear();
     }else {
          cout << "Unable to open file"; 
     }    

     

     // use k-means to separate clusters of which the size is larger than expectation size
     multimap<int,string> clid_seq_kmeans_part;
     for (int i=expectation_int+1; i<=cluster_seed_size.size();i++)
     {
          int num = cluster_seed_size.count(i);
          if(num > 0){
               multimap<int, int>::iterator iter_lower,iter_upper,iter_output;
               iter_lower = cluster_seed_size.lower_bound(i);
               iter_upper = cluster_seed_size.upper_bound(i);

               for(iter_output = iter_lower; iter_output != iter_upper; iter_output++){
                    int cluster_id = iter_output->second;

                    pair <multimap<int,string>::iterator, multimap<int,string>::iterator> ret;
                    ret = clid_seq.equal_range(cluster_id);

                    for (multimap<int,string>::iterator it=ret.first; it!=ret.second; ++it){
                         clid_seq_kmeans_part.insert(pair<int,string>(it->first,it->second));
                    }

                    multimap<int,string> clid_seq_after_kmeans;
                    if (KM_flag == 1){
                         clid_seq_after_kmeans = K_Means_1(clid_seq_kmeans_part,expectation_int);
                    } else(KM_flag == 2){
                         clid_seq_after_kmeans = K_Means_2(clid_seq_kmeans_part,expectation_int);
                    }
                    
                    clid_seq.erase(cluster_id);

                    multimap<int, string>::iterator iter;
                    int clid;
                    for (iter = clid_seq_after_kmeans.begin(); iter != clid_seq_after_kmeans.end(); iter++)
                    {
                         if(iter == clid_seq_after_kmeans.begin()){
                              clid = iter->first;
                         }

                         int new_clid = iter->first;
                         if(new_clid != clid){
                              clid = new_clid;
                              max_clid_clid_seq_snd++;
                         }
                         clid_seq.insert(pair<int,string>(max_clid_clid_seq + max_clid_clid_seq_snd,iter->second));
                    }

                    max_clid_clid_seq = max_clid_clid_seq + max_clid_clid_seq_snd;
                    max_clid_clid_seq_snd = 1;

                    clid_seq_kmeans_part.clear();
                    clid_seq_after_kmeans.clear();
               }
          }    
     }

     statis_data << "after k-means" << endl;
     cluster_seed_size.clear();
     

     for (int i=0; i<=max_clid_clid_seq + max_clid_clid_seq_snd; i++)
     {
          int num = clid_seq.count(i);
          int size_per_cluster = 0;
          if (num > 0)
          {
               multimap<int, string>::iterator iter_lower,iter_upper,iter;
               iter_lower = clid_seq.lower_bound(i);
               iter_upper = clid_seq.upper_bound(i);
     
               for (iter = iter_lower; iter!=iter_upper; iter++)
               {
                    
                    size_per_cluster++;
               }

               cluster_seed_size.insert(pair<int,int>(size_per_cluster,iter_lower->first));
               size_per_cluster = 0;
          }
          
     }

     for (int i=1; i<cluster_seed_size.size(); i++)
     {
          int num = cluster_seed_size.count(i);

          if (num > 0)
          {
               statis_data << "size\t" << i << "\tnum\t" << num << endl;
          }

     }

     // 比对bowtie2

     // 原始seed j-index  GXX
     
     clid_seq.clear();
     
     multimap<int, int>::iterator iter_cluster_seed_output;

     for (iter_cluster_seed_output = cluster_seed_output.begin(); iter_cluster_seed_output!= cluster_seed_output.end(); iter_cluster_seed_output++)
     {
          int clid = iter_cluster_seed_output->first;
          int seqid = iter_cluster_seed_output->second;
          string seq = sid_seq.find(seqid)->second;
          clid_seq.insert(pair<int,string>(clid,seq));
     }
     
     

     ifstream iFile("bowtie2.txt");
     multimap<int,int> sid_offset;
     multimap<int,string> sid_ref;

     while (getline(iFile,line))
     {
          int space[3];
          int space_index=0;
          for (int i=0; i<line.size(); i++)
          {
               if (line[i] == ' ')
               {
                    space[space_index++] = i;
               }

               if (space_index == 3)
               {
                    break;
               }
          }

          int sid_int = stoi(line.substr(0,space[0])) -1 ;
          int offset = stoi(line.substr(space[2]+1,space[3]-space[2]));
          string ref = line.substr(space[1]+1,space[2]-space[1]);

          sid_offset.insert(pair<int,int>(sid_int,offset));
          sid_ref.insert(pair<int,string>(sid_int,ref));
     }

     multimap<int,int> clid_seqid;
     multimap<int, string>::iterator iter_clid_seq;
     ofile.close();
     ofile.open ("output_after_kmeans.txt", ios::out|ios::trunc);
     ofile << "CLID\t" << "SeqID" << endl;

     multimap<string, int>::iterator iter_seq_sid;
     for (iter_seq_sid = seq_sid.begin(); iter_seq_sid!= seq_sid.end(); iter_seq_sid++)
     {
          ofile_test << iter_seq_sid->first << "\t" << iter_seq_sid->second << endl;
     }

     //clid_seq  cluster_seed_output
     for (iter_clid_seq = clid_seq.begin(); iter_clid_seq!= clid_seq.end(); iter_clid_seq++)
     {
          int first = iter_clid_seq->first;
          string second = iter_clid_seq->second;
          int seqid = seq_sid.find(iter_clid_seq->second)->second;
          clid_seqid.insert(pair<int,int>(iter_clid_seq->first,seqid));


          ofile << iter_clid_seq->first << '\t' << seqid;

          if (sid_offset.find(seqid) != sid_offset.end())
          {
               multimap<int, int>::iterator iter_lower,iter_upper,iter_output;
               iter_lower = sid_offset.lower_bound(seqid);
               iter_upper = sid_offset.upper_bound(seqid);
               ofile << "\tbowtie2\t";
               for(iter_output = iter_lower; iter_output != iter_upper; iter_output++){
                    ofile <<iter_output->second << "\t";
               }
          }

          ofile << endl;

     }

     // jaccard index
     int union_jaccard = 0;
     int join_jaccard = 0;

     multimap<int, int>::iterator iter_clid_seqid;
     multimap<int,int> clid_sid_bowtie;

     for (iter_clid_seqid = clid_seqid.begin(); iter_clid_seqid!= clid_seqid.end(); iter_clid_seqid++)
     {
          if (sid_offset.find(iter_clid_seqid->second) != sid_offset.end())
          {
               if (clid_seqid.count(iter_clid_seqid->first) == 1)
               {
                    union_jaccard++;
                    join_jaccard++;
               } else if (clid_seqid.count(iter_clid_seqid->first) > 1)
               {
                    clid_sid_bowtie.insert(pair<int,int>(iter_clid_seqid->first,sid_offset.find(iter_clid_seqid->second)->first));


               }

          }
     }
     // 2004
     multimap<int, int>::iterator iter_clid_sid_bowtie;
     int min = 100;
     for (iter_clid_sid_bowtie = clid_sid_bowtie.begin(); iter_clid_sid_bowtie!= clid_sid_bowtie.end(); iter_clid_sid_bowtie++)
     {
          int tttt = iter_clid_sid_bowtie->first;
          int length = clid_sid_bowtie.count(iter_clid_sid_bowtie->first);
          if (length == 1)
          {
               continue;
          }

          multimap<int, int>::iterator iter_lower,iter_upper,iter_output;
          iter_lower = clid_sid_bowtie.lower_bound(iter_clid_sid_bowtie->first);
          iter_upper = clid_sid_bowtie.upper_bound(iter_clid_sid_bowtie->first);

          
          int* bowtie_pos = new int[length];
          string* bowtie_ref = new string[length];
          int index=0;
          for(iter_output = iter_lower; iter_output != iter_upper; iter_output++){        
               bowtie_pos[index] = sid_offset.find(iter_output->second)->second ;
               bowtie_ref[index] = sid_ref.find(iter_output->second)->second;
               
               index++;
               iter_clid_sid_bowtie++;
               union_jaccard++;
          }
          iter_clid_sid_bowtie--;


          for (int i=0; i<length ; i++)
          {
               for (int j=i+1; j<length; j++)
               {
                    int dis = abs(bowtie_pos[i] - bowtie_pos[j]);
                    int left = bowtie_pos[i];
                    int right = bowtie_pos[j];
                    string lefts = bowtie_ref[i];
                    string rights = bowtie_ref[j];
                    if (dis <= BOWTIE_DIFFERENCE && dis > 0 && bowtie_ref[i] == bowtie_ref[j])
                    {
                         join_jaccard++;  //15751
                         cout << left << " " << right << " " << dis << endl;
                    }
               }
          }
          delete[] bowtie_pos;
          delete[] bowtie_ref;
     }

     double jaccard_index = (double) join_jaccard / union_jaccard;
     cout << "jaccard index\t" << jaccard_index << "=" << join_jaccard << " / " << union_jaccard << endl;
     

     ofile.close();

     sid_seq.clear();
     seq_sid.clear();
     //K_Means(clid_seq);
     statis_data.close();
     
     return 1;
}
