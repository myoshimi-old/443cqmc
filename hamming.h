#ifndef _INC_HAMMING
#define _INC_HAMMING

#include"443qmc.h"

typedef struct _hammingData{
  int     flag;
  UINT32  mask;
  int arraysize;
  UINT32* array;
  struct _hammingData* next;
} hammingData;

typedef struct _hamming{
  int len;
  struct _hammingData *head, *tail;
} hamming;

hamming* hammingNew();
void hammingView(hamming* h, UINT32* term);
void hammingAppend(hamming* hin, int flag, UINT32 mask);
void hammingMatch(hamming* h, hammingData *a, hammingData *b, UINT32 *term);
int bitcount(UINT32 n);
int hammingDistance(UINT32 a, UINT32 b, UINT32 ma, UINT32 mb);
void hammingFree(hamming* h);

#endif
