#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"hamming.h"
#include"unit.h"

int unitBitcount(UINT8* d, int cnum){
  int i;
  int b, s, n;
  
  b = (int)ceil((REAL)cnum / (sizeof(UINT8)*8));

  s = 0;
  for(i=b-1;i>=0;i--){
    n = d[i];
    n = (UINT8)( ((n & 0xAA) >> 1) + (n & 0x55) );
    n = (UINT8)( ((n & 0xCC) >> 2) + (n & 0x33) );
    n = (UINT8)( ((n & 0xF0) >> 4) + (n & 0x0F) );
    s += n;
  }
  return s;
}

void unitSet(unit* u, int aidx, int d){
  int a, m;
  a = d / (sizeof(UINT8)*8);
  m = d % (sizeof(UINT8)*8);
  u->array[aidx][a] |= (1 << m);
}

UINT32 udata2uint(UINT8* a, int n){
  int i;
  int b;
  UINT32 ret;
  
  b = (int)ceil((REAL)n / (sizeof(UINT8)*8));
  ret = 0;
  
  for(i=b-1;i>=0;i--){
    ret *= 256;
    ret += a[i];
  }
  
  return ret;
}

unit* unitClone(unit* uin, int cnum){
  int i, j, b;
  unit* ret;
  
  b = (int)ceil((REAL)cnum / (sizeof(UINT8)*8));
  
  ret = (unit*)malloc(sizeof(unit));
  ret->num = uin->num;
  ret->array = (UINT8**)malloc(sizeof(UINT8*)*ret->num);

  for(i=0; i<ret->num;i++){
    ret->array[i] = (UINT8*)malloc(sizeof(UINT8)*b);
    for(j=b-1;j>=0;j--) ret->array[i][j] = uin->array[i][j];
  }
  return ret;
}

void unitFree(unit* u){
  int i;
  for(i=0;i<u->num;i++) free(u->array[i]);
  free(u->array);
  free(u);
}

unit* unitMerge(unit* u1, unit* u2, int n){
  unit* ret;
  int i, j, k;
  int b;
  int num, nn;
  UINT8** array;

  b = (int)ceil((REAL)n / (sizeof(UINT8)*8));
  
  //ret = (unit*)malloc(sizeof(unit));
  //ret->num = u1->num * u2->num;
  //ret->array = (UINT8**)malloc(sizeof(UINT8*)*ret->num);

  num = u1->num * u2->num;
  array = (UINT8**)malloc(sizeof(UINT8*)*num);
  
  for(i=0;i<num;i++){
    array[i] = (UINT8*)malloc(sizeof(UINT8)*b);
    for(j=b-1;j>=0;j--) array[i][j] = 0;
  }
  for(i=0;i<u1->num;i++){
    for(j=0;j<u2->num;j++){
      for(k=b-1;k>=0;k--){
        array[i*u2->num+j][k] = u1->array[i][k] | u2->array[j][k];
      }
    }
  }
  
  // ソート
  printf("S:%d b:%d ", num, b);
  uintSort(array, num, b);

  // 重複の削除＋返り値への代入
  nn = 1;
  for(i=0;i<num-1;i++){
    if(uintCmp(array[i], array[i+1], b) != EQ) nn++;
  }
  
  printf(" n:%d\n", nn);

  ret = (unit*)malloc(sizeof(unit));
  ret->num = nn;
  ret->array = (UINT8**)malloc(sizeof(UINT8*)*nn);
  
  k = 0;
  for(i=0;i<nn;i++){
    ret->array[i] = (UINT8*)malloc(sizeof(UINT8)*b);
  }

  k=0;
  uintCpy(ret->array[k], array[0], b);
  k++;
  for(i=1;i<num;i++){
    if(uintCmp(array[i-1], array[i], b) != EQ){
      uintCpy(ret->array[k], array[i], b);
      k++;
    }
  }
  
  return ret;
}

int uintCmp(UINT8* a1, UINT8* a2, int b){
  int i;
  int t=0;
  int ret=0;
  
  t=0; i=b;
  while((t == 0) & (i>=0)){
    t = (int)(a1[i] - a2[i]);
    i--;
  };
  /*
  printf("t %d\n", t);
  for(i=b-1;i>=0;i--){
    printf("%X", a1[i]);
  }
  printf("\n");
  for(i=b-1;i>=0;i--){
    printf("%X", a2[i]);
  }
  printf("\n");
  */
  if(t<0) ret = LT;
  else if(t>0) ret = GT;
  else ret = EQ;
  
  return ret;
}

UINT8* uintMed3(UINT8* x, UINT8* y, UINT8* z, int b){
  UINT8* ret;
  if(uintCmp(x, y, b) == LT){
    // x < y
    if(uintCmp(y, z, b) == LT)      ret = y;
    else if(uintCmp(z, x, b) == LT) ret = x;
    else                            ret = z;
  }
  else{
    // y < x
    if(uintCmp(z, y, b) == LT)      ret = z;
    else if(uintCmp(x, z, b) == LT) ret = x;
    else                            ret = y;
  }
  return ret;
}

void uintQuicksort(UINT8** a, int left, int right, int b){
  int i, j;
  UINT8 *tmp, *pivot;
  if(left < right){
    i = left; j = right;
    pivot = uintMed3(a[i], a[i+(j-i)/2], a[j], b);
    while(1){
      while(uintCmp(a[i], pivot, b) == LT) i++;
      while(uintCmp(pivot, a[j], b) == LT) j--;
      if(i>=j) break;
      tmp = a[i];
      a[i] = a[j];
      a[j] = tmp;
      i++; j--;
    }

    //if(right - left > 100000)
    //printf("sort %d %d\n", left, right);
    
    uintQuicksort(a, left,  i-1, b);
    uintQuicksort(a, j+1, right, b);
  }
}

void uintSort(UINT8** array, int num, int b){
  /*
    int i, j;
  UINT8* t;
  for(i=0;i<num;i++){
    for(j=i+1;j<num;j++){
      if(uintCmp(array[i], array[j], b) == GT){
        t = array[i];
        array[i] = array[j];
        array[j] = t;
      }
    }
  }
  */

  uintQuicksort(array, 0, num-1, b);
}

void uintCpy(UINT8* d, UINT8* s, int b){
  int i;
  for(i=b-1;i>=0;i--) d[i] = s[i];
}

/*
void unitFree(unit* u){
  int i;
  for(i=0;i<u->num;i++){
    free(u->array[i]);
  }
  free(u->array);
  free(u);
}
*/

/*
unit* unitMerge(unit* u1, unit* u2, int n){
  unit* ret;
  int i, j, k;
  int b;

  b = (int)ceil((REAL)n / (sizeof(UINT8)*8));
  
  ret = (unit*)malloc(sizeof(unit));
  ret->num = u1->num * u2->num;
  ret->array = (UINT8**)malloc(sizeof(UINT8*)*ret->num);
  for(i=0;i<ret->num;i++){
    ret->array[i] = (UINT8*)malloc(sizeof(UINT8)*b);
    for(j=b-1;j>=0;j--) ret->array[i][j] = 0;
  }
  for(i=0;i<u1->num;i++){
    for(j=0;j<u2->num;j++){
      for(k=b-1;k>=0;k--){
        ret->array[i*u2->num+j][k] = u1->array[i][k] | u2->array[j][k];
      }
    }
  }
  
  return ret;
}
*/

/*
unit* unitOR(unit* u1, unit* u2, int n){
  int b, i;
  unit* r;
  
  b = (int)ceil((REAL)n / (sizeof(UINT8)*8));
  r = (unit*)malloc(sizeof(unit));
  r->array = (UINT8*)malloc(sizeof(UINT8)*b);
  
  for(i=0;i<b;i++){
    r->array[i] = u1->array[i] | u2->array[i];
  }

  return r;
}

// n : 要素数
unit* unitMalloc(int n){
  unit* u;
  int b;
  b = (int)ceil((REAL)n / (sizeof(UINT8)*8));
  
  u = (unit*)malloc(sizeof(unit));
  u->array = (UINT8*)malloc(sizeof(UINT8)*b);
  
  return u;
}

// n : 番号
void unitAdd(unit* u, int n){
  int idx;
  idx = n / (sizeof(UINT8)*8);
  u->array[idx] |= (1 << (n % (sizeof(UINT8)*8)));
}

void unitView(unit* u, int n){
  int b, i;
  b = (int)ceil((REAL)n / (sizeof(UINT8)*8));
  for(i=b-1;i>=0;i--){
    printf("%01X", u->array[i]);
  }
  printf("\n");
}


int unitBitcount(unit* u, int n){
  int b, num, i;
  b = (int)ceil((REAL)n / (sizeof(UINT8)*8));
  num = 0;
  for(i = 0;i<b;i++){
    num+=bitcount(u->array[i]);
  }
  return num;
}
*/


