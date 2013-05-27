#include<stdio.h>
#include<stdlib.h>
#include<limits.h>
#include"hamming.h"

void hammingView(hamming* h, UINT32* term){
  hammingData *lp;
  int i;
  
  printf(" n: %d \n", h->len);
  lp = h->head;
  while(lp != NULL){
    for(i=0;i<lp->arraysize;i++){
      printf("[%d] %04X ", lp->array[i], term[lp->array[i]]);
    }
    printf("\n");
    lp = lp->next;
  }
}

hamming* hammingNew(){
  hamming* h;
  h = (hamming*)malloc(sizeof(hamming));
  h->len = 0;
  h->head = NULL;
  h->tail = NULL;
  return h;
}

void hammingFree(hamming* h){
  hammingData *p, *d;
  p = h->head;
  while(p != NULL){
    d = p;
    p = p->next;
    
    free(d->array);
    free(d);
  }
  free(h);
}

int bitcount(UINT32 n){
    n -= n >> 1 & 0x55555555;
    n = (n & 0x33333333) + (n >> 2 & 0x33333333);
    n = (n + (n >> 4)) & 0x0f0f0f0f;
    return (n * 0x01010101) >> 24;
}

int hammingDistance(UINT32 a, UINT32 b, UINT32 ma, UINT32 mb){
  UINT32 d, dm;
  d  = (a|ma|mb) ^ (b|ma|mb);
  dm = ma ^ mb;
  return bitcount(d)+bitcount(dm);
}

void hammingAppend(hamming* hin, int flag, UINT32 mask){
  hammingData *h;

  h = (hammingData*)malloc(sizeof(hammingData));
  h->flag = flag;
  h->mask = mask;
  h->next = NULL;
  
  if(hin->len == 0) hin->tail = h;
  else              h->next = hin->head;
  hin->len++;
  hin->head = h;
}

void bitExpression(UINT32* term, hammingData* hd, int blen){
  int i;
  UINT32 d;
  UINT32 m;
  char c[8];
  
  d = term[hd->array[0]];
  m = hd->mask;

  for(i=0;i<blen;i++){
    if(m % 2) c[i] = 'X';
    else      c[i] = '0'+d%2;
    m /= 2;
    d /= 2;
  }
  
  for(i=blen-1;i>=0;i--){
    printf("%c", c[i]);
  }


}

void hammingMatch(hamming* h, hammingData *a, hammingData *b, UINT32 *term){
  int i;
  int midx, nidx;
  UINT32 m, n;
  int dist;
  UINT32 dc, cb;
  hammingData *hd;
  hammingData *p;
  UINT32 diff;
  int flag;

  dist = hammingDistance(term[a->array[0]], term[b->array[0]],
                         a->mask, b->mask);
  printf("dist: %d  a:[%4d]%04X b:[%4d]%04X\n",
         dist, a->array[0], term[a->array[0]],
         b->array[0], term[b->array[0]]);

  if(dist == 1){
    dc = ((term[a->array[0]]|a->mask) ^ (term[b->array[0]]|b->mask))
      | a->mask | b->mask;
    cb = bitcount(term[a->array[0]]|dc)-bitcount(dc);
    printf("  DC:%04X, cb:%d\n",dc, cb);
    // 既にハミング配列に存在していたら処理しない
    // 存在チェック
    
    //hammingAppend(h, 0, dc);
    hd = (hammingData*)malloc(sizeof(hammingData));
    //hd = h->head;
    hd->flag = 0; hd->mask = dc; hd->next=NULL;
    hd->arraysize = a->arraysize + b->arraysize;
    hd->array = (UINT32*)malloc(sizeof(UINT32)*hd->arraysize);
    printf("  ArraySize : %d\n", hd->arraysize);
    // マージソート
    midx=0; nidx=0;
    for(i=0;i<hd->arraysize;i++){
      if(midx < a->arraysize) m = a->array[midx];
      else m = UINT_MAX;
      if(nidx < b->arraysize) n = b->array[nidx];
      else n = UINT_MAX;
      if(m < n){ hd->array[i] = m; midx++; }
      else{      hd->array[i] = n; nidx++; }
    }
    //for(i=0;i<a->arraysize;i++) hd->array[i] = a->array[i];
    //for(i=0;i<b->arraysize;i++) hd->array[a->arraysize+i] = b->array[i];
    flag=0; p=h->head;
    while((p != NULL) && (flag == 0)){
      if((p->mask == hd->mask) && (hd->arraysize == p->arraysize)){
        diff = 0;
        for(i=0;i<hd->arraysize;i++){
          diff += p->array[i] - hd->array[i];
          if(diff != 0) break;
        }
        if(diff == 0) flag = 1;
      }
      p=p->next;
    }
    // 同じものは発見されなかった
    printf("  Flag: %d\n", flag);
    if(flag == 0){
      if(h->len == 0) h->tail  = hd;
      else            hd->next = h->head;
      h->head = hd;
      h->len++;
    }
    a->flag++; b->flag++;
  }
}

