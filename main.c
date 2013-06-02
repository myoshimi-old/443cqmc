#include<stdio.h>
#include<unistd.h>
#include<stdlib.h>
#include<math.h>

#include"443qmc.h"
#include"unit.h"
#include"unit32.h"
#include"hamming.h"

#include <sys/time.h>

double gettimeofday_sec(){
  struct timeval t;
  gettimeofday(&t, NULL);
  return (double)t.tv_sec + (double)t.tv_usec * 1e-6;
}


void qm(int num_true, UINT32* minterm_true,
        int num_dc,   UINT32* minterm_dc){
  int i, j, k, m;
  int num_minterm, max_value=0;
  UINT32* minterm;
  int blength=0;
  int un;
  unit32 **ua;
  int idx;
  int bc;
  int cycle;
  hamming **hp, **hp2;
  hammingData *lp, *rp;
  hamming *primes;
  int pnum;
  hammingData** primearray;
  int num;
  //UINT8** unit;
  UINT8** table;
  UINT8* a;
  UINT8* p;
  
  unit**  u;
  unit *uu, *uu2;
  int b;
  double t1, t2;

  int* cnum;
  int maxidx, maxnum;

  printf("======================== Computing =========================\n");
  printf("# of 1  : %d\n", num_true);
  printf("# of DC : %d\n", num_dc);
  num_minterm = num_true + num_dc;
  
  minterm = (UINT32*)malloc(sizeof(UINT32)*num_minterm);
  for(i=0;i<num_true;i++) minterm[i] = minterm_true[i];
  for(i=0;i<num_dc;i++)   minterm[num_true+i] = minterm_dc[i];
  // trueおよびdcの最大値(入力値の最大値)を求める
  for(i=0;i<num_minterm;i++)
    if(minterm[i] > max_value) max_value = minterm[i];
  while(max_value >= (int)pow(2.0,(REAL)blength)) blength++;
  printf("Max Value: %d (%d [bits])\n", max_value, blength);
  
  printf("======================== Input Data ========================\n");
  for(i=0;i<num_minterm;i++) printf("[%5d]: %8d\n", i, minterm[i]);
  printf("------------------------------------------------------------\n");
  // ハミング重みごとにまとめる
  hp = (hamming**)malloc(sizeof(hamming*)*(blength+1));
  for(i=0;i<blength+1;i++){
    hp[i] = hammingNew();
  }
  for(i=0;i<num_minterm;i++){
    bc = bitcount(minterm[i]);
    hammingAppend(hp[bc], 0, 0);
    hp[bc]->head->array = (UINT32*)malloc(sizeof(UINT32));
    hp[bc]->head->arraysize = 1;
    hp[bc]->head->array[0] = i;
  }
  
  // 主項リストの初期化
  primes = hammingNew();
  
  // 状態を表示
  for(i=0;i<blength+1;i++){
    printf("*** %d ***\n", i);
    hammingView(hp[i], minterm);
  }

  cycle = 0;
  do{
    printf("--------------------- Cycles %4d --------------------------\n", cycle);
    hp2 = (hamming**)malloc(sizeof(hamming*)*(blength+1));
    for(i=0;i<blength+1;i++) hp2[i] = hammingNew();
    for(i=0;i<blength;i++){
      printf("[%d] [%d]\n", i, i+1);
      lp = hp[i]->head;
      while(lp != NULL){
        rp = hp[i+1]->head;
        while(rp != NULL){
          hammingMatch(hp2[i], lp, rp, minterm);
          rp=rp->next;
        }
        lp=lp->next;
      }
    }

    // hpのフラグを調査して主項を抽出
    for(i=0;i<blength+1;i++){
      printf("Flag **[%d]**\n", i);
      lp = hp[i]->head;
      while(lp !=NULL){
        for(j=0;j<lp->arraysize;j++){
          if(j==0) printf(" %d(", lp->flag);
          printf("%d", lp->array[j]);
          if(j != lp->arraysize-1) printf(",");
          else printf(")");
        }
        
        if(lp->flag == 0){
          rp = (hammingData*)malloc(sizeof(hammingData));
          rp->flag = 0; rp->mask=lp->mask;
          rp->arraysize = lp->arraysize;
          rp->array = (UINT32*)malloc(sizeof(UINT32)*lp->arraysize);
          for(j=0;j<lp->arraysize;j++) rp->array[j] = lp->array[j];
          rp->next = primes->head;
          primes->head = rp;
          primes->len++;
          printf(" Add primes");
        }
        printf("\n");
        lp = lp->next;
      }
    }
    // Primes
    printf("*** Primes ***\n");
    //hammingView(primes, minterm);
    lp = primes->head;
    while(lp != NULL){
      bitExpression(minterm, lp, blength);
      printf("\n");
      lp = lp->next;
    }
    printf("------------------------------------------------------------\n");

    // hpを開放する
    for(i=0;i<blength+1;i++) hammingFree(hp[i]);
    free(hp);
    hp = hp2;
    
    // 状態を表示
    for(i=0;i<blength+1;i++){
      printf("*** %d ***\n", i);
      hammingView(hp[i], minterm);
    }
    
    // 残りの式の数を調査
    num = 0;
    for(i=0;i<blength+1;i++){
      num += hp[i]->len;
    }
    printf("Remain: %d\n", num);
    printf("------------------------------------------------------------\n");
    cycle++;
  }while(num != 0);

  for(i=0;i<blength+1;i++) hammingFree(hp[i]);
  free(hp);

  // Primes
  printf("*** Primes ***\n");
  hammingView(primes, minterm);
  printf("------------------------------------------------------------\n");
  // 扱いやすいようにコピー
  pnum = primes->len;
  primearray = (hammingData**)malloc(sizeof(hammingData*)*pnum);
  lp = primes->head;
  for(i=0;i<pnum;i++){
    primearray[i] = lp;
    lp = lp->next;
  }
  free(primes);
  
  // Tables
  // 主項の数だけ要素を持つ配列を，最小項の数だけ作る
  table = (UINT8**)malloc(sizeof(UINT8*)*pnum);
  for(i=0;i<pnum;i++){
    table[i] = (UINT8*)malloc(sizeof(UINT8)*num_true);
    k=0;
    for(j=0;j<num_true;j++){
      table[i][j] = 0;
      if(k < primearray[i]->arraysize){
        if(j == primearray[i]->array[k]){
          table[i][j] = 1;
          k++;
        }
      }
    }
  }

  p = (UINT8*)malloc(sizeof(UINT8)*pnum);
  for(i=0;i<pnum;i++) p[i] = 0;
  
  // 必須項を見つける
  a = (UINT8*)malloc(sizeof(UINT8)*num_true);
  for(j=0;j<num_true;j++) a[j] = 0;
  
  cnum = (int*)malloc(sizeof(int)*pnum);
  
  for(j=0;j<num_true;j++){
    num = 0;
    for(i=0;i<pnum;i++){
      if(num == 0 && table[i][j]==1) m = i;
      num += table[i][j];
    }
    if(num == 1){
      for(k=0;k<num_true;k++) a[k] |= table[m][k];
      p[m] = 1;
    }
    
    printf("m[%d] %d", j, num);
    if(num == 1){
      printf(" * e[%d]", m);
    }
    printf("\n");
  }

  m = 0;
  for(j=0;j<num_true;j++)
    if(a[j] == 0) m++;
  printf("Remains Zero:%d\n", m);

  while(m != 0){
    for(i=0;i<pnum;i++) cnum[i] = 0;
    for(i=0;i<pnum;i++){
      if(p[i] == 0){
        for(j=0;j<num_true;j++)
          if(a[j] == 0 && table[i][j] == 1) cnum[i]++;
      }
    }
    
    maxidx = 0; maxnum = 0;
    for(i=0;i<pnum;i++){
      if(maxnum < cnum[i]){
        maxidx = i;
        maxnum = cnum[i];
      }
    }
    
    printf("maximum coverage e[%d]: %d bits\n", maxidx, maxnum);
    for(j=0;j<num_true;j++) a[j] |= table[maxidx][j];
    p[maxidx] = 1;

    m = 0;
    for(j=0;j<num_true;j++) if(a[j] == 0) m++;
    printf("Remains Zero:%d\n", m);
  }

  m = 0;
  for(i=0;i<pnum;i++) if(p[i] != 0) m++;
  printf("%d terms:\n", m);

  for(i=0;i<pnum;i++){
    if(p[i] != 0){
      bitExpression(minterm, primearray[i], blength);
      printf("\n");
    }
  }


  /*
  #ifdef DEBUG
  for(i=0;i<uu->num;i++){
    printf("[%d](%d) ", i, ua[i]->num);
    for(j=0;j<ua[i]->num;j++) printf("%u ", ua[i]->array[j]);
    printf("\n");
  }
  #endif

  printf("Bits: %d\n", uu->num);
  
  {
    int min, min_num;
    min = 0;
    min_num = ua[min]->num;
    for(i=0;i<uu->num;i++){
      if(ua[i]->num < min_num){
        min = i;
        min_num = ua[i]->num;
      }
    }
    printf("[%d](%d) ", min, min_num);
    for(j=0;j<ua[min]->num;j++){
      bitExpression(minterm, primearray[ua[min]->array[j]], blength);
      printf(" ");
    }
    printf("\n");
  }
  */
}

int main(int argc, char** argv){
  UINT32 *minterm_true, *minterm_dc;
  int num_true = 6;
  int num_dc   = 2;
  int i;

  FILE *fp_true, *fp_dc, *fp;
  char line[256];
  /*
  minterm_true = (UINT32*)malloc(sizeof(UINT32)*num_true);
  minterm_dc   = (UINT32*)malloc(sizeof(UINT32)*num_dc);
  
  minterm_true[ 0] = 4;
  minterm_true[ 1] = 12;
  minterm_true[ 2] = 15;
  minterm_true[ 3] = 10;
  minterm_true[ 4] = 11;
  minterm_true[ 5] = 8;
  
  minterm_dc[0] = 14;
  minterm_dc[1] = 9;
  
  qm(num_true, minterm_true, num_dc, minterm_dc);
  */

  if(argc != 3){
    printf("%s [file_true] [file_dc]\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  if((fp_true = fopen(argv[1], "r")) == NULL){
    printf("file open error.\n");
    exit(EXIT_FAILURE);
  }
  if((fp_dc = fopen(argv[2], "r")) == NULL){
    printf("file open error.\n");
    exit(EXIT_FAILURE);
  }
  fp = fp_true;
  num_true = 0;
  while(fgets(line,sizeof(line),fp) != NULL) num_true++;
  fp = fp_dc;
  num_dc = 0;
  while(fgets(line,sizeof(line),fp) != NULL) num_dc++;
  fclose(fp_true);
  fclose(fp_dc);


  minterm_true = (UINT32*)malloc(sizeof(UINT32)*num_true);
  minterm_dc   = (UINT32*)malloc(sizeof(UINT32)*num_dc);
  
  if((fp_true = fopen(argv[1], "r")) == NULL){
    printf("file open error.\n");
    exit(EXIT_FAILURE);
  }
  if((fp_dc = fopen(argv[2], "r")) == NULL){
    printf("file open error.\n");
    exit(EXIT_FAILURE);
  }
  fp=fp_true; i=0;
  while(fgets(line,sizeof(line),fp) != NULL){
    minterm_true[i] = atoi(line);
    i++;
  }
  fp=fp_dc; i=0;
  while(fgets(line,sizeof(line),fp) != NULL){
    minterm_dc[i] = atoi(line);
    i++;
  }
  fclose(fp_true);
  fclose(fp_dc);

  qm(num_true, minterm_true, num_dc, minterm_dc);

  printf("Free\n");
  free(minterm_true);
  free(minterm_dc);
  
  return 0;
}


