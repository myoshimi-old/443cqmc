#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include"443qmc.h"
#include"unit.h"
#include"hamming.h"

void qm(int num_true, UINT32* minterm_true,
        int num_dc,   UINT32* minterm_dc){
  int i, j, k, m;
  int num_minterm, max_value=0;
  UINT32* minterm;
  int blength=0;
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
  unit**  u;
  unit *uu, *uu2;
  int b;
  UINT8 tp;

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
          else printf(")\n");
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
        }
        lp = lp->next;
      }
    }
    // Primes
    printf("*** Primes ***\n");
    hammingView(primes, minterm);
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
      if(k < primearray[i]->arraysize){
        if(j == primearray[i]->array[k]){
          table[i][j] = 1;
          k++;
        }
      }
      else table[i][j] = 0;
    }
  }
  
  for(i=0;i<pnum;i++){
    printf("[%3d] ", i);
    for(j=0;j<num_true;j++)
      printf("%1d ", table[i][j]);
    printf("\n");
  }
  
  u = (unit**)malloc(sizeof(unit*)*num_true);
  b = (int)ceil((REAL)pnum / (sizeof(UINT8)*8));
  
  for(j=0;j<num_true;j++){
    num = 0;
    for(i=0;i<pnum;i++) num += table[i][j];
    u[j] = (unit*)malloc(sizeof(unit));
    u[j]->num = num;
    u[j]->array = (UINT8**)malloc(sizeof(UINT8*)*num);
    for(i=0;i<u[j]->num;i++){
      u[j]->array[i] = (UINT8*)malloc(sizeof(UINT8)*b);
      for(k=0;k<b;k++) u[j]->array[i][k] = 0;
    }
    
    k = 0;
    for(i=pnum-1;i>=0;i--){
      if(table[i][j] == 1){
        printf("m[%d] %d\n", j, i);
        unitSet(u[j], k, i);
        k++;
      }
    }
  }
  for(j=0;j<num_true;j++){
    printf("unum[%d] %d\n", j, u[j]->num);
  }
  // tableはもう使わないのでfreeする
  for(i=0;i<pnum;i++) free(table[i]);
  free(table);

  
  uu = unitClone(u[0], pnum);
  for(i=1;i<num_true;i++){
    uu2 = unitMerge(uu, u[i], pnum);
    for(j=0;j<uu2->num;j++){
      printf("(%d) ", j);
      for(k=b-1;k>=0;k--) printf("%1X", uu2->array[j][k]);
      printf("\n");
    }
    
    for(j=0;j<uu->num;j++) free(uu->array[j]);
    free(uu->array);
    free(uu);
    uu = uu2;
  }
  
  // uはもう使わないのでfreeする
  printf("Free u\n");
  for(i=0;i<num_true;i++){
    for(j=0;j<u[i]->num;j++) free(u[i]->array[j]);
    free(u[i]->array);
    free(u[i]);
  }
  free(u);
  
  for(i=0;i<uu->num;i++){
    m = 0;
    for(j=0;j<b;j++){
      tp = uu->array[i][j];
      
      for(k=0;k<8;k++){
        if(tp%2==1) printf("%d ", m);
        tp = tp / 2;
        m++;
      }
      printf("\n");
    }
  }


}






int main(int argc, char** argv){
  UINT32 *minterm_true, *minterm_dc;
  int num_true = 6;
  int num_dc   = 2;
  
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

  free(minterm_true);
  free(minterm_dc);

  return 0;
}


