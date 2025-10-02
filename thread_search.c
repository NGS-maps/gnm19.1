#include <ncurses.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <limits.h>
#include <math.h>
#include <unistd.h>
#include "gnm.h"

extern int  w, h;              // ウィンドウの縦h 横w サイズ 行数と列数
extern int show_origin;        // 表示画面左上端のゲノム塩基位置番号
extern int cpos;               // カーソル塩基位置
extern int cx;                 // カーソルX座標
extern int cy;                 // カーソルY座標
extern char message[512];      // 表示メッセージ
extern int moveto;             // 移動先塩基位置
extern int n_thread;
extern long int n_spot;
extern int mt;                 // mismatch tolerance 許容するミスマッチの数
extern FILE *logfile;
extern int included_chars[256];
extern hit_pos_info fhits;
extern hit_pos_info rhits;
extern int current_chr;
extern int n_ref;              // 入力ファスタに含まれる配列数（染色体数・スキャッフォルド数）
extern chromosome_info *chrs;
extern int n_hit_pos;
extern int n_fpos;
extern int n_rpos;
extern hpos_info *hpis;
extern int hit_point;
extern pthread_t *handle;                                    ///
extern thread_arg_t *targ;                                   ///
extern int logout;

///////////////////////////////////////////////////////////////////////////
int thread_search(void *arg)
{
int i,j,k;
thread_arg_t *targ = (thread_arg_t *)arg;
int tid;
int ni;
tid = targ->thread_no;
int len = strlen(targ->str);
hit_pos_info *tp;
hit_pos_info *p;
                                     // 各スレッドの守備範囲を探索
for(i=targ->toch;i>=targ->frch;i--)  // 染色体ループ リストへは前から加わっていくのでサフィックスは後ろから前に走る
 {
 if(i < 0)    // なぜかforループではiが負になるのを止められないので負になったらブレイク
  break;
 int pfr = 0;
 int pto = chrs[i].seq_len;
 if(i==targ->frch)              // スレッド最初の染色体
  pfr = targ->frpos;
 if(i==targ->toch)              // スレッド最後の染色体
  pto = targ->topos;
 if(pto == chrs[i].seq_len)     // 染色体終末端までptoが走るなら
  pto=pto-len+1;                // 探索文字数手前で止める
 for(j=pto-1;j>=pfr;j--)        // リストへは前から加わっていくのでサフィックスは後ろから前に走る
  {
  int flag = 0;
  int flag2 = 0;
  for(k=0;k<len;k++)
   {
   if((targ->str[k] != chrs[i].seq[j+k]) && (targ->str[k] != chrs[i].seq[j+k]-32))  // 小文字を区別しない
    {
    flag ++;
    if(flag > mt)                            // ミスマッチ許容 mt以下
     {
     flag2 = 1;
     break;
     }
    }
   }
  if(flag2 == 0)        // 一致があればヒットのリストに追加
   {
   tp = (hit_pos_info *)malloc(sizeof(hit_pos_info));
   tp->chr = i;
   tp->pos = j;
   p = targ->head;      // 今までのヘッドを退避
   targ->head = tp;     // 新しいヒットを頭に入れて
   tp -> next = p;      // 今までのヘッドを頭の次に
   if(targ->nh == 0)    // 最初のヒットを尻尾とする
    targ->tail =tp;
   targ->nh ++;
   }
  }
 }
return 1;
}

/////////////////////////////////////////////////////////////////////////////
int seq_search(char *sstr,int forb)   // forb 1:fwd -1:bwd
{
int i,j,k;
int hval;
int len;
hit_pos_info *p;
hit_pos_info *temp;
char rstr[512];
hit_pos_info *fp;
hit_pos_info *rp;
hit_pos_info *tp;
int cc = 0;

reverse_complement(sstr,rstr);

if(hpis != NULL)
 {
 free(hpis);
 hpis = NULL;
 }
len = strlen(sstr);

/////////////////////////////////////////////////////
 
for(i=0;i<n_thread;i++)                              // 順鎖サーチ
 {
 targ[i].thread_no = i;
 targ[i].nh = 0;
 targ[i].head = NULL;
 targ[i].tail = NULL;
 strcpy(targ[i].str,sstr);
 pthread_create(&handle[i],NULL,(void *)thread_search,(void *)&targ[i]);  // スレッド分岐
 }

for(i=0;i<n_thread;i++)
 {
 pthread_join(handle[i],NULL);                       // スレッド合流
 }

///////////////////////////////////////////////////// 各スレッドのヒットを一本のスレッドにまとめる

n_fpos = 0;
fhits.next =  NULL;
for(i=n_thread-1;i>=0;i--)
 {
 if(targ[i].head != NULL)     // そのスレッドにヒットがあるなら
  {
  p = fhits.next;
  fhits.next = targ[i].head;
  targ[i].tail -> next = p;
  n_fpos += targ[i].nh;  // 前向き総ヒット数にそのスレッドのヒット数を加える
  }
 }

/////////////////////////////////////////////////////

for(i=0;i<n_thread;i++)                             // 逆鎖サーチ
 {
 targ[i].thread_no = i;
 targ[i].nh = 0;
 targ[i].head = NULL;
 targ[i].tail = NULL;
 strcpy(targ[i].str,rstr);
 pthread_create(&handle[i],NULL,(void *)thread_search,(void *)&targ[i]); // スレッド分岐
 }

for(i=0;i<n_thread;i++)
 {
 pthread_join(handle[i],NULL);                      // スレッド合流
 }

///////////////////////////////////////////////////// 各スレッドのヒットを一本のリストにまとめる

n_rpos = 0;
rhits.next =  NULL;
p = &rhits;
for(i=n_thread-1;i>=0;i--)
 {
 if(targ[i].head != NULL)                           // そのスレッドにヒットがあるなら
  {
  p = rhits.next;
  rhits.next = targ[i].head;
  targ[i].tail -> next = p;
  n_rpos += targ[i].nh;                             // 後向き総ヒット数にそのスレッドのヒット数を加える
  }
 }

/////////////////////////////////////////////////////

 int judge=0;                                       // リストから構造体の配列にデータを移行
 fp = fhits.next;
 rp = rhits.next;
 n_hit_pos = n_fpos + n_rpos;
 hpis = (hpos_info *)malloc(sizeof(hpos_info) * (n_fpos + n_rpos + 2));
 while(1)
  {
  if(fp == NULL)
   if(rp == NULL)
    break;
   else
    judge = 2;
  else
   if(rp == NULL)
    judge = 1;
   else
    judge = chr_compare(fp->chr,fp->pos,rp->chr,rp->pos);   // fp が前か同じなら1 後なら2 
  if(judge == 1)
   {
   hpis[cc].chr = fp->chr;
   hpis[cc].pos = fp->pos;
   hpis[cc].dir = 1;                              // 順方向ヒット
   tp = fp->next;
   free(fp);
   fp = tp;
   cc++;
   }
  else
   {
   hpis[cc].chr = rp->chr;
   hpis[cc].pos = rp->pos;
   hpis[cc].dir = 2;                              // 相補逆鎖
   tp = rp->next;
   free(rp);
   rp = tp;
   cc++;
   }
  }

/////////////////////////////////////////////////////

for(i=0;i<n_hit_pos;i++)
 {
 if(hpis[i].dir == 1)                             // 順方向
  {
  for(j=0;j<len;j++)
   {
   if((chrs[hpis[i].chr].seq[hpis[i].pos+j] == sstr[j]) || (chrs[hpis[i].chr].seq[hpis[i].pos+j] == sstr[j] + 32))
    {
    chrs[hpis[i].chr].atr[hpis[i].pos+j+1] = chrs[hpis[i].chr].atr[hpis[i].pos+j+1] | BIT_FLAG_1;
    }
   else
    {
    chrs[hpis[i].chr].atr[hpis[i].pos+j+1] = chrs[hpis[i].chr].atr[hpis[i].pos+j+1] | BIT_FLAG_2;
    }
   }
  }
 else                                             // 相補逆鎖
  {
  for(j=0;j<len;j++)
   {
   if((chrs[hpis[i].chr].seq[hpis[i].pos+j] == rstr[j]) || (chrs[hpis[i].chr].seq[hpis[i].pos+j] == rstr[j] + 32))
    {
    chrs[hpis[i].chr].atr[hpis[i].pos+j+1] = chrs[hpis[i].chr].atr[hpis[i].pos+j+1] | BIT_FLAG_3;
    }
   else
    {
    chrs[hpis[i].chr].atr[hpis[i].pos+j+1] = chrs[hpis[i].chr].atr[hpis[i].pos+j+1] | BIT_FLAG_4;
    }
   }
  }
 }

if(cc != 0)                                       // 検索結果のヒットが０でなければ
 {
 int fflag = 0;
 for(i=0;i<cc-1;i++)
  {
  int fchr = hpis[i].chr;
  int fpos = hpis[i].pos;
  int tchr = hpis[i+1].chr;
  int tpos = hpis[i+1].pos;
  int cchr = current_chr;
  cpos = show_origin + (cy-1) * w + cx;           // cpos       : カーソル位置
  if(cchr > fchr)                                 // 現在の染色体番号がiの染色体番号より大きい
   {
   if(cchr < tchr)                                // 現在の染色体番号がiの染色体番号より小さい
    {
    if(forb == 1)
     hit_point = i+1;
    else
     hit_point = i;
    fflag = 1;
    break;
    }  
   else 
    {
    if(cchr == tchr)
     {
     if(cpos <= tpos)
      {
      if(forb == 1)
       hit_point = i+1;
      else
       hit_point = i;
      fflag = 1;
      break;
      }
     }
    }
   }
  else
   {
   if(cchr == fchr)
    {
    if(cpos >= fpos)
     {
     if(cchr < tchr)
      {
      if(forb == 1)
       hit_point = i+1;
      else
       hit_point = i;
      fflag = 1;
      break;
      }
     else
      {
      if(cpos <= tpos)
       {
       if(forb == 1)
        hit_point = i+1;
       else
        hit_point = i;
       fflag = 1;
       break;
       }
      }
     }
    }
   }
  }
 if(fflag == 0)
  {
  hit_point = 0;
  }

 current_chr = hpis[hit_point].chr;
 moveto = hpis[hit_point].pos+1;
 move_to_pos(moveto);
 }

sprintf(message,"MATCH UPTO %d Mismatches for %s is %10d FWD:%8d REV:%8d",mt,sstr,n_fpos+n_rpos,n_fpos,n_rpos);
noecho();
cbreak();
return 0;
}
