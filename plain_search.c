#include <ncurses.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <limits.h>
#include <math.h>
#include <unistd.h>
#include "gnm.h"

extern int w;
extern int h;
extern int show_origin;            // 表示画面左上端のゲノム塩基位置番号 
extern int cpos;                   // カーソル塩基位置 
extern int cx;                     // カーソルX座標
extern int cy;                     // カーソルY座標
extern char message[512];          // 表示メッセージ 
extern int moveto;                 // 移動先塩基位置 
extern int mt;                     // mismatch tolerance 許容するミスマッチの数 
extern FILE *logfile;
extern int current_chr; 
extern int n_ref;                  // 入力ファスタに含まれる配列数 
extern chromosome_info *chrs;
extern int *included_chars;
extern int logout;

int plain_fsearch(char *sstr) // FWD
{
int i,j,k;
char rstr[512];
reverse_complement(sstr,rstr);
int cchr = current_chr;
cpos = show_origin + (cy-1) * w + cx;  // cpos : カーソル塩基位置
int flag = 0;
int mchr;
int mpos;
int mdir;
int len = strlen(sstr);
if(logout == 1)
 fprintf(logfile,"CCHR %d CPOS %d SSTR %s RSTR %s\n",cchr,cpos,sstr,rstr);
for(i=cchr;i<n_ref;i++)                // 現在のカーソル位置から配列の末尾まで
 {
 int from;
 if(i==cchr)
  from=cpos+1;
 else
  from=0;
 for(j=from;j<chrs[i].seq_len;j++)
  {
  int mmc = 0;
  int flag2= 0;
  for(k=0;k<len;k++)
   {
   if((sstr[k] != chrs[i].seq[j+k]) && (sstr[k] != chrs[i].seq[j+k]-32))  // 小文字を区別しない
    {
    mmc ++;
    if(mmc > mt)
     {
     flag2 = 1;
     break;
     }
    }
   }
  if(flag2 == 0) // 順方向マッチ
   {
   mchr = i; mpos = j; mdir = 1; flag = 1;
   break;
   }
  mmc = 0;
  flag2= 0;
  for(k=0;k<len;k++)
   {
   if((rstr[k] != chrs[i].seq[j+k]) && (rstr[k] != chrs[i].seq[j+k]-32))  // 小文字を区別しない
    {
    mmc ++;
    if(mmc > mt)
     {
     flag2 = 1;
     break;
     }
    }
   }
  if(flag2 == 0) // 逆方向マッチ
   {
   mchr = i; mpos = j; mdir = 2; flag = 1;
   break;
   }
  }
 if(flag == 1)
  break;
 }
if(flag == 0)
 {
 for(i=0;i<=cchr;i++)
  {
  for(j=0;j<chrs[i].seq_len;j++)
   {
   int mmc = 0;  
   int flag2= 0;
   for(k=0;k<len;k++)
    {     
    if((sstr[k] != chrs[i].seq[j+k]) && (sstr[k] != chrs[i].seq[j+k]-32))  // 小文字を区別しない 
     {     
     mmc ++;
     if(mmc > mt)
      {      
      flag2 = 1; 
      break;
      } 
     } 
    }     
   if(flag2 == 0) // 順方向マッチ 
    {     
    mchr = i; 
    mdir = 1; 
    mpos = j; 
    flag = 1; 
    break;
    }     
   mmc = 0;  
   flag2= 0;
   for(k=0;k<len;k++)
    {     
    if((rstr[k] != chrs[i].seq[j+k]) && (rstr[k] != chrs[i].seq[j+k]-32))  // 小文字を区別しない 
     {     
     mmc ++;
     if(mmc > mt)
      {      
      flag2 = 1; 
      break;
      } 
     } 
    }     
   if(flag2 == 0) // 逆方向マッチ 
    {     
    mchr = i; 
    mpos = j; 
    mdir = 2; 
    flag = 1; 
    break;
    } 
   }
  if(flag == 1)
   break;
  }
 }
if(flag == 1)
 {
 current_chr = mchr;
 moveto = mpos+1;
    
 if(mdir == 1)
  {
  for(i=0;i<len;i++)
   {
   if((sstr[i] == chrs[mchr].seq[mpos+i]) || (sstr[i] == chrs[mchr].seq[mpos+i]-32))
    {
    chrs[mchr].atr[mpos+i+1] = chrs[mchr].atr[mpos+i+1] | BIT_FLAG_1;
    }
   else
    {
    chrs[mchr].atr[mpos+i+1] = chrs[mchr].atr[mpos+i+1] | BIT_FLAG_2;
    }
   }
  }
 if(mdir == 2)
  {
  for(i=0;i<len;i++)
   {
   if((rstr[i] == chrs[mchr].seq[mpos+i]) || (rstr[i] == chrs[mchr].seq[mpos+i]-32))
    {
    chrs[mchr].atr[mpos+i+1] = chrs[mchr].atr[mpos+i+1] | BIT_FLAG_3;
    }
   else
    {
    chrs[mchr].atr[mpos+i+1] = chrs[mchr].atr[mpos+i+1] | BIT_FLAG_4;
    }
   }
  }
 move_to_pos(moveto);
    
 sprintf(message,"MATCH UPTO %d MisMatches for %s FOUND",mt,sstr);
 noecho();
 cbreak();
 return 0;
 }
return 0;
}

int plain_bsearch(char *sstr) // BWD
{
int i,j,k;
char rstr[512];
reverse_complement(sstr,rstr);
int cchr = current_chr;
cpos = show_origin + (cy-1) * w + cx - 2;        // cpos       : カーソル位置
int flag = 0;
int mchr;
int mpos;
int mdir;
int len = strlen(sstr);
if(logout == 1)
 fprintf(logfile,"CCHR %d CPOS %d SSTR %s RSTR %s\n",cchr,cpos,sstr,rstr);
for(i=cchr;i>=0;i--)
 {
 int from;
 if(i==cchr)
  from=cpos;
 else
  from=chrs[i].seq_len-1;
 for(j=from;j>=0;j--)
  {
  int mmc = 0;
  int flag2= 0;
  for(k=0;k<len;k++)
   {
   if((sstr[k] != chrs[i].seq[j+k]) && (sstr[k] != chrs[i].seq[j+k]-32))  // 小文字を区別しない
    {
    mmc ++;
    if(mmc > mt)
     {
     flag2 = 1;
     break;
     }
    }
   }
  if(flag2 == 0) // 順方向マッチ
   {
   mchr = i; mpos = j; mdir = 1; flag = 1;
   break;
   }
  mmc = 0;
  flag2= 0;
  for(k=0;k<len;k++)
   {
   if((rstr[k] != chrs[i].seq[j+k]) && (rstr[k] != chrs[i].seq[j+k]-32))  // 小文字を区別しない
    {
    mmc ++;
    if(mmc > mt)
     {
     flag2 = 1;
     break;
     }
    }
   }
  if(flag2 == 0) // 逆方向マッチ
   {
   mchr = i; mpos = j; mdir = 2; flag = 1;
   break;
   }
  }
 if(flag == 1)
  break;
 }
if(flag == 0)
 {
 for(i=n_ref-1;i>=cchr;i--)
  {
  for(j=chrs[i].seq_len-1;j>=0;j--)
   {
   int mmc = 0;  
   int flag2= 0;
   for(k=0;k<len;k++)
    {     
    if((sstr[k] != chrs[i].seq[j+k]) && (sstr[k] != chrs[i].seq[j+k]-32))  // 小文字を区別しない 
     {     
     mmc ++;
     if(mmc > mt)
      {      
      flag2 = 1; 
      break;
      } 
     } 
    }     
   if(flag2 == 0) // 順方向マッチ 
    {     
    mchr = i; 
    mdir = 1; 
    mpos = j; 
    flag = 1; 
    break;
    }     
   mmc = 0;  
   flag2= 0;
   for(k=0;k<len;k++)
    {     
    if((rstr[k] != chrs[i].seq[j+k]) && (rstr[k] != chrs[i].seq[j+k]-32))  // 小文字を区別しない 
     {     
     mmc ++;
     if(mmc > mt)
      {      
      flag2 = 1; 
      break;
      } 
     } 
    }     
   if(flag2 == 0) // 逆方向マッチ 
    {     
    mchr = i; 
    mpos = j; 
    mdir = 2; 
    flag = 1; 
    break;
    } 
   }
  if(flag == 1)
   break;
  }
 }
if(flag == 1)
 {
 current_chr = mchr;
 moveto = mpos+1;
    
 if(mdir == 1)
  {
  for(i=0;i<len;i++)
   {
   if((sstr[i] == chrs[mchr].seq[mpos+i]) || (sstr[i] == chrs[mchr].seq[mpos+i]-32))
    {
    chrs[mchr].atr[mpos+i+1] = chrs[mchr].atr[mpos+i+1] | BIT_FLAG_1;
    }
   else
    {
    chrs[mchr].atr[mpos+i+1] = chrs[mchr].atr[mpos+i+1] | BIT_FLAG_2;
    }
   }
  }
 if(mdir == 2)
  {
  for(i=0;i<len;i++)
   {
   if((rstr[i] == chrs[mchr].seq[mpos+i]) || (rstr[i] == chrs[mchr].seq[mpos+i]-32))
    {
    chrs[mchr].atr[mpos+i+1] = chrs[mchr].atr[mpos+i+1] | BIT_FLAG_3;
    }
   else
    {
    chrs[mchr].atr[mpos+i+1] = chrs[mchr].atr[mpos+i+1] | BIT_FLAG_4;
    }
   }
  }
 move_to_pos(moveto);
    
 sprintf(message,"MATCH UPTO %d Mismatches for %s FOUND",mt,sstr);
 noecho();
 cbreak();
 return 0;
 }
return 0;
}
