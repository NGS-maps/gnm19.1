#include <ncurses.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <limits.h>
#include <math.h>
#include <unistd.h>
#include "gnm.h"

extern int n_ref;
extern int ref_len;
extern int *ref_lens;
extern int genomelen;
extern unsigned int *ref_offset;
extern chromosome_info *chrs;
extern FILE *logfile;
extern int  w, h;              // ウィンドウの縦h横wサイズ 行数と列数
extern int show_origin;        // 表示画面左上端のゲノム塩基位置番号 
extern int cpos;               // カーソル塩基位置 
extern int cx;                 // カーソルX座標
extern int cy;                 // カーソルY座標
extern int logout;

extern int n_thread;
extern pthread_t *handle;
extern thread_arg_t *targ;

///////////////////////////////////////////////////////////////////////////////
int isnumber(char c)
 {
 if((c >= '0') && (c <= '9'))
  {
  return 1;
  }
 return 0;
 }
///////////////////////////////////////////////////////////////////////////////
int arenumbers(char *str)
 {
 int i,len;
 len = strlen(str);
 int flag = 0; 
 for(i=0;i<len;i++)
  {
  if(!((str[i] >= '0') && (str[i] <= '9')))
   {
   flag = 1; 
   break;
   }
  }
 if(flag == 0)
  return(atoi(str));
 else 
  return (-99999);
 }
///////////////////////////////////////////////////////////////////////////////
int chr_compare(int fc, int fp, int rc, int rp)
 {
 if(fc < rc)       // fc < rc
  return 1;
 else
  {
  if(fc > rc)      // fc > rc
   {
   return 2;
   }
  else             // fc = rc
   {
   if(fp <= rp)
    return 1;
   else
    return 2;
   }
  }
 }
///////////////////////////////////////////////////////////////////////////////
int reverse_complement(char *instring,char *outstring) // 相補逆鎖を返す関数
{
int i,j,k; 
int len;
char c;
len = strlen(instring);
for(i=0;i<len;i++)
 {
 switch (instring[len-i-1])    
  {
  case 'A': outstring[i] = 'T';  break;
  case 'T': outstring[i] = 'A';  break;
  case 'G': outstring[i] = 'C';  break;
  case 'C': outstring[i] = 'G';  break;
  case 'a': outstring[i] = 't';  break;
  case 't': outstring[i] = 'a';  break;
  case 'g': outstring[i] = 'c';  break;
  case 'c': outstring[i] = 'g';  break;
  case 'n': outstring[i] = 'N';  break;
  case '.': outstring[i] = 'N';  break;
  case '-': outstring[i] = 'N';  break;
  default : outstring[i] = 'N';  break;
  }
 }
outstring[len] = '\0';
return 0;
}
///////////////////////////////////////////////////////////////////////////////
int splitnum(char *str)
 {
 int len = strlen(str);
 int num = 0;
 int flag = 0;
 
 if(isnumber(str[1]))
  {
  return 999;               // エラーコード　数値は1桁9まで
  }

 if(isnumber(str[0]))       // １文字目が数字なら返す数値（許容するミスマッチ数）を読んで文字列を修正
  {
  num=(int)str[0] - 48;
  for(int i=0;i<len-1;i++)
   str[i] = str[i+1];
  flag = 0;
  for(int i=0;i<len;i++)
   {
   if(isincluded(str[i]))
    {
    flag = 1;
    break;
    }
   }
  if(flag == 1)
   return 9999;               // エラーコード　探索ワードはアルファベットのみ
  else
   return num;
  }
 else                       // １文字目が数字でなければ　０と元の文字列をそのまま返す
  {
  flag = 0;
  for(int i=0;i<len;i++)
   {
   if(isincluded(str[i]))
    {
    flag = 1;
    break;
    }
   }
  if(flag == 1)
   return 9999;               // エラーコード　探索ワードはアルファベットのみ
  else
   return 0;
  }
 }
///////////////////////////////////////////////////////////////////////////////
void handle_chr_pos(int *newchr,int *newpos,char *str)
 {
 int i;
 char chrstr[512];
 char posstr[512];
 char *p; 
 if(strstr(str,":") != NULL)       // :を挟んで、染色体:塩基位置
  {
  strcpy(chrstr,str);
  p=strstr(chrstr,":");
  *p = '\0';
  strcpy(posstr,strstr(str,":")+1);
  *newpos = arenumbers(posstr);
  }
 else                              // :なしで、染色体のみ
  {
  strcpy(chrstr,str);
  *newpos = 1;
  }
 *newchr = arenumbers(chrstr);
 if(*newchr == -99999)             // 染色体を表すストリングが数値でない
  {
  for(i=0;i<n_ref;i++)
   {   
   if(strstr(chrs[i].title,chrstr) != NULL)
    {   
    break;
    }   
   }   
  *newchr = i+1;
  }
 if(logout==1)
  fprintf(logfile,"### %s %s %s %5d %10d\n",str,chrstr,posstr,*newchr,*newpos);
 }
///////////////////////////////////////////////////////////////////////////////
unsigned int gch(long int gp) 
 {
 unsigned int i;
 for(i = 0;i< n_ref;i++)
  {
  if((gp > ref_offset[i]) && (gp <= ref_offset[i+1]))
   break;
  }
 return i;
 }
///////////////////////////////////////////////////////////////////////////////
unsigned int gpos(long int gp) 
 {
 unsigned int i;
 unsigned int num;
 for(i = 0;i< n_ref;i++)
  {
  if((gp > ref_offset[i]) && (gp <= ref_offset[i+1]))
   {   
   num = gp - ref_offset[i];
   break;
   }   
  }
 return num;
 }
///////////////////////////////////////////////////////////////////////////
void atr_reset()
 {
 int i,j;
 for(i=0;i<n_ref;i++)
  {
  memset(chrs[i].atr,0,chrs[i].seq_len);
  }
 }
///////////////////////////////////////////////////////////////////////////
int move_to_pos(int pos)                // 指定塩基位置への移動（座標計算）show_originとcx,cyを決めて返す
{
//fprintf(logfile,"MOVETO %d\n",pos);
if(pos <= 0)
 pos = 1;
if(pos > genomelen)
 pos = genomelen;
int bottom_margin = h/2;
if((pos >= show_origin) && (pos < show_origin +((h-4)*w))) // pos へ移動 飛び先が表示画面内 
 {    
 cy = ((pos-show_origin) / w)+1;
 if(pos%w == 0)
  cx = w-1; 
 else 
  cx = (pos%w)-1;
 }    
else                                    // 飛び先が表示画面外
 {    
 if(pos > cpos)                         // 現在位置より下流（後ろ）
  {    
  show_origin=((pos/w)*w)+1 - (h-bottom_margin-2)*w;
  cy=h-bottom_margin-1;
  }    
 else                                   // 現在位置より上流（前）
  {    
  if(pos%w == 0)
   show_origin = ((pos/w) * w)+1-(bottom_margin-1)*w;
  else 
   show_origin = ((pos/w) * w)+1-(bottom_margin-2)*w;
  cy=1+bottom_margin-2;
  if(pos%w == 0)
   show_origin = ((pos/w) * w)+1-w;
  else 
   show_origin = ((pos/w) * w)+1;
  cy=1;
  }    
 if(pos%w == 0)
  cx = w-1; 
 else 
  cx = (pos%w)-1;
 }    
return 1;
}
/////////////////////////////////////////////////////////////////////////////
