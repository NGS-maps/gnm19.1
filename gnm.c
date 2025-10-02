#include <ncurses.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <limits.h>
#include <math.h>
#include <unistd.h>
#include "gnm.h"

///////////////////////////////////////////////////////////////////////////////
///////////////////// Grobal Variables ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
int  w, h;                  // ウィンドウの縦h横wサイズ 行数と列数
unsigned int  genomelen=1;  // 読み込fasta配列長
int show_origin = 1;        // 表示画面左上端のゲノム塩基位置番号
int cpos;                   // カーソル塩基位置
int cx=0;                   // カーソルX座標
int cy=1;                   // カーソルY座標
char message[512];          // 表示メッセージ
int moveto=1;               // 移動先塩基位置
int cs_len;                 // コマンド入力文字列の長さ
int n_thread = 0;
long int n_spot = 0;
int mt = 0;                 // mismatch tolerance 許容するミスマッチの数
char com_str[512];          // コマンド入力文字列  シーケンス
char search_str[512];       // コマンド入力文字列  シーケンス
char arg1[256];
char arg2[256];
FILE *logfile;
int search_mode = 0;        // 0: befoer first search 1: plain_fsearch, 2:plain_bsearch, 
                            // 3: seq_search(fwd) 4: seq_search(bwd)
int mark_chr[256];
int mark_pos[256];
int included_chars[256];
int logout = 0;

int vmode = 0;
int rmode = 0;
int vichr;
int vipos;
int vepos;
int vic[256];
int vip[256];
int vep[256];

int  read_gff = 0;
char gff_flnm[256];
char       tstr[256];
hit_pos_info fhits;
hit_pos_info rhits;
int current_chr = 0;
int n_ref;                  // 入力ファスタに含まれる配列数
int word_len = WORD_LEN;
chromosome_info *chrs;
int n_hit_pos = 0;
int n_fpos = 0;
int n_rpos = 0;
hpos_info *hpis=NULL;
int hit_point = 0;
long int     ref_len = 0;
int          *ref_lens;
unsigned int *ref_offset;
word_info    *titles;
char         *refseq;
pthread_t *handle;
thread_arg_t *targ;

int max_so;              // ゲノム表示位置の最大値（ウィンドウサイズに依存
int endpos;
int cc=0;
char buff[1024];          // 読み込みバッファ
char in_flnm[128];       // 入力ファイル名

int pos;                 // 現在位置
char  show_char;         // 表示するキャラクタ

int  key;                // キーボード入力一時保管char
char k2;
int fflag;
///////////////////////////////////////////////////////////////////////////////
///////////////////// Functions ///////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

int isincluded(char c)
 {
 if((int)included_chars[c] == 0)
  {
  return 1;
  }   
 return 0;
 }
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
/// MAIN FUNCTION                                            
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
int i,j,k,l;             // ループカウンタ
FILE *genomefile;        // 入力ファイルハンドル
FILE *bf2infile;         // 入力bf2ファイルハンドル
char c;                  // 読み込みキャラクタバッファ
char title[1024];        // 入力fastaタイトル行
int sc=0;
int len;
int bf2in = 0;

if(logout == 1)
 logfile = fopen("log.log","w");

readargs(argc,argv);

sprintf(message,"                             \n");

word_len=WORD_LEN;

n_thread = sysconf(_SC_NPROCESSORS_CONF);
printf("Number of threads on this computer is %d\n",n_thread);

for(i=0;i<256;i++)
 {
 mark_chr[i] = -1;
 mark_pos[i] = -1;
 vic[i]      = -1;
 vip[i]      = -1;
 vep[i]      = -1;
 }

strcpy(in_flnm,arg1);

//////////////////////////////////////////////////////////////　FASTA読み込み
 genomefile = fopen(in_flnm,"r");
 if(genomefile == NULL)
  {
  printf("Failed to open input file %s\n",in_flnm);
  exit(1);
  }
 
 //////////////////////////////////////////////////// FIRST SCAN
 n_ref = 0;
 while(fgets(buff,1024,genomefile))
  {
  if(buff[0] == '>')
   n_ref ++;
  }
 printf("Number of sequences(chromosomes) in %s is %d\n",in_flnm,n_ref);
 chrs = (chromosome_info *)malloc(sizeof(chromosome_info)*(n_ref+2));
 if(chrs == NULL)
  {
  printf("Failed to allocate memory for chromosome informations\n");
  exit(1);
  }
 for(i=0;i<n_ref;i++)
  {
  chrs[i].seq_len = 0;
  chrs[i].title_len = 0;
  }
 rewind(genomefile);

 ref_lens = (int *)malloc(sizeof(int)*(n_ref+2));                            //
 ref_offset = (unsigned int *)malloc(sizeof(unsigned int)*(n_ref+2));        //
 titles = (word_info *)malloc(sizeof(word_info)*n_ref);                      //

 //////////////////////////////////////////////////// SECOND SCAN
 cc = 0;
 while(fgets(buff,1024,genomefile))
  {
  if(buff[0] == '>')
   {
   len = strlen(buff);
   buff[len-1] = '\0';
   chrs[cc].title_len = len;
   chrs[cc].title = (char *)malloc(sizeof(char)*(len+5));
   if(chrs[cc].title == NULL)
    {
    printf("Failed to allocate memory\n");
    exit(1);
    }
   strcpy(chrs[cc].title,buff);
   cc ++;
   }
  else
   {
   chrs[cc-1].seq_len += (strlen(buff)-1);
   }
  }
 ref_offset[0] = 0;                                  //
 for(i=0;i<n_ref;i++)                                //
  {                                                  //
  ref_lens[i] = chrs[i].seq_len;                     //
  if(i != 0)                                         //
   ref_offset[i]  = ref_len;                         //
  ref_len     += chrs[i].seq_len;                    //
  }                                                  //
 ref_offset[n_ref] = ref_len;                        //
 ////////////////////////////////////////////////////
 for(i=0;i<n_ref;i++)
  {
  chrs[i].seq = (char *)malloc(sizeof(char) * (chrs[i].seq_len + 5));
  if(chrs[i].seq == NULL)
   {
   printf("Failed to allocate memory\n");
   exit(1);
   }
  chrs[i].atr = (unsigned char *)malloc(sizeof(unsigned char) * (chrs[i].seq_len + 5));
  if(chrs[i].atr == NULL)
   {
   printf("Failed to allocate memory\n");
   exit(1);
   }
  }
 rewind(genomefile);
 
 //////////////////////////////////////////////////// THIRD SCAN
 cc = 0;
 sc = 0;
 while(fgets(buff,1024,genomefile))
  {
  if(buff[0] == '>')
   {
   cc ++;
   sc = 0;
   }
  else
   {
   len = strlen(buff)-1;
   buff[len] = '\0';
   strcpy(&chrs[cc-1].seq[sc],buff);
   sc += len;
   }
  }

for(i=0;i<256;i++)
 included_chars[i] = 0;
for(i=0;i<n_ref;i++)
 {
 for(j=0;j<chrs[i].seq_len;j++)
  {
  included_chars[chrs[i].seq[j]] = 1;
  }
 }
if(logout == 1)
 {
 for(i=0;i<256;i++)
  {
  if(included_chars[i] == 1)
   {
   fprintf(logfile,"%c\n",(char)i);
   }
  }
 }

////////////////////////////////////////////////////////////// fasta/bf2 読み込み終了
//////////////////////////////////////////////////////////////////////////////////////
////  文字列探索用のスレッド関数の準備
//////////////////////////////////////////////////////////////////////////////////////
handle = (pthread_t *)malloc(sizeof(pthread_t) * n_thread);
targ   = (thread_arg_t *)malloc(sizeof(thread_arg_t) * n_thread);
long int asegment = ref_len / n_thread;               // スレッドあたりの受け持ち塩基長
//////////////////////////////////////////////////////// スレッドあたりの守備範囲の設定
targ[0].gfr = 0;
targ[0].gto = asegment;
for(i=1;i<n_thread;i++)
 {
 targ[i].gfr = targ[i-1].gto + 1;
 targ[i].gto = targ[i-1].gto + asegment;
 }
targ[0].frch  = 0;
targ[0].frpos = 0;
targ[0].toch  = gch(targ[0].gto);
targ[0].topos = gpos(targ[0].gto);
targ[0].head = NULL;
targ[0].tail = NULL;
for(i=1;i<n_thread-1;i++)
 {
 targ[i].frch  = gch(targ[i].gfr);
 targ[i].frpos = gpos(targ[i].gfr);
 targ[i].toch  = gch(targ[i].gto);
 targ[i].topos = gpos(targ[i].gto);
 targ[i].head = NULL;
 targ[i].tail = NULL;
 }
targ[n_thread-1].frch  = gch(targ[n_thread-1].gfr);
targ[n_thread-1].frpos = gpos(targ[n_thread-1].gfr);
targ[n_thread-1].toch  = n_ref-1;
targ[n_thread-1].topos = chrs[n_ref-1].seq_len;
targ[n_thread-1].head = NULL;
targ[n_thread-1].tail = NULL;
//////////////////////////////////////////////////////// スレッドあたりの守備範囲の設定

//////////////////////////////////////
event_loop();
//////////////////////////////////////

endwin();
return (0);
}
