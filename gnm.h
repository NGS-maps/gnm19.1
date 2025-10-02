#include <ncurses.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <limits.h>
#include <math.h>
#include <unistd.h>

#define WORD_LEN          512  
#define BUFF_LEN          256  
#define MAX_THREAD_NUM     64   
#define SIGMA             256  

#define BIT_FLAG_1  0x01
#define BIT_FLAG_2  0x02
#define BIT_FLAG_3  0x04
#define BIT_FLAG_4  0x08
#define BIT_FLAG_5  0x10
#define BIT_FLAG_6  0x20
#define BIT_FLAG_7  0x40
#define BIT_FLAG_8  0x80


///////////////////////////////////////////////////////////////////////////////
////////////////////////////////STRUCTURE DEFINITION 構造体定義 
///////////////////////////////////////////////////////////////////////////////
typedef struct attr_info
{
unsigned char fwd_match    : 1;
unsigned char fwd_mismatch : 1;
unsigned char bwd_match    : 1;
unsigned char bwd_mismatch : 1;
unsigned char marked       : 1;
unsigned char visual       : 1;
}  attr_info;

typedef struct word_info             //   
 {
 char word[WORD_LEN];
 } word_info;

typedef struct hit_pos_info          //   
 {
 unsigned int chr; 
 unsigned int pos; 
 struct hit_pos_info *next;
 } hit_pos_info;

typedef struct thread_arg_t          ///////////////////////////////////  スレッド用構造体 
 {
 long int gfr; 
 long int gto; 
 unsigned int frch;
 unsigned int toch;
 unsigned int frpos;
 unsigned int topos;
 unsigned int thread_no;
 hit_pos_info *head;
 hit_pos_info *tail;
 int nh;
 char str[1024];
 } thread_arg_t;

typedef struct chromosome_info      //   
 {
 char *title;
 char *seq;
 unsigned char *atr;       // 描画上の各塩基位置でのアトリビュート 
 unsigned short int title_len;
 unsigned int seq_len;
 } chromosome_info;

typedef struct hpos_info            //   
 {
 unsigned int chr; 
 unsigned int pos; 
 unsigned int dir; // 1:fwd, 2:bwd
 } hpos_info;
///////////////////////////////////////////////////////////////////////////////
//////////////////////////////// STRUCTURE DEFINITION END  
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//////////////////////////////// FUNCTION PROTOTYPES       
///////////////////////////////////////////////////////////////////////////////
void readargs(int argc, char **argv);

int isnumber(char c);
int isincluded(char c);
int arenumbers(char *str);
int reverse_complement(char *instring,char *outstring);
int splitnum(char *str);
int chr_compare(int fc, int fp, int rc, int rp);

int plain_fsearch(char *sstr);
int plain_bsearch(char *sstr);
int thread_search(void *arg);
int seq_search(char *sstr,int forb);
int event_loop(void);
void handle_chr_pos(int *newchr,int *newpos,char *str);
unsigned int gch(long int gp);
unsigned int gpos(long int gp);
void atr_reset();

int move_to_pos(int pos);
