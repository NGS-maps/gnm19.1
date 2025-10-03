#include <ncurses.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <limits.h>
#include <math.h>
#include <unistd.h>
#include "gnm.h"

#define ctrl(x) (x & 0x1F)

//////// ncurses.h で 0:COLOR_BLACK,1:COLOR_RED,2:COLOR_GREEN,3:COLOR_YELLOW,4:COLOR_BLUE
////////              5:COLOR_MAGENTA,6:COLOR_CYAN,7:COLOR_WHITEまで定義されている  
#define COLOR_LIGHT_RED    8    
#define COLOR_LIGHT_GREEN  9
#define COLOR_LIGHT_BLUE  10  
#define COLOR_GRAY        11   
#define COLOR_LIGHT_GRAY  12  
#define COLOR_LIGHT_GRAY2 13
#define COLOR_DARK_GRAY   14   
#define COLOR_LIGHT_PINK  15  

///////////////////////////////////////////////////////////////////////////////
///////////////////// Grobal Variables ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
extern int  w, h;                  // ウィンドウの縦h横wサイズ 行数と列数
extern unsigned int  genomelen;    // 読み込fasta配列長
extern int show_origin;            // 表示画面左上端のゲノム塩基位置番号
extern int cpos;                   // カーソル塩基位置
extern int cx;                     // カーソルX座標
extern int cy;                     // カーソルY座標
extern char message[712];          // 表示メッセージ
extern int moveto;                 // 移動先塩基位置
extern char com_str[512];          // コマンド入力文字列  シーケンス
extern char search_str[512];       // コマンド入力文字列  シーケンス
extern FILE *logfile;
extern int search_mode;            // 0: before the first search 1: plain_fsearch(fwd), 2:plain_search(bwd), 
                                   // 3: seq_search(fwd) 4: seq_search(bwd)
extern int logout;
extern int mark_chr[256];
extern int mark_pos[256];
extern int vmode;
extern int vichr;
extern int vipos;
extern int vepos;
extern int vic[256];
extern int vip[256];
extern int vep[256];
extern int current_chr;
extern int n_ref;                  // 入力ファスタに含まれる配列数
extern chromosome_info *chrs;
extern int n_hit_pos;
extern hpos_info *hpis;
extern int hit_point;
extern int max_so;                 // ゲノム表示位置の最大値（ウィンドウサイズに依存
extern int endpos;
extern int cc;
extern char buff[1024];            // 読み込みバッファ
extern char in_flnm[128];          // 入力ファイル名
extern int pos;                    // 現在位置
extern char show_char;             // 表示するキャラクタ
extern int  key;                   // キーボード入力一時保管char
extern char k2; 
extern int fflag;
extern int mt;
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
int event_loop(void)
{
int i,j,k,l,m,n;

int newchr;
int newpos;
//////////////////////////////////////////////////////////////
initscr();                                    ///////////スクリーン初期化
start_color();
init_color(COLOR_WHITE,1000,1000,1000);
init_color(COLOR_LIGHT_RED,1000,900,900);
init_color(COLOR_LIGHT_GREEN,900,1000,900);
init_color(COLOR_LIGHT_BLUE,800,1000,1000);
init_color(COLOR_LIGHT_PINK,1000,900,900);

init_color(COLOR_GRAY,500,500,500);
init_color(COLOR_LIGHT_GRAY,900,900,900);
init_color(COLOR_LIGHT_GRAY2,800,800,800);
init_color(COLOR_DARK_GRAY,300,300,300);

init_pair(0,COLOR_BLACK,COLOR_WHITE);         // Basic Black on White


init_pair(1,COLOR_RED,COLOR_WHITE);
init_pair(2,COLOR_GREEN,COLOR_WHITE);
init_pair(3,COLOR_YELLOW,COLOR_WHITE);
init_pair(4,COLOR_BLUE,COLOR_WHITE);
init_pair(5,COLOR_MAGENTA,COLOR_WHITE);
init_pair(6,COLOR_CYAN,COLOR_WHITE);
init_pair(7,COLOR_WHITE,COLOR_WHITE);
init_pair(10,COLOR_BLACK,COLOR_RED);
init_pair(11,COLOR_RED,COLOR_RED);
init_pair(12,COLOR_GREEN,COLOR_RED);
init_pair(13,COLOR_YELLOW,COLOR_RED);
init_pair(14,COLOR_BLUE,COLOR_RED);
init_pair(15,COLOR_MAGENTA,COLOR_RED);
init_pair(16,COLOR_CYAN,COLOR_RED);
init_pair(17,COLOR_WHITE,COLOR_RED);
init_pair(20,COLOR_BLACK,COLOR_GREEN);
init_pair(21,COLOR_RED,COLOR_GREEN);          //  2
init_pair(22,COLOR_GREEN,COLOR_GREEN);
init_pair(23,COLOR_YELLOW,COLOR_GREEN);
init_pair(24,COLOR_BLUE,COLOR_GREEN);
init_pair(25,COLOR_MAGENTA,COLOR_GREEN);
init_pair(26,COLOR_CYAN,COLOR_GREEN);
init_pair(27,COLOR_WHITE,COLOR_GREEN);
init_pair(30,COLOR_BLACK,COLOR_YELLOW);
init_pair(31,COLOR_RED,COLOR_YELLOW);
init_pair(32,COLOR_GREEN,COLOR_YELLOW);
init_pair(33,COLOR_YELLOW,COLOR_YELLOW);
init_pair(34,COLOR_BLUE,COLOR_YELLOW);
init_pair(35,COLOR_MAGENTA,COLOR_YELLOW);
init_pair(36,COLOR_CYAN,COLOR_YELLOW);
init_pair(37,COLOR_WHITE,COLOR_YELLOW);
init_pair(40,COLOR_BLACK,COLOR_BLUE);
init_pair(41,COLOR_RED,COLOR_BLUE);
init_pair(42,COLOR_GREEN,COLOR_BLUE);
init_pair(43,COLOR_YELLOW,COLOR_BLUE);
init_pair(44,COLOR_BLUE,COLOR_BLUE);
init_pair(45,COLOR_MAGENTA,COLOR_BLUE);
init_pair(46,COLOR_CYAN,COLOR_BLUE);
init_pair(47,COLOR_WHITE,COLOR_BLUE);
init_pair(50,COLOR_BLACK,COLOR_MAGENTA);
init_pair(51,COLOR_RED,COLOR_MAGENTA);
init_pair(52,COLOR_GREEN,COLOR_MAGENTA);
init_pair(53,COLOR_YELLOW,COLOR_MAGENTA);
init_pair(54,COLOR_BLUE,COLOR_MAGENTA);
init_pair(55,COLOR_MAGENTA,COLOR_MAGENTA);
init_pair(56,COLOR_CYAN,COLOR_MAGENTA);
init_pair(57,COLOR_WHITE,COLOR_MAGENTA);      //  3 bwd_mismatch
init_pair(60,COLOR_BLACK,COLOR_CYAN);
init_pair(61,COLOR_RED,COLOR_LIGHT_BLUE);     //  1 fwd_match
init_pair(62,COLOR_GREEN,COLOR_CYAN);
init_pair(63,COLOR_YELLOW,COLOR_CYAN);
init_pair(64,COLOR_BLUE,COLOR_CYAN);
init_pair(65,COLOR_MAGENTA,COLOR_CYAN);
init_pair(66,COLOR_CYAN,COLOR_CYAN);
init_pair(67,COLOR_WHITE,COLOR_CYAN);
init_pair(100,COLOR_BLACK,COLOR_WHITE);
init_pair(101,COLOR_RED,COLOR_LIGHT_PINK);    //  2 bwd_match
init_pair(102,COLOR_BLACK,COLOR_LIGHT_GRAY);  // Basic Black on Light Gray
init_pair(103,COLOR_WHITE,COLOR_GRAY);        // Basic Black on Dark Gray
init_pair(104,COLOR_BLACK,COLOR_GRAY);        // Basic Black on Gray
init_pair(105,COLOR_BLACK,COLOR_DARK_GRAY);   // Basic Black on Dark Gray
init_pair(106,COLOR_BLACK,COLOR_LIGHT_GRAY2); // Basic Black on Light Gray2 (darker)
init_pair(107,COLOR_WHITE,COLOR_BLUE);        //  4 fwd_mismatch

bkgd(COLOR_PAIR(100));

noecho();
cbreak();
keypad(stdscr, TRUE);

//////////////////////////////////////////////// MAIN LOOP
//////////////////////////////////////////////// MAIN LOOP
while (1) 
//////////////////////////////////////////////// MAIN LOOP
//////////////////////////////////////////////// MAIN LOOP
 {                                             ///////////イベントループ
 getmaxyx(stdscr, h, w);                       // ウィンドウサイズ取得
 genomelen = chrs[current_chr].seq_len;        // 表示中の染色体塩基数
 max_so=genomelen-((genomelen%w)+w*(h-5));     // 表示位置の終点　（ゲノム長の終わりがある場合）
 erase();                                      // 画面消去

///////////////////////////////////////////////////////////////////画面表示
                                               // show_origin: 表示ページの左上角の塩基位置
 cpos = show_origin + (cy-1) * w + cx;         // cpos       : カーソル塩基位置  cx:カーソルの画面上x座標 
                                               //                                cy:カーソルの画面上y座標
 endpos = show_origin+w*(h-4)-1;               // 表示範囲の末尾の塩基位置
 if(endpos > genomelen)                        // endposがgenomelenより大きければgenomelenに揃える
  endpos=genomelen;
 attrset(A_REVERSE);                           // 表示反転 
 //////////////////////////////////////////////////////////// 染色体情報表示
 move(0,0);                                    // 天辺行
 cc = w - 35;                                  // 右端から左に35文字戻る
 if(vmode == 0)
  sprintf(buff,"  Chr: %8d/%8d             %*s ",current_chr+1,n_ref,cc,chrs[current_chr].title);
 else
  sprintf(buff,"V Chr: %8d/%8d             %*s ",current_chr+1,n_ref,cc,chrs[current_chr].title);
 buff[w] = '\0';
 printw("%s",buff);
 //////////////////////////////////////////////////////////// ゲノム位置情報表示
 move(h-3,0);                                  // 最下行の三つ上
 cc = w - 73;                                  // 右端から左に73文字戻る
 sprintf(buff," Cursor at:%9d  Showing:%10d~%10d  of: %10d     %*s ",
              cpos,show_origin,endpos,genomelen,cc,in_flnm);
 buff[w] = '\0';
 printw("%s",buff);

 attrset(A_NORMAL);                            // それ以外は通常表示: 反転解除
 //////////////////////////////////////////////////////////// メッセージ表示
 move(h-1,0);                                  // 最下行
 sprintf(buff," %s",message);
 buff[w] = '\0';
 printw("%s",buff);
////////////////////////////////////////////////////////////// 画面表示 DRAW
 pos = show_origin;
 for(i=1;i<h-3;i++)                            // 描画域  行数   （縦）    
  {
  for(j=0;j<w;j++)                             // 描画域　列数   （横）
   {
   attrset(A_NORMAL);
   if(pos <= genomelen)
    {
    if((i==cy) && (j==cx))                     // カーソル位置なら
     {
     attrset(A_REVERSE);                       // 表示反転 
     }
    else
     {
     attrset(A_NORMAL);                        // それ以外は通常表示
     if((vmode == 1) && (vichr == current_chr))// ビジュアルモード 選択中
      {
      if(vipos < cpos)
       {
       if((vipos <= pos) && (pos <= cpos))
        {
        attrset(A_BOLD);
        attrset(COLOR_PAIR(102));
        }
       }
      else
       {
       if((vipos >= pos) && (pos >= cpos))
        {
        attrset(A_BOLD);
        attrset(COLOR_PAIR(102));
        }
       }
      }                                        // ここまでビジュアルモード

     if((chrs[current_chr].atr[pos] & BIT_FLAG_1) != 0)  // fwd_match
      attrset(COLOR_PAIR(61));
     if((chrs[current_chr].atr[pos] & BIT_FLAG_2) != 0)  // fwd_mismatch
      attrset(COLOR_PAIR(107));
     if((chrs[current_chr].atr[pos] & BIT_FLAG_3) != 0)  // bwd_match
      attrset(COLOR_PAIR(101));
     if((chrs[current_chr].atr[pos] & BIT_FLAG_4) != 0)  // bwd_mismatch
      attrset(COLOR_PAIR(57));
     if((chrs[current_chr].atr[pos] & BIT_FLAG_5) != 0)  // visual
      attrset(COLOR_PAIR(102));
     if((chrs[current_chr].atr[pos] & BIT_FLAG_6) != 0)  // marked
      attrset(COLOR_PAIR(103));
     }
    }
   if(pos > genomelen)                         // ゲノム範囲外（末端以降）には ~ を表示
    show_char = '~';
   else
    if(pos < 1)
     show_char = '~';
    else
     show_char = chrs[current_chr].seq[pos-1]; // 各ゲノム位置の塩基を表示
   mvaddch(i,j,show_char);
   pos++;
   }
  }
 refresh();
 ///////////////////////////////////////////////////////////////////画面表示  DRAWここまで
 
 ////////////////////////////////////////////////////////////////////////キー入力へのレスポンス
 key = getch();                               // キーボード入力の文字をkeyへ代入

// if (key == 'q')                              // q ならプログラム終了
//  break;

 if((key >= '0') && (key <= '9'))             // キー入力が数値: 続くキー入力の繰り返し数（複数桁あり） 
                                              // i.e. Nh で上流へN塩基　移動
  { //////////////////////////////////////////// INPUT STARTS WITH NUMBERS ////////////////////////////////////
  int count=0;
  char num[32];
  int number;
  int flag = 0;
  move(h-2,0);
  printw("%c",(char)key);
  num[count]=(char)key;
  echo();
  while((k2=getch()))                        // 複数桁の数値と繰り返す実行内容を入力するためのループ
   {   
   if((k2 < '0') || (k2 > '9'))              // 数値に続いて数値以外  // MARK
    {   
    if(k2 == ':')                            // 数値に続いて: なら左：染色体番号、右：延期位置　に移動
     {                                       // i.e. 5:100000 ５番染色体の100000塩基目
     count++;                                //       :300    今いる染色体の300塩基目
     num[count]='\0';                        //      6:       ６番染色体の１塩基目
     newchr=atoi(num);                       // 移動先染色体
     getstr(com_str);
     int an = arenumbers(com_str);
     if(an == -99999)
      break;
     else
      if(an <= 0)
       newpos = 1;
      else
       newpos = an;                           // 移動先塩基位置
     if(logout == 1)
      fprintf(logfile,"Move to Chr %10d Pos %10d \n",newchr,newpos);

      if(current_chr == newchr)                  // 同一染色体内
       {   
       move_to_pos(newpos);
       }   
      else
       {   
       vmode=0;
       moveto=newchr;
       if(moveto >= n_ref)
        moveto = n_ref;
       if(moveto < 1)
        moveto = 1;
       noecho();
       cbreak();
       current_chr = moveto-1;
       genomelen = chrs[current_chr].seq_len;  //　塩基長を移動先染色体の長さにアップデート
       move_to_pos(newpos); 
       }  
     break;
     }
    if((k2 == 'l')||(k2 == (char)KEY_RIGHT)) // 数値に続いてl なら回数右
     {   
     count++;
     num[count]='\0';
     number=atoi(num);
     if(logout == 1)
      fprintf(logfile,"CPOS %10d %10d base Right\n",cpos,number);
     moveto = cpos + number;
     if(moveto > genomelen)
      moveto = genomelen+1;
     if(logout == 1)
      fprintf(logfile,"moveto %10d genomelen %10d\n",moveto, genomelen);
     flag = 1;
     }   
    if((k2 == 'j')||(k2 == (char)KEY_DOWN))  // 数値に続いてj なら回数下
     {   
     count++;
     num[count]='\0';
     number=atoi(num);
     if(logout == 1)
      fprintf(logfile,"CPOS %10d %10d base Down\n",cpos,number);
     moveto = cpos + (number * w); 
     if(moveto > genomelen)                  
      {
      int step = (genomelen - cpos) / w;
      moveto = cpos + (step * w);
      if(logout == 1)
       moveto = genomelen+1;             // MARKMV
      }
     if(logout == 1)
      fprintf(logfile,"moveto %10d genomelen %10d\n",moveto, genomelen);
     flag = 1;
     }   
    if((k2 == 'k')||(k2 == (char)KEY_UP))    // 数値に続いてk なら回数上
     {   
     count++;
     num[count]='\0';
     number=atoi(num);
     if(logout == 1)
      fprintf(logfile,"CPOS %10d %10d base Up\n",cpos,number);
     moveto = cpos - (number * w);
     if(moveto < 1)
      {
      int step = (cpos-1) / w;
      moveto = cpos - (step * w);
      }
     if(logout == 1)
      fprintf(logfile,"moveto %10d genomelen %10d w %8d\n",moveto, genomelen,w);
     flag = 1;
     }   
    if((k2 == 'h')||(k2 == (char)KEY_LEFT))  // 数値に続いてh なら回数左
     {   
     count++;
     num[count]='\0';
     number=atoi(num);
     if(logout == 1)
      fprintf(logfile,"CPOS %10d %10d base Left\n",cpos,number);
     moveto = cpos - number;
     if(logout == 1)
      fprintf(logfile,"moveto %10d genomelen %10d\n",moveto, genomelen);
     flag = 1;
     }   
    break;                                   // 数値に続いて数値以外の処理はここでブレーク-> 下のFLAGへ抜ける
    }   
   else                                      // 数値に続いて数字なら桁数を増やしてループを回す
    {   
    count++;
    num[count]=(char)k2;
    }   
   }   
  if(flag == 1)                              // FLAG
   {   
   if(moveto <= 0)
    moveto=1;
   if(moveto > genomelen+1)
    moveto=genomelen+1;
   move_to_pos(moveto);                      // movetoの塩基位置へカーソル移動
   refresh();
   }   
  noecho();
  } /////////////////////////////////////////////////////////////////////////////////////////
 else
  {
  switch (key)                               // 数値以外のキーストロークへの対応
   {
   case ':': move(h-2,0);                    // 飛び先塩基位置指定  染色体番号を省略  :100  同一染色体内　塩基位置100へ
             printw(":");
             echo();
             nocbreak();
             getstr(com_str);
             if(com_str[0]=='q')             // :q = プログラムの終了
              return 0;
             moveto=atoi(com_str);
             if(moveto < 1)
              {
              moveto = 1;
              }
             else
              {
              if(moveto > genomelen)
               {
               moveto = genomelen;
               }
              }
             noecho();
             cbreak();
             move_to_pos(moveto);
             refresh();
             break;
  
   case '\'':                             // ' シングルクォート(+alphabet)　マーク位置へジャンプ
             k2 = getch();
             if(mark_chr[k2] == -1)
              break;
  
             newchr = mark_chr[k2];
             newpos = mark_pos[k2];
  
             if(current_chr == newchr)
              {
              moveto=newchr;
              if(moveto >= n_ref)
               moveto = n_ref;
              if(moveto < 1)
               moveto = 1;
              noecho();
              cbreak();
              current_chr = moveto-1;
              moveto = newpos;
              move_to_pos(moveto);
              }
             else
              {
              vmode=0;
              moveto=newchr;
              if(moveto >= n_ref)
               moveto = n_ref;
              if(moveto < 1)
               moveto = 1;
              noecho();
              cbreak();
              current_chr = moveto-1;
              moveto = newpos;
              show_origin = ((moveto/w) * w)+1;
              cy = 1;
              if(moveto%w == 0)
               cx = w-1;
              else
              cx = (moveto%w)-1;
              }
             refresh();
             break;
  
   case '`':                             // ' バッククォート(+alphabet)　ビジュアルマーク位置開始点へジャンプ
             k2 = getch();
             if(vic[k2] == -1)
              break;
  
             newchr = vic[k2];
             newpos = vip[k2];
  
             if(current_chr == newchr)
              {
              moveto=newchr;
              if(moveto >= n_ref)
               moveto = n_ref;
              if(moveto < 1)
               moveto = 1;
              noecho();
              cbreak();
              current_chr = moveto-1;
              moveto = newpos;
              move_to_pos(moveto);
              }
             else
              {
              vmode=0;
              moveto=newchr;
              if(moveto >= n_ref)
               moveto = n_ref;
              if(moveto < 1)
               moveto = 1;
              noecho();
              cbreak();
              current_chr = moveto-1;
              moveto = newpos;
              show_origin = ((moveto/w) * w)+1;
              cy = 1;
              if(moveto%w == 0)
               cx = w-1;
              else
              cx = (moveto%w)-1;
              }
             refresh();
             break;

   case '/': move(h-2,0);                    // 全体サーチ  全体の初めの位置に飛ぶ
             printw("/");
             echo();
             nocbreak();
             getstr(com_str);

             if(com_str[0] == 'v')           // MARKV
              {
              char vmark;
              if((com_str[1] >= '0') && (com_str[1] <= '9'))
               {
               mt = com_str[1] - '0';
               vmark = com_str[2];
               }
              else
               {
               mt = 0;
               vmark = com_str[1];
               }
              if(vic[vmark] == -1)
               {
               sprintf(message,"You specified undefined Visual Region %c\n",vmark);
               noecho();
               cbreak();
               refresh();
               break;
               }
              else
               {
               int vc = 0;
               for(i=vip[vmark]-1;i<vep[vmark];i++)
                {
                com_str[vc++] = chrs[vic[vmark]].seq[i];
                }
               com_str[vc] = '\0';
               }
              }
             else
              {
              mt = splitnum(com_str);         // 数値が１桁1~9なら曖昧検索のミスマッチ許容数として対応
  
              if(mt == 999)
               {
               strcpy(message,"Mismatch number should be 1~9\n");
               noecho();
               cbreak();
               refresh();
               break;
               }
              if(mt == 9999)                  // 検索文字列がアルファベット以外を含む
               {
               sprintf(message,"Query includes characters not included in Reference fasta: %s\n",com_str);
               noecho();
               cbreak();
               refresh();
               break;
               }
              }
             if(strlen(com_str) <= 10)        // 探索文字列が10塩基以下であれば単純サーチ、10以上ならマルチスレッドで全体サーチ
              {
              plain_fsearch(com_str);
              search_mode = 1;
              }
             else
              {
              seq_search(com_str,1);
              search_mode = 3;
              }
             strcpy(search_str,com_str);
             noecho();
             cbreak();
             refresh();
             break;
  
   case '?': move(h-2,0);                    // 全体サーチ  全体の初めの位置に飛ぶ
             printw("?");
             echo();
             nocbreak();
             getstr(com_str);
  
             if(com_str[0] == 'v')
              {
              char vmark;
              if((com_str[1] >= '0') && (com_str[1] <= '9'))
               {
               mt = com_str[1] - '0';
               vmark = com_str[2];
               }
              else
               {
               mt = 0;
               vmark = com_str[1];
               }
              if(vic[vmark] == -1)
               {
               sprintf(message,"You specified undefined Visual Region %c\n",vmark);
               noecho();
               cbreak();
               refresh();
               break;
               }
              else
               {
               int vc = 0;
               for(i=vip[vmark]-1;i<vep[vmark];i++)
                {
                com_str[vc++] = chrs[vic[vmark]].seq[i];
                }
               com_str[vc] = '\0';
               }
              }
             else
              {
              mt = splitnum(com_str);         // 数値が１桁1~9なら曖昧検索のミスマッチ許容数として対応
  
              if(mt == 999)
               {
               strcpy(message,"Mismatch number should be 1~9\n");
               noecho();
               cbreak();
               refresh();
               break;
               }
              if(mt == 9999)                  // 検索文字列がアルファベット以外を含む
               {
               sprintf(message,"Search string includes non-alphabet character: %s\n",com_str);
               noecho();
               cbreak();
               refresh();
               break;
               }
              }
             if(strlen(com_str) <= 10)        // 探索文字列が10塩基以下であれば単純サーチ、
                                              // 10以上ならマルチスレッドで全体サーチ
              {
              plain_bsearch(com_str);
              search_mode = 2;
              }
             else
              {
              seq_search(com_str,-1);
              search_mode = 4;
              }
  
             strcpy(search_str,com_str);
             noecho();
             cbreak();
             refresh();
             break;
  
   case 'n':                                // サーチ繰り返し 前方
             if(search_mode == 0)
              break;
             if((search_mode == 1) || (search_mode == 2))
              {
              plain_fsearch(search_str);
              break;
              }
  
             fflag = 0; 
             for(i=0;i<n_hit_pos-1;i++)
              {
              int fchr = hpis[i].chr;
              int fpos = hpis[i].pos;
              int tchr = hpis[i+1].chr;
              int tpos = hpis[i+1].pos;
              int cchr = current_chr;
              cpos = show_origin + (cy-1) * w + cx;                           // cpos       : カーソル位置 
              if(cchr > fchr)        // 現在の染色体番号がiの染色体番号より大きい 
               {
               if(cchr < tchr)       // 現在の染色体番号がiの染色体番号より小さい 
                {    
                hit_point = i+1; 
                fflag = 1; 
                break;
                }    
               else 
                {    
                if(cchr == tchr)
                 {    
                 if(cpos <= tpos)
                  {    
                  hit_point = i+1; 
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
                  hit_point = i+1; 
                  fflag = 1; 
                  break;
                  }    
                 else 
                  {    
                  if(cpos <= tpos)
                   {    
                   hit_point = i+1; 
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
             if(current_chr != hpis[hit_point].chr)
              vmode = 0;
             current_chr = hpis[hit_point].chr;
             moveto = hpis[hit_point].pos;
  
             move_to_pos(moveto+1);
             refresh();
             break;
  
   case 'N':                                // サーチ繰り返し 後方MARK
             if(search_mode == 0)
              break;
             if((search_mode == 1) || (search_mode == 2))
              {
              plain_bsearch(search_str);
              break;
              }
  
             fflag = 0; 
             for(i=0;i<n_hit_pos-1;i++)
              {    
              int fchr = hpis[i].chr;
              int fpos = hpis[i].pos;
              int tchr = hpis[i+1].chr;
              int tpos = hpis[i+1].pos;
              int cchr = current_chr;
              cpos = show_origin + (cy-1) * w + cx;  // cpos       : カーソル位置 
              if(cchr > fchr)        // 現在の染色体番号がiの染色体番号より大きい 
               {    
               if(cchr < tchr)       // 現在の染色体番号がi+1の染色体番号より小さい 
                {    
                hit_point = i; 
                fflag = 1; 
                break;
                }    
               else 
                {    
                if(cchr == tchr)     // 現在の染色体番号とi+1の染色体番号が同じ
                 {    
                 if(cpos <= tpos+1)   // 現在の塩基位置がi+1の塩基位置より低い
                  {    
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
                if(cpos > fpos+1)    // +1
                 {    
                 if(cchr < tchr)
                  {    
                  hit_point = i; 
                  fflag = 1; 
                  break;
                  }    
                 else               // i+1の染色体番号とiの染色体番号と現在の染色体番号が同じ 
                  {    
                  if(cpos <= tpos+1)  // 
                   {    
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
              hit_point = n_hit_pos-1; 
              }    
  
             if(current_chr != hpis[hit_point].chr)
              vmode = 0;
             current_chr = hpis[hit_point].chr;
             moveto = hpis[hit_point].pos+1;
  
             move_to_pos(moveto);

             refresh();
             break;
  
   case 'r':                                      // 検索結果等のアトリビュート（強調表示）をリセット
             attrset(A_REVERSE);                           // 表示反転 
             move(0,0);                                    // 天辺行
             cc = w - 35;                                  // 右端から左に35文字戻る
             sprintf(buff,"R Chr: %8d/%8d             %*s ",current_chr+1,n_ref,cc,chrs[current_chr].title);
             buff[w] = '\0';
             printw("%s",buff);
             attrset(A_NORMAL);                           // 表示反転 

             k2 = getch();
             if((k2 == 'r') || (k2 == 'a'))
              {
              if(hpis != NULL)
               {
               free(hpis);
               hpis = NULL;
               }
              atr_reset();
              refresh();
              break;
              }
             if(k2 == 'v')                        // visualマーク位置 とそのアトリビュートをリセット
              {
              for(i='A';i<='z';i++)
               {
               if(vic[i] != -1)              // visualマーク されているアルファベットであれば
                {
                for(j=vip[i];j<=vep[i];j++)
                 {
                 if((chrs[vic[i]].atr[j] & BIT_FLAG_5) != 0)
                  chrs[vic[i]].atr[j] -= BIT_FLAG_5;
                 }
                vic[i] = -1;
                vip[i] = -1;
                vep[i] = -1;
                }
               }
              refresh();
              break;
              }
             if(k2 == 'm')                        // MARKマーク位置 とそのアトリビュートをリセット
              {
              for(i='A';i<='z';i++)
               {
               if(mark_chr[i] != -1)                                // mark されているアルファベットであれば
                if((chrs[mark_chr[i]].atr[mark_pos[i]] & BIT_FLAG_6) != 0)
                 chrs[mark_chr[i]].atr[mark_pos[i]] -= BIT_FLAG_6;
               mark_chr[i] = -1;
               mark_pos[i] = -1;
               }
              refresh();
              break;
              }
             if((k2 == 's') || (k2 == 'f'))       // 検索結果 Search/Find のアトリビュートをリセット
              {
              if(hpis != NULL)
               {
               free(hpis);
               hpis = NULL;
               }
              atr_reset();
              for(i='A';i<='z';i++)
               {
               if(mark_chr[i] != -1)                                // mark されているアルファベットであれば
                chrs[mark_chr[i]].atr[mark_pos[i]] = chrs[mark_chr[i]].atr[mark_pos[i]] | BIT_FLAG_6;
               if(vic[i] != -1)                                     // visual設定 されているアルファベットであれば
                for(j=vip[i];j<=vep[i];j++)
                 chrs[vic[i]].atr[j] = chrs[vic[i]].atr[j] | BIT_FLAG_5;
               }
              refresh();
              break;
              }
             refresh();
             break;
  
   case 'G': moveto=genomelen;              // G 染色体内最後尾への移動
             move_to_pos(moveto);
             break;
  
   case 'g':                                       // gg 染色体内１塩基目への移動
              k2 = getch();
              if(k2 == 'g')
               {
               moveto = 1;
               show_origin = 1;
               cx = 0;
               cy = 1;
               }
             break;
  
   case 'm':                         // マーク　MARK Position   ma~mz, mA~mZ 
                                     //can jump back to the position by 'a~'z, 'A~'Z
             attrset(A_REVERSE);                           // 表示反転 
             move(0,0);                                    // 天辺行
             cc = w - 35;                                  // 右端から左に35文字戻る
             sprintf(buff,"M Chr: %8d/%8d             %*s ",current_chr+1,n_ref,cc,chrs[current_chr].title);
             buff[w] = '\0';
             printw("%s",buff);
             attrset(A_NORMAL);                           // 表示反転 

             k2 = getch();
             if(((k2 >= 'a') && (k2 <= 'z')) ||
                ((k2 >= 'A') && (k2 <= 'Z')))
              {
              if(mark_chr[k2] != -1)                                // mark されているアルファベットであれば
               if((chrs[mark_chr[k2]].atr[mark_pos[k2]] & BIT_FLAG_6) != 0)
                chrs[mark_chr[k2]].atr[mark_pos[k2]] -= BIT_FLAG_6;  // もとの位置のアトリビュートはリセット
              mark_chr[k2] = current_chr;
              mark_pos[k2] = show_origin + ((cy-1)*w) + cx;
              chrs[mark_chr[k2]].atr[mark_pos[k2]] = chrs[mark_chr[k2]].atr[mark_pos[k2]] | BIT_FLAG_6;  
                                                                                      // BIT_FLAG_6 を立てる
              }
            break;
  
   case 'v':                         // ビジュアルモード　領域選択 Enter Visual Mode
            if(vmode == 0)
              {
              vmode = 1;
              vichr  = current_chr;
              vipos  = show_origin + ((cy-1)*w) + cx;
              }
            else
             {
             if(vmode == 1)
              {
              if(current_chr == vichr)
               {
               k2 = getch();
               if(((k2 >= 'a') && (k2 <= 'z')) ||
                  ((k2 >= 'A') && (k2 <= 'Z')))
                {
                if(vic[k2] != -1)          // すでにビジュアルで指定されているアルファベットの場合
                 {
                 for(i=vip[k2];i<=vep[k2];i++)   // 元のアトリビュートをリセット
                  {
                  chrs[vic[k2]].atr[i] -= BIT_FLAG_5;
                  }
                 }
                int pospos = show_origin + ((cy-1)*w) + cx;
                if(pospos > vipos)
                 {
                 for(i=vipos;i<=pospos;i++)
                  chrs[vichr].atr[i] = chrs[vichr].atr[i] | BIT_FLAG_5;
                 vepos = pospos;
                 }
                else
                 {
                 for(i=pospos;i<=vipos;i++)
                  chrs[vichr].atr[i] = chrs[vichr].atr[i] | BIT_FLAG_5;
                 vepos = vipos;
                 vipos = pospos;
                 }
                vic[k2] = vichr;
                vip[k2] = vipos;
                vep[k2] = vepos;
                }
               }
              vmode = 0;
              }
             }
            break;
  
   case KEY_SRIGHT :
   case 'L' : current_chr++;                      // 染色体番号インクリメント
              if(current_chr>= n_ref)
               current_chr= n_ref-1;
              move_to_pos(1);
              vmode=0;
              break;
  
   case KEY_SLEFT :
   case 'H' : current_chr--;                      // 染色体番号デクリメント
              if(current_chr< 0)
               current_chr= 0;
              move_to_pos(1);
              vmode=0;
              break;
  
//   case 'q' : return(0);                       // gnmを終了

   int apage;
   case ctrl('b') : apage = h-5;            // Page Up
                    moveto = cpos - (apage*w);
                    if (moveto < 1)
                     moveto = 1;
                    move_to_pos(moveto);
                    break;

   case ctrl('f') : apage = h-5;            // Page Down
                    moveto = cpos + (apage*w);
                    if (moveto > genomelen)
                     moveto = genomelen;
                    move_to_pos(moveto);
                    break;

   int half;
   case ctrl('u') : half = (h-4)/2;            // HALF Page Up
                    moveto = cpos - (half*w);
                    if (moveto < 1)
                     moveto = 1;
                    move_to_pos(moveto);
                    break;

   case ctrl('d') : half = (h-4)/2;            // HALF Page Down
                    moveto = cpos + (half*w);
                    if (moveto > genomelen)
                     moveto = genomelen;
                    move_to_pos(moveto);
                    break;
  
   case KEY_LEFT:
   case 'h' : cx --;                      // カーソル左へ
             if(cx < 0)             // 左端
              {
               if(cy == 1)          // 上端   
                {
                show_origin -= w;
                if(show_origin < 1)
                 {
                 show_origin += w;
                 cx = 0;
                 cy = 1;
                 }
                else
                 {
                 cx = w-1;
                 }
                }
               else
                {
                cx = w-1; 
                cy --;
                }
              }
             break;
   
   case KEY_RIGHT:
   case 'l': cx ++;                      // カーソル右へ
             if(show_origin + ((cy-1)*w) + cx > genomelen)  // ゲノム最後尾より右には行けない
              {
              cx --;
              break;
              }
             if(cx > w-1)         // 右端
              {
              if(cy > h-5)        // 下端
               {
               show_origin += w;
               cx = 0;
               }
              else
               {
               cx = 0; 
               cy ++;
               }
              }
              break;
  
   case KEY_UP:
   case 'k': cy--;                       // カーソル上へ
             if(cy < 1)            // 上端
              {
              show_origin -= w; 
              cy = 1;
              if(show_origin < 1)  // スクロールできない
               {
               show_origin += w;
               }
              } 
             break;
  
   case KEY_DOWN:
   case 'j': cy++;                       // カーソル下へ
             if((show_origin + ((cy-1)*w) + cx) > genomelen)
              {
               {
               cy--;
               }
             // cx = (genomelen % w)-1;
              }
             if(cy > (h-4))        // 下端
              {
              show_origin += w; 
              cy=h-4;
              if(show_origin > max_so)
               {
               show_origin -= w;
               }
              }
             break;

   default: 
             break;
   }
  }
////////////////////////////////////////////////////////////////////////キー入力レスポンス
//////////////////////////////////////////////// MAIN LOOP
//////////////////////////////////////////////// MAIN LOOP
 }
//////////////////////////////////////////////// MAIN LOOP
//////////////////////////////////////////////// MAIN LOOP
return (0);
}
