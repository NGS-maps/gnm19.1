#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

extern char arg1[];
extern char arg2[];

extern int  read_gff;
extern char gff_flnm[];


void errorm(int en);

void readargs(int argc, char **argv)
{
int i;
int temp;
int arg;
int flag = 0;
char tempchar[80];

if(argc==1)
 {
 errorm(0);
 }

for(arg=1;arg<argc;arg++)
 {
 if(argv[arg][0] == '-')
  {
///////////////////////////  Options
  if((argv[arg][1] == 'g') && (argv[arg][2] == 'f') &&(argv[arg][3] == 'f'))    
                                   // -gff gfffilename : specify gff file name
   {
   if(argv[arg][4] != '\0')
    errorm(1);
   else
    {
    read_gff  = 1;

    arg ++;
    if(argv[arg][0] == '-')
     errorm(2);
    strcpy(gff_flnm,argv[arg]);
    }
   continue;
   }

  }
 else
  {
  if(flag == 0)
   {                               // read first argument without - as the reference fasta file name
   strcpy(arg1,argv[arg]);
   }
  else
   {                               // read second argument without - as the read fastaq file name
   if(flag == 1)
    strcpy(arg2,argv[arg]);
   }
  flag ++;
  if(flag >= 3)
   errorm(2);
  }
 }
}

void errorm(int en)
 {
 printf("ERROR %d\n",en);
 printf("Usage:   gnm sequence.fasta\n");
 printf("             a nucleic acid multi fasta file\n");
 exit(1);
 }
