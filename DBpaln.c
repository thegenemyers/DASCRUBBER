/*******************************************************************************************
 *
 *  Scan every read looking for palindromes and record the resulting intervals in a
 *    "paln" track.
 *
 *  Author:  Gene Myers
 *  Date  :  March 2015
 *
 *******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "DB.h"
#include "align.h"
#include "seqpats.h"

#ifdef HIDE_FILES
#define PATHSEP "/."
#else
#define PATHSEP "/"
#endif

#undef  PAL_DEBUG

static char *Usage = " [-v] <source:db|dam>";

static int     VERBOSE;

static HITS_DB _DB, *DB  = &_DB;   //  Data base

static FILE   *PAL_AFILE;   //  .paln.anno
static FILE   *PAL_DFILE;   //  .paln.data
static int64   PAL_INDEX;   //  Current index into .paln.data file

//  Statistics

static int64 npals, nbps;


int main(int argc, char *argv[])
{ char *root, *pwd;
  int   status;

  //  Process arguments

  { int  i, j, k;
    int  flags[128];
    // char *eptr;

    ARG_INIT("DBpaln")

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("v")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];

    if (argc != 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }

    status = Open_DB(argv[1],DB);
    if (status < 0)
      exit (1);
  }

  //   Set up QV and trimming tracks

  pwd   = PathTo(argv[1]);
  if (status)
    root  = Root(argv[1],".dam");
  else
    root  = Root(argv[1],".db");

  PAL_AFILE = Fopen(Catenate(pwd,PATHSEP,root,".paln.anno"),"w");
  PAL_DFILE = Fopen(Catenate(pwd,PATHSEP,root,".paln.data"),"w");
  if (PAL_AFILE == NULL || PAL_DFILE == NULL)
    exit (1);

  free(pwd);
  free(root);

  { int size;

    size = 0;
    fwrite(&(DB->nreads),sizeof(int),1,PAL_AFILE);
    fwrite(&size,sizeof(int),1,PAL_AFILE);
    PAL_INDEX = 0;
    fwrite(&PAL_INDEX,sizeof(int64),1,PAL_AFILE);
  }

  //  Initialize statistics gathering

  if (VERBOSE)
    { npals = 0;
      nbps  = 0;
    }

  //  Process each read

  { int      i, beg, end;
    Pal_Hit *hit;

    for (i = 0; i < DB->nreads; i++)
      {
#ifdef PAL_DEBUG
        int rlen = DB->reads[i].rlen;

        printf("Analyzing %d(%x) %d\n",i,(DB->reads[i].flags & DB_BEST) != 0,rlen);
#endif

        end = -1;
        while ((hit = Pal_Finder(DB,i)) != NULL)
          { beg = hit->midpt - hit->width/2;
            if (beg > end)
              { if (end >= 0)
                  { fwrite(&end,sizeof(int),1,PAL_DFILE);
                    PAL_INDEX += sizeof(int);
                  }
                fwrite(&beg,sizeof(int),1,PAL_DFILE);
                PAL_INDEX += sizeof(int);
              }
            if (hit->midpt + hit->width/2 > end)
              end = hit->midpt + hit->width/2;

            npals += 1;
            nbps  += hit->width;
#ifdef PAL_DEBUG
            printf("  @ %5d: %4d = %2d%%",hit->midpt,hit->width,hit->qvmat);
            if (hit->gap > 0)
              printf(" [%2d]",hit->gap);
            printf("\n");
#endif
          }
        if (end >= 0)
          { fwrite(&end,sizeof(int),1,PAL_DFILE);
            PAL_INDEX += sizeof(int);
          }
        fwrite(&PAL_INDEX,sizeof(int64),1,PAL_AFILE);
      }
  }

  //  If verbose output statistics summary to stdout

  if (VERBOSE)
    { printf("  Found %lld palindromes, totaling ",npals);
      Print_Number(nbps,0,stdout);
      printf(" base pairs\n");
    }

  //  Clean up

  fclose(PAL_AFILE);
  fclose(PAL_DFILE);

  free(Prog_Name);

  exit (0);
}
