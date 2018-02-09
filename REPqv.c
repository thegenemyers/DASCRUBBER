/*******************************************************************************************
 *
 *  Read the .qual track of each db or db block on the command line and output a histogram
 *    of the intrinsic QV's therein.
 *
 *  Author:  Gene Myers
 *  Date  :  August 2017
 *
 *******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include "DB.h"

//  Command format and global parameter variables

static char *Usage = "<source:db> ...";

//  Global Variables

#define  MAXQV   50     //  Max QV score is 50
#define  MAXQV1  51

static DAZZ_DB _DB, *DB  = &_DB;   //  Data base

static int64  *QV_IDX;     //  qual track index
static uint8  *QV;         //  qual track values

int main(int argc, char *argv[])
{ int c;

  //  Process arguments

  Prog_Name = Strdup("REPqv","");

  if (argc < 2)
    { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
      exit (1);
    }

  //  Open trimmed DB and the qual-track

  for (c = 1; c < argc; c++)
    { int         status;
      char       *root;
      int         i, bval, gval, cover, hgap_min;
      DAZZ_TRACK *track;
      int64       nreads, totlen;
      int64       qgram[MAXQV1];
      int64       qsum, qtotal;

      status = Open_DB(argv[c],DB);
      if (status < 0)
        exit (1);
      if (status == 1)
        { fprintf(stderr,"%s: Cannot be called on a .dam index: %s\n",Prog_Name,argv[1]);
          exit (1);
        }
      Trim_DB(DB);

      track = Load_Track(DB,"qual");
      if (track != NULL)
        { FILE *afile;
          int   size, tracklen, extra;

          QV_IDX = (int64 *) track->anno;
          QV     = (uint8 *) track->data;

          //  if newer .qual tracks with -c meta data, grab it

          if (DB->part)
            { afile = fopen(Catenate(DB->path,
                                  Numbered_Suffix(".",DB->part,"."),"qual",".anno"),"r");
              if (afile == NULL)
                afile = fopen(Catenate(DB->path,".","qual",".anno"),"r");
            }
          else
            afile = fopen(Catenate(DB->path,".","qual",".anno"),"r");
          fread(&tracklen,sizeof(int),1,afile);
          fread(&size,sizeof(int),1,afile);
          fseeko(afile,0,SEEK_END);
          extra = ftell(afile) - (size*(tracklen+1) + 2*sizeof(int));
          fseeko(afile,-extra,SEEK_END);
          if (extra == 2*sizeof(int))
            { fread(&cover,sizeof(int),1,afile);
              fread(&hgap_min,sizeof(int),1,afile);
            }
          else if (extra == sizeof(int))
            { fread(&cover,sizeof(int),1,afile);
              hgap_min = 0;
            }
          else
            { cover = -1;
              hgap_min = 0;
            }
          fclose(afile);
        }
      else
        { fprintf(stderr,"%s: Must have a 'qual' track, run DASqv\n",Prog_Name);
          exit (1);
        }

      root   = Root(argv[c],".db");
      nreads = DB->nreads;
      totlen = DB->totlen;

      for (i = 0; i <= MAXQV; i++)
        qgram[i] = 0;

      for (i = 0; i < QV_IDX[nreads]; i++)
        qgram[QV[i]] += 1;

      printf("\nDASqv");
      if (cover >= 0)
        printf(" -c%d",cover);
      else
        printf(" -c??");
      if (hgap_min > 0)
        printf(" [-H%d]",hgap_min);
      printf(" %s\n\n",root);

      if (hgap_min > 0) 
        { for (i = 0; i < DB->nreads; i++)
            if (DB->reads[i].rlen < hgap_min)
              { nreads -= 1;
                totlen -= DB->reads[i].rlen;
              }
        }
      printf("  Input:  ");
      Print_Number(nreads,7,stdout);
      printf("reads,  ");
      Print_Number(totlen,12,stdout);
      printf(" bases");
      if (hgap_min > 0) 
        { printf(" (another ");
          Print_Number(DB->nreads-nreads,0,stdout);
          printf(" were < H-length)");
        }
      printf("\n");

      if (cover >= 0)
        { int qv_deep;

          if (cover >= 40)
            qv_deep = cover/8;
          else if (cover >= 20)
            qv_deep = 5;
          else
            qv_deep = cover/4;

          printf("\n  Histogram of q-values (average %d best)\n",2*qv_deep);
        }
      else
        printf("\n  Histogram of q-values\n");

      qtotal = 0;
      for (i = 0; i <= MAXQV; i++)
        qtotal += qgram[i];

      qsum = qgram[MAXQV];
      printf("\n      %2d:  %9lld  %5.1f%%\n\n",MAXQV,qgram[MAXQV],(100.*qsum)/qtotal);

      bval    = gval = -1;
      qtotal -= qsum;
      qsum    = 0;
      for (i = MAXQV-1; i >= 0; i--)
        if (qgram[i] > 0)
          { qsum += qgram[i];
            printf("      %2d:  %9lld  %5.1f%%\n",i,qgram[i],(100.*qsum)/qtotal);
            if ((100.*qsum)/qtotal > 7. && bval < 0)
              bval = i;
            if ((100.*qsum)/qtotal > 20. && gval < 0)
              gval = i;
          }

      printf("\n  Recommend \'DAStrim -g%d -b%d'\n\n",gval,bval);

      free(root);
      Close_DB(DB);
    }

  free(Prog_Name);
  exit (0);
}
