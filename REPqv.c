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

    { DAZZ_DB    _DB, *DB  = &_DB;
      DAZZ_EXTRA ex_hgap, ex_cest, ex_qvs, ex_dif;

      //  Load DB

      { int status;

        status = Open_DB(argv[c],DB);
        if (status < 0)
          exit (1);
        if (status == 1)
          { fprintf(stderr,"%s: Cannot be called on a .dam index: %s\n",Prog_Name,argv[1]);
            exit (1);
          }
        Trim_DB(DB);
      }

      //  Get .qual track extras

      { FILE *afile;
        char *aname;
        int   extra, tracklen, size;

        afile = NULL;
        if (DB->part)
          { aname = Strdup(Catenate(DB->path,Numbered_Suffix(".",DB->part,"."),"qual",".anno"),
                           "Allocating anno file");
            if (aname == NULL)
              exit (1);
            afile = fopen(aname,"r");
            if (afile == NULL)
              { fprintf(stderr,"%s: Must have a 'qual.%d' track, run DASqv\n",Prog_Name,DB->part);
                exit (1);
              }
          }
        else
          { aname = Strdup(Catenate(DB->path,".","qual",".anno"),"Allocating anno file");
            if (aname == NULL)
              exit (1);
            afile  = fopen(aname,"r");
            if (afile == NULL)
              { fprintf(stderr,"%s: Must have a 'qual' track, run DASqv\n",Prog_Name);
                exit (1);
              }
          }

        fread(&tracklen,sizeof(int),1,afile);
        fread(&size,sizeof(int),1,afile);
        fseeko(afile,0,SEEK_END);
        extra = ftell(afile) - (size*(tracklen+1) + 2*sizeof(int));
        fseeko(afile,-extra,SEEK_END);
        ex_hgap.nelem = 0;
        if (Read_Extra(afile,aname,&ex_hgap) != 0)
          { fprintf(stderr,"%s: Hgap threshold extra missing from .qual track?\n",Prog_Name);
            exit (1);
          }
        ex_cest.nelem = 0;
        if (Read_Extra(afile,aname,&ex_cest) != 0)
          { fprintf(stderr,"%s: Coverage estimate extra missing from .qual track?\n",Prog_Name);
            exit (1);
          }
        ex_qvs.nelem = 0;
        if (Read_Extra(afile,aname,&ex_qvs) != 0)
          { fprintf(stderr,"%s: QV histogram extra missing from .qual track?\n",Prog_Name);
            exit (1);
          }
        ex_dif.nelem = 0;
        if (Read_Extra(afile,aname,&ex_dif) != 0)
          { fprintf(stderr,"%s: Differences histogram extra missing from .qual track?\n",Prog_Name);
            exit (1);
          }
        fclose(afile);
      }

      //  Generate display

      { char       *root;
        int64       nreads, totlen;
        int         hgap_min, cover;
        int64      *qgram, *sgram;
        int         maxqv;

        //  Get relevant variables

        root   = Root(argv[c],".db");
        nreads = DB->nreads;
        totlen = DB->totlen;

        hgap_min = (int) ((int64 *) (ex_hgap.value))[0];
        cover    = (int) ((int64 *) (ex_cest.value))[0];
        qgram    = (int64 *) (ex_qvs.value);
        maxqv    = ex_qvs.nelem - 1;
        sgram    = (int64 *) (ex_dif.value);

        printf("\nDASqv");
        if (hgap_min > 0)
          printf(" -H%d",hgap_min);
        printf(" -c%d %s\n\n",cover,root);

        if (hgap_min > 0) 
          { int i;

            for (i = 0; i < DB->nreads; i++)
              if (DB->reads[i].rlen < hgap_min)
                { nreads -= 1;
                  totlen -= DB->reads[i].rlen;
                }
          }

        //  Display histograms

        printf("\n  Input:  ");
        Print_Number(nreads,7,stdout);
        printf("reads,  ");
        Print_Number(totlen,12,stdout);
        printf(" bases");
        if (hgap_min > 0)
          { printf(" (another ");
            Print_Number(DB->nreads - nreads,0,stdout);
            printf(" were < H-length)");
          }
        printf("\n");

        { int64 ssum, qsum;
          int64 stotal, qtotal;
          int   qv_deep;
          int   gval, bval;
          int   i;

          stotal = qtotal = 0;
          for (i = 0; i <= maxqv; i++)
            { stotal += sgram[i];
              qtotal += qgram[i];
            }

          if (cover >= 40)
            qv_deep = cover/8;
          else if (cover >= 20)
            qv_deep = 5;
          else
            qv_deep = cover/4;

          printf("\n  Histogram of q-values (average %d best)\n",qv_deep);
          printf("\n                 Input                 QV\n");
          qsum = qgram[maxqv];
          ssum = sgram[maxqv];
          printf("\n    %2d:  %9lld  %5.1f%%    %9lld  %5.1f%%\n\n",
                 maxqv,sgram[maxqv],(100.*ssum)/stotal,qgram[maxqv],(100.*qsum)/qtotal);

          qtotal -= qsum;
          stotal -= ssum;
          ssum = qsum = 0;
          for (i = maxqv-1; i >= 0; i--)
            if (qgram[i] > 0)
              { ssum += sgram[i];
                qsum += qgram[i];
                printf("    %2d:  %9lld  %5.1f%%    %9lld  %5.1f%%\n",
                       i,sgram[i],(100.*ssum)/stotal,
                         qgram[i],(100.*qsum)/qtotal);
              }

          //  Estimate -g and -b parameters

          bval = gval = -1;
          qsum = 0;
          for (i = maxqv-1; i >= 0; i--)
            if (qgram[i] > 0)
              { qsum += qgram[i];
                if ((100.*qsum)/qtotal > 7. && bval < 0)
                  bval = i+1;
                if ((100.*qsum)/qtotal > 20. && gval < 0)
                  gval = i+1;
              }
  
          printf("\n  Recommend \'DAStrim -g%d -b%d'\n\n",gval,bval);
        }

        free(root);
        Close_DB(DB);
      }
    }

  free(Prog_Name);
  exit (0);
}
