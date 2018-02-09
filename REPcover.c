/*******************************************************************************************
 *
 *  Read the .covr track of each db or db block on the command line and output a histogram
 *    of the coverage of the unmasked portions and an estimate of the coverage of the
 *    underlying genome
 *
 *  Author:  Gene Myers
 *  Date  :  January 2017
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

static DAZZ_DB _DB, *DB  = &_DB;   //  Data base

int main(int argc, char *argv[])
{ int c;

  //  Process arguments

  Prog_Name = Strdup("REPcover","");

  if (argc < 2)
    { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
      exit (1);
    }

  //  Open trimmed DB and the qual-track

  for (c = 1; c < argc; c++)
    { int         status;
      char       *root;
      int         i, cmax, hgap_min, cover;
      int64       nreads, totlen;
      int64      *cgram;
      int64       ssum, stotal;
      DAZZ_EXTRA  ex_hgap, ex_covr;

      status = Open_DB(argv[c],DB);
      if (status < 0)
        exit (1);
      if (status == 1)
        { fprintf(stderr,"%s: Cannot be called on a .dam index: %s\n",Prog_Name,argv[1]);
          exit (1);
        }
      Trim_DB(DB);

      { FILE *afile;
        char *aname;
        int   extra;

        afile = NULL;
        if (DB->part)
          { aname = Strdup(Catenate(DB->path,Numbered_Suffix(".",DB->part,"."),"covr",".anno"),
                           "Allocating anno file");
            if (aname == NULL)
              exit (1);
            afile  = fopen(aname,"r");
          }
        if (afile == NULL)
          { aname = Strdup(Catenate(DB->path,".","covr",".anno"),"Allocating anno file");
            if (aname == NULL)
              exit (1);
            afile  = fopen(aname,"r");
          }
        if (afile == NULL)
          { fprintf(stderr,"%s: Must have a 'covr' track, run DAScover\n",Prog_Name);
            exit (1);
          }

        fseeko(afile,0,SEEK_END);
        extra = ftell(afile) - sizeof(int)*2;
        fseeko(afile,-extra,SEEK_END);
        ex_covr.nelem = 0;
        if (Read_Extra(afile,aname,&ex_covr) != 0)
          { fprintf(stderr,"%s: Histogram extra missing from .covr track?\n",Prog_Name);
            exit (1);
          }
        ex_hgap.nelem = 0;
        if (Read_Extra(afile,aname,&ex_hgap) != 0)
          { fprintf(stderr,"%s: Hgap threshold extra missing from .covr track?\n",Prog_Name);
            exit (1);
          }
        fclose(afile);
      }

      root   = Root(argv[c],".db");
      nreads = DB->nreads;
      totlen = DB->totlen;

      hgap_min = (int) ((int64 *) (ex_hgap.value))[0];
      cgram    = (int64 *) (ex_covr.value);
      cmax     = ex_covr.nelem - 1;

      printf("\nDAScover");
      if (hgap_min > 0)
        printf(" -H%d",hgap_min);
      printf(" %s\n\n",root);

      if (hgap_min > 0) 
        { for (i = 0; i < DB->nreads; i++)
            if (DB->reads[i].rlen < hgap_min)
              { nreads -= 1;
                totlen -= DB->reads[i].rlen;
              }
        }

      printf("\nInput:  ");
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

      stotal = 0;
      for (i = 0; i <= cmax; i++)
        stotal += cgram[i];

      printf("\nCoverage Histogram\n\n");
      ssum = cgram[cmax];
      if (ssum > 0)
        printf("    %4d:  %9lld  %5.1f%%\n\n",
               cmax,cgram[cmax],(100.*ssum)/stotal);
      stotal -= ssum;
      ssum    = 0;
      for (i = cmax-1; i >= 0; i--)
        if (cgram[i] > 0)
          { ssum += cgram[i];
            printf("    %4d:  %9lld  %5.1f%%\n",
                   i,cgram[i],(100.*ssum)/stotal);
          }

      i = 0;
      while (cgram[i+1] < cgram[i])
        i += 1;
      for (cover = i++; i < cmax; i++)
        if (cgram[cover] < cgram[i])
          cover = i;
        
      printf("\n  Coverage is estimated at %d\n\n",cover);

      free(root);
      Close_DB(DB);
    }

  free(Prog_Name);
  exit (0);
}
