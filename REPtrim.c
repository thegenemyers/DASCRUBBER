/*******************************************************************************************
 *
 *  Read the .trim track of each db or db block on the command line and output
 *    a summary of the scrubbing that took place on that db or block.
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
#include "align.h"

//  Command format and global parameter variables

static char *Usage = "<source:db> ...";

//  Gap states

#define LOWQ  0   //  Gap is spanned by many LAs and patchable
#define SPAN  1   //  Gap has many paired LAs and patchable
#define SPLIT 2   //  Gap is a chimer or an unpatchable gap
#define ADPRE 3   //  Gap is due to adaptemer, trim prefix interval to left
#define ADSUF 4   //  Gap is due to adaptemer, trim suffix interval to right



//  Global Variables (must exist across the processing of each pile)

static HITS_DB _DB, *DB  = &_DB;   //  Data base

static int64  *TRIM_IDX;   //  trim track index
static int    *TRIM;       //  trim track values

int main(int argc, char *argv[])
{ int c;

  //  Process arguments

  Prog_Name = Strdup("REPtrim","");

  if (argc < 2)
    { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
      exit (1);
    }

  //  Open trimmed DB and .qual and .trim tracks

  for (c = 1; c < argc; c++)
    { int         status;
      char       *root;
      int         i, a, tb, te;
      int         alen;
      HITS_TRACK *track;
      int64       nreads, totlen;
      int64       nelim, nelimbp;
      int64       n5trm, n5trmbp;
      int64       n3trm, n3trmbp;
      int64       natrm, natrmbp;
      int64       ngaps, ngapsbp;
      int64       nlowq, nlowqbp;
      int64       nspan, nspanbp;
      int64       nchim, nchimbp;
      int         BAD_QV, GOOD_QV, COVERAGE, HGAP_MIN;

      status = Open_DB(argv[c],DB);
      if (status < 0)
        exit (1);
      if (status == 1)
        { fprintf(stderr,"%s: Cannot be called on a .dam index: %s\n",Prog_Name,argv[1]);
          exit (1);
        }
      Trim_DB(DB);

      track = Load_Track(DB,"trim");
      if (track != NULL)
        { FILE *afile;
          int   size, tracklen, extra;

          TRIM_IDX = (int64 *) track->anno;
          TRIM     = (int *) track->data;
          for (i = 0; i <= DB->nreads; i++)
            TRIM_IDX[i] /= sizeof(int);

          if (DB->part)
            { afile = fopen(Catenate(DB->path,
                                  Numbered_Suffix(".",DB->part,"."),"trim",".anno"),"r");
              if (afile == NULL)
                afile = fopen(Catenate(DB->path,".","trim",".anno"),"r");
            }
          else
            afile = fopen(Catenate(DB->path,".","trim",".anno"),"r");
          fread(&tracklen,sizeof(int),1,afile);
          fread(&size,sizeof(int),1,afile);
          fseeko(afile,0,SEEK_END);
          extra = ftell(afile) - (size*(tracklen+1) + 2*sizeof(int));
          fseeko(afile,-extra,SEEK_END);
          if (extra == 4*sizeof(int))
            { fread(&GOOD_QV,sizeof(int),1,afile);
              fread(&BAD_QV,sizeof(int),1,afile);
              fread(&COVERAGE,sizeof(int),1,afile);
              fread(&HGAP_MIN,sizeof(int),1,afile);
            }
          else if (extra == 3*sizeof(int))
            { fread(&GOOD_QV,sizeof(int),1,afile);
              fread(&BAD_QV,sizeof(int),1,afile);
              fread(&COVERAGE,sizeof(int),1,afile);
              HGAP_MIN = 0;
            }
          else if (extra == 2*sizeof(int))
            { fread(&GOOD_QV,sizeof(int),1,afile);
              fread(&BAD_QV,sizeof(int),1,afile);
              COVERAGE = -1;
              HGAP_MIN = 0;
            }
          else
            { GOOD_QV  = -1;
              BAD_QV   = -1;
              COVERAGE = -1;
              HGAP_MIN = 0;
            }
          fclose(afile);
        }
      else
        { fprintf(stderr,"%s: Must have a 'trim' track, run DAStrim\n",Prog_Name);
          exit (1);
        }

      root   = Root(argv[c],".db");
      nreads = DB->nreads;
      totlen = DB->totlen;

      nelim   = 0;
      n5trm   = 0;
      n3trm   = 0;
      natrm   = 0;
      nelimbp = 0;
      n5trmbp = 0;
      n3trmbp = 0;
      natrmbp = 0;

      ngaps   = 0;
      nlowq   = 0;
      nspan   = 0;
      nchim   = 0;
      ngapsbp = 0;
      nlowqbp = 0;
      nspanbp = 0;
      nchimbp = 0;

      for (a = 0; a < DB->nreads; a++)
        { tb = TRIM_IDX[a];
          te = TRIM_IDX[a+1];
          alen  = DB->reads[a].rlen;
          if (alen < HGAP_MIN)
            { nreads -= 1;
              totlen -= alen;
            }
          else if (tb >= te)
            { nelim += 1;
              nelimbp += alen;
            }
          else
            { if (TRIM[tb] > 0)
                { n5trm += 1;
                  n5trmbp += TRIM[tb];
                } 
              if (TRIM[te-1] < alen)
                { n3trm += 1;
                  n3trmbp += alen - TRIM[te-1];
                } 
              while (tb + 3 < te)
                { ngaps += 1;
                  ngapsbp += TRIM[tb+3] - TRIM[tb+1];
                  if (TRIM[tb+2] == LOWQ)
                    { nlowq += 1;
                      nlowqbp += TRIM[tb+3] - TRIM[tb+1];
                    }
                  else if (TRIM[tb+2] == SPAN)
                    { nspan += 1;
                      nspanbp += TRIM[tb+3] - TRIM[tb+1];
                    }
                  else if (TRIM[tb+2] == ADPRE)
                    { natrm += 1;
                      natrmbp += TRIM[tb+3] - TRIM[tb];
                    }
                  else if (TRIM[tb+2] == ADSUF)
                    { natrm += 1;
                      natrmbp += TRIM[tb+4] - TRIM[tb+1];
                    }
                  else
                    { nchim += 1;
                      nchimbp += TRIM[tb+3] - TRIM[tb+1];
                    }
                  tb += 3;
                }
            }
        }

      printf("\nStatistics for DAStrim");
      if (GOOD_QV >= 0)
        printf(" -g%d",GOOD_QV);
      else
        printf(" -g??");
      if (BAD_QV >= 0)
        printf(" -b%d",BAD_QV);
      else
        printf(" -b??");
      if (HGAP_MIN > 0)
        printf(" [-H%d]",HGAP_MIN);
      printf(" %s\n\n",root);

      printf("  Input:    ");
      Print_Number((int64) nreads,7,stdout);
      printf(" (100.0%%) reads     ");
      Print_Number(totlen,12,stdout);
      printf(" (100.0%%) bases");
      if (HGAP_MIN > 0)
        { printf(" (another ");
          Print_Number((int64) (DB->nreads-nreads),0,stdout);
          printf(" were < H-length)");
        }
      printf("\n");

      printf("  Trimmed:  ");
      Print_Number(nelim,7,stdout);
      printf(" (%5.1f%%) reads     ",(100.*nelim)/nreads);
      Print_Number(nelimbp,12,stdout);
      printf(" (%5.1f%%) bases\n",(100.*nelimbp)/totlen);

      printf("  5' trim:  ");
      Print_Number(n5trm,7,stdout);
      printf(" (%5.1f%%) reads     ",(100.*n5trm)/nreads);
      Print_Number(n5trmbp,12,stdout);
      printf(" (%5.1f%%) bases\n",(100.*n5trmbp)/totlen);

      printf("  3' trim:  ");
      Print_Number(n3trm,7,stdout);
      printf(" (%5.1f%%) reads     ",(100.*n3trm)/nreads);
      Print_Number(n3trmbp,12,stdout);
      printf(" (%5.1f%%) bases\n",(100.*n3trmbp)/totlen);

      if (natrm > 0)
        { printf("  Adapter:  ");
          Print_Number(natrm,7,stdout);
          printf(" (%5.1f%%) reads     ",(100.*natrm)/nreads);
          Print_Number(natrmbp,12,stdout);
          printf(" (%5.1f%%) bases\n",(100.*natrmbp)/totlen);
        }

      printf("\n");

      printf("  Gaps:     ");
      Print_Number(ngaps,7,stdout);
      printf(" (%5.1f%%) gaps      ",(100.*(ngaps))/nreads);
      Print_Number(ngapsbp,12,stdout);
      printf(" (%5.1f%%) bases\n",(100.*(ngapsbp))/totlen);

      printf("    Low QV: ");
      Print_Number(nlowq,7,stdout);
      printf(" (%5.1f%%) gaps      ",(100.*(nlowq))/nreads);
      Print_Number(nlowqbp,12,stdout);
      printf(" (%5.1f%%) bases\n",(100.*(nlowqbp))/totlen);

      printf("    Span'd: ");
      Print_Number(nspan,7,stdout);
      printf(" (%5.1f%%) gaps      ",(100.*(nspan))/nreads);
      Print_Number(nspanbp,12,stdout);
      printf(" (%5.1f%%) bases\n",(100.*(nspanbp))/totlen);

      printf("    Break:  ");
      Print_Number(nchim,7,stdout);
      printf(" (%5.1f%%) gaps      ",(100.*(nchim))/nreads);
      Print_Number(nchimbp,12,stdout);
      printf(" (%5.1f%%) bases\n",(100.*(nchimbp))/totlen);

      printf("\n");

      printf("  Clipped:  ");
      Print_Number(n5trm+n3trm+nelim+nchim,7,stdout);
      printf(" clips              ");
      Print_Number(n5trmbp+n3trmbp+nelimbp+nchimbp,12,stdout);
      printf(" (%5.1f%%) bases\n",(100.*(n5trmbp+n3trmbp+nelimbp+nchimbp))/totlen);

      printf("  Patched:  ");
      Print_Number(nlowq+nspan,7,stdout);
      printf(" patches            ");
      Print_Number(nlowqbp+nspanbp,12,stdout);
      printf(" (%5.1f%%) bases\n",(100.*(nlowqbp+nspanbp))/totlen);

      free(root);
      Close_DB(DB);
    }

  free(Prog_Name);
  exit (0);
}
