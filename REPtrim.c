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

static int64  *QV_IDX;     //  qual track index
static uint8  *QV;         //  qual track values

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
      int         bval, gval, alen;
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
        { QV_IDX = (int64 *) track->anno;
          QV     = (uint8 *) track->data;
        }
      else
        { fprintf(stderr,"%s: Must have a 'qual' track, run DASqv\n",Prog_Name);
          exit (1);
        }

      track = Load_Track(DB,"trim");
      if (track != NULL)
        { FILE *afile;

          TRIM_IDX = (int64 *) track->anno;
          TRIM     = (int *) track->data;
          for (i = 0; i <= DB->nreads; i++)
            TRIM_IDX[i] /= sizeof(int);

          afile = fopen(Catenate(DB->path,".","trim",".anno"),"r");
          fseeko(afile,-2*sizeof(int),SEEK_END);
          fread(&gval,sizeof(int),1,afile);
          fread(&bval,sizeof(int),1,afile);
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

      for (a = 0; a < nreads; a++)
        { tb = TRIM_IDX[a];
          te = TRIM_IDX[a+1];
          alen  = DB->reads[a].rlen;
          if (tb >= te)
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

      printf("\nStatistics for DAStrim -g%d -b%d %s\n\n",gval,bval,root);

      printf("  Input:    ");
      Print_Number((int64) nreads,7,stdout);
      printf(" (100.0%%) reads     ");
      Print_Number(totlen,12,stdout);
      printf(" (100.0%%) bases\n");

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
