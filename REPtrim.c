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

    { DAZZ_DB    _DB, *DB  = &_DB;
      DAZZ_EXTRA ex_hgap, ex_cest, ex_good, ex_bad, ex_trim;

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

      //  Get .trim track extras

      { FILE *afile;
        char *aname;
        int   extra, tracklen, size;

        afile = NULL;
        if (DB->part)
          { aname = Strdup(Catenate(DB->path,Numbered_Suffix(".",DB->part,"."),"trim",".anno"),
                           "Allocating anno file");
            if (aname == NULL)
              exit (1);
            afile  = fopen(aname,"r");
            if (afile == NULL)
              { fprintf(stderr,"%s: Must have a 'trim.%d' track, run DAStrim\n",Prog_Name,DB->part);
                exit (1);
              }
          }
        else
          { aname = Strdup(Catenate(DB->path,".","trim",".anno"),"Allocating anno file");
            if (aname == NULL)
              exit (1);
            afile = fopen(aname,"r");
            if (afile == NULL)
              { fprintf(stderr,"%s: Must have a 'trim' track, run DAStrim\n",Prog_Name);
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
          { fprintf(stderr,"%s: Hgap threshold extra missing from .trim track?\n",Prog_Name);
            exit (1);
          }
        ex_cest.nelem = 0;
        if (Read_Extra(afile,aname,&ex_cest) != 0)
          { fprintf(stderr,"%s: Coverage estimate extra missing from .trim track?\n",Prog_Name);
            exit (1);
          }
        ex_good.nelem = 0;
        if (Read_Extra(afile,aname,&ex_good) != 0)
          { fprintf(stderr,"%s: Good QV threshold extra missing from .trim track?\n",Prog_Name);
            exit (1);
          }
        ex_bad.nelem = 0;
        if (Read_Extra(afile,aname,&ex_bad) != 0)
          { fprintf(stderr,"%s: Bad QV threshdold extra missing from .trim track?\n",Prog_Name);
            exit (1);
          }
        ex_trim.nelem = 0;
        if (Read_Extra(afile,aname,&ex_trim) != 0)
          { fprintf(stderr,"%s: Trimming statistics extra missing from .trim track?\n",Prog_Name);
            exit (1);
          }
        fclose(afile);
      }

      //  Generate Display

      { char       *root;
        int64       nreads, totlen;
        int64       nelim, nelimbp;
        int64       n5trm, n5trmbp;
        int64       n3trm, n3trmbp;
        int64       natrm, natrmbp;
        int64       ngaps, ngapsbp;
        int64       nlowq, nlowqbp;
        int64       nspan, nspanbp;
        int64       nchim, nchimbp;
        int         rlog, blog;
        int         cover, hgap_min;
        int         bad_qv, good_qv;
        int64      *tstats;

        //  Get relevant variables

        root   = Root(argv[c],".db");
        nreads = DB->nreads;
        totlen = DB->totlen;

        hgap_min = (int) ((int64 *) (ex_hgap.value))[0];
        cover    = (int) ((int64 *) (ex_cest.value))[0];
        good_qv  = (int) ((int64 *) (ex_good.value))[0];
        bad_qv   = (int) ((int64 *) (ex_bad.value))[0];
        tstats   = (int64 *) (ex_trim.value);

        nelim   = tstats[0];
        n5trm   = tstats[1];
        n3trm   = tstats[2];
        natrm   = tstats[3];
        nelimbp = tstats[4];
        n5trmbp = tstats[5];
        n3trmbp = tstats[6];
        natrmbp = tstats[7];

        ngaps   = tstats[8];
        nlowq   = tstats[9];
        nspan   = tstats[10];
        nchim   = tstats[11];
        ngapsbp = tstats[12];
        nlowqbp = tstats[13];
        nspanbp = tstats[14];
        nchimbp = tstats[15];

        printf("\nDAStrim");
        if (hgap_min > 0)
          printf(" [-H%d]",hgap_min);
        printf(" -c%d -g%d -b%d %s\n\n",cover,good_qv,bad_qv,root);

        //  Compensate for HGAP

        if (hgap_min > 0)
          { int i;

            for (i = 0; i < DB->nreads; i++)
              if (DB->reads[i].rlen < hgap_min)
                { nreads -= 1;
                  totlen -= DB->reads[i].rlen;
                }
          }

        //  Compute maximum field widths of statistics

        { int64 mult;

          rlog = 0;
          mult = 1;
          while (mult <= nreads || mult <= ngaps)
            { mult *= 10;
              rlog += 1;
            }
          if (rlog <= 3)
            rlog = 3;
          else
            rlog += (rlog-1)/3;
  
          blog = 0;
          mult = 1;
          while (mult <= totlen)
            { mult *= 10;
              blog += 1;
            }
          if (blog <= 3)
            blog = 3;
          else
            blog += (blog-1)/3;
        }

        //  Display the statistices

        printf("  Input:    ");
        Print_Number((int64) nreads,rlog,stdout);
        printf(" (100.0%%) reads     ");
        Print_Number(totlen,blog,stdout);
        printf(" (100.0%%) bases");
        if (hgap_min > 0)
          { printf(" (another ");
            Print_Number((int64) (DB->nreads-nreads),0,stdout);
            printf(" were < H-length)");
          }
        printf("\n");

        printf("  Trimmed:  ");
        Print_Number(nelim,rlog,stdout);
        printf(" (%5.1f%%) reads     ",(100.*nelim)/nreads);
        Print_Number(nelimbp,blog,stdout);
        printf(" (%5.1f%%) bases\n",(100.*nelimbp)/totlen);

        printf("  5' trim:  ");
        Print_Number(n5trm,rlog,stdout);
        printf(" (%5.1f%%) reads     ",(100.*n5trm)/nreads);
        Print_Number(n5trmbp,blog,stdout);
        printf(" (%5.1f%%) bases\n",(100.*n5trmbp)/totlen);

        printf("  3' trim:  ");
        Print_Number(n3trm,rlog,stdout);
        printf(" (%5.1f%%) reads     ",(100.*n3trm)/nreads);
        Print_Number(n3trmbp,blog,stdout);
        printf(" (%5.1f%%) bases\n",(100.*n3trmbp)/totlen);

        printf("  Adapter:  ");
        Print_Number(natrm,rlog,stdout);
        printf(" (%5.1f%%) reads     ",(100.*natrm)/nreads);
        Print_Number(natrmbp,blog,stdout);
        printf(" (%5.1f%%) bases\n",(100.*natrmbp)/totlen);

        printf("\n");

        printf("  Gaps:     ");
        Print_Number(ngaps,rlog,stdout);
        printf(" (%5.1f%%) gaps      ",(100.*(ngaps))/nreads);
        Print_Number(ngapsbp,blog,stdout);
        printf(" (%5.1f%%) bases\n",(100.*(ngapsbp))/totlen);

        printf("    Low QV: ");
        Print_Number(nlowq,rlog,stdout);
        printf(" (%5.1f%%) gaps      ",(100.*(nlowq))/nreads);
        Print_Number(nlowqbp,blog,stdout);
        printf(" (%5.1f%%) bases\n",(100.*(nlowqbp))/totlen);

        printf("    Span'd: ");
        Print_Number(nspan,rlog,stdout);
        printf(" (%5.1f%%) gaps      ",(100.*(nspan))/nreads);
        Print_Number(nspanbp,blog,stdout);
        printf(" (%5.1f%%) bases\n",(100.*(nspanbp))/totlen);

        printf("    Break:  ");
        Print_Number(nchim,rlog,stdout);
        printf(" (%5.1f%%) gaps      ",(100.*(nchim))/nreads);
        Print_Number(nchimbp,blog,stdout);
        printf(" (%5.1f%%) bases\n",(100.*(nchimbp))/totlen);

        printf("\n");

        printf("  Clipped:  ");
        Print_Number(n5trm+n3trm+nelim+nchim,rlog,stdout);
        printf(" clips              ");
        Print_Number(n5trmbp+n3trmbp+nelimbp+nchimbp,blog,stdout);
        printf(" (%5.1f%%) bases\n",(100.*(n5trmbp+n3trmbp+nelimbp+nchimbp))/totlen);

        printf("  Patched:  ");
        Print_Number(nlowq+nspan,rlog,stdout);
        printf(" patches            ");
        Print_Number(nlowqbp+nspanbp,blog,stdout);
        printf(" (%5.1f%%) bases\n",(100.*(nlowqbp+nspanbp))/totlen);

        free(root);
        Close_DB(DB);
      }
    }

  free(Prog_Name);
  exit (0);
}
