/************************************************************************************\
*                                                                                    *
* Copyright (c) 2015, Dr. Eugene W. Myers (EWM). All rights reserved.                *
*                                                                                    *
* Redistribution and use in source and binary forms, with or without modification,   *
* are permitted provided that the following conditions are met:                      *
*                                                                                    *
*  · Redistributions of source code must retain the above copyright notice, this     *
*    list of conditions and the following disclaimer.                                *
*                                                                                    *
*  · Redistributions in binary form must reproduce the above copyright notice, this  *
*    list of conditions and the following disclaimer in the documentation and/or     *
*    other materials provided with the distribution.                                 *
*                                                                                    *
*  · The name of EWM may not be used to endorse or promote products derived from     *
*    this software without specific prior written permission.                        *
*                                                                                    *
* THIS SOFTWARE IS PROVIDED BY EWM ”AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES,    *
* INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND       *
* FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL EWM BE LIABLE   *
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES *
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS  *
* OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY      *
* THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING     *
* NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN  *
* IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                                      *
*                                                                                    *
* For any issues regarding this software and its use, contact EWM at:                *
*                                                                                    *
*   Eugene W. Myers Jr.                                                              *
*   Bautzner Str. 122e                                                               *
*   01099 Dresden                                                                    *
*   GERMANY                                                                          *
*   Email: gene.myers@gmail.com                                                      *
*                                                                                    *
\************************************************************************************/

/*******************************************************************************************
 *
 *  Given the "trim" track patching directions, produce a new database <X>.E.db from
 *    <X>.db containing all the patched and cut images of the original reads.
 *
 *  Author:  Gene Myers
 *  Date  :  August 2016
 *
 *******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "DB.h"
#include "align.h"

#ifdef HIDE_FILES
#define PATHSEP "/."
#else
#define PATHSEP "/"
#endif

#undef   DEBUG_PATCHING
#undef    SHOW_PIECES
#undef   DEBUG_MAPPING
#undef   DEBUG_INDEX

static char *Usage = "[-v] [-x<int>] <source:db> <target:db>";

//  Gap states

#define LOWQ  0   //  Gap is spanned by many LAs and patchable
#define SPAN  1   //  Gap has many paired LAs and patchable
#define SPLIT 2   //  Gap is a chimer or an unpatchable gap

//  Global Variables (must exist across the processing of each pile)

  //  Read-only

static DAZZ_DB _DB, *DB = &_DB;    //  Data base

static int64   *TRIM_IDX;
static int     *TRIM;

static int64   *PATCH_IDX;
static int     *PATCH;


void Print_Seq(char *target, int tlen)
{ static int letter[4] = { 'a', 'c', 'g', 't' };
  int i;

  for (i = 0; i < tlen; i++)
    { if (i%100 == 0 && i != 0)
        printf("\n");
      printf("%c",letter[(int) target[i]]);
    }
  printf("\n");
}

#define STACK_SIZE 50

static char *BSTACK[STACK_SIZE];

static int Load_Model(int *patch, char *target, int depth)
{ int   bread, bbeg, bend;
  int   lend, rend;
  int   gb, ge, pb;
  int   tlen, plen;
  char *bseq;

  if (BSTACK[depth] == NULL)
    BSTACK[depth] = New_Read_Buffer(DB);

  bread = patch[0];
  if (bread < 0)
    bread = - (bread+1);
  else
    bread -= 1;
  bbeg   = patch[1];
  bend   = patch[2];

#ifdef DEBUG_PATCHING
  printf("%*sPatching %d%c[%d,%d]->[%d,%d]\n",
         2*depth,"",bread,(patch[0]<0?'c':'n'),patch[1],patch[2],bbeg,bend);
#endif

  bseq = Load_Subread(DB,bread,bbeg,bend,BSTACK[depth],0) - bbeg;

  gb = TRIM_IDX[bread];
  ge = TRIM_IDX[bread+1];
  pb = PATCH_IDX[bread];

  tlen = 0;
  rend = -1;
  for (; gb < ge; gb += 3)
    { lend = TRIM[gb];
      if (lend > bbeg)
        { if (rend < bbeg || lend > bend || PATCH[pb] == 0)
            return (-1);
          plen = Load_Model(PATCH+pb,target+tlen,depth+1);
          if (plen < 0)
            return (-1);
          tlen += plen;

          rend = TRIM[gb+1];
          if (rend > bend)
            rend = bend;
          memmove(target+tlen,bseq+lend,rend-lend);
          tlen += rend-lend; 
#ifdef DEBUG_PATCHING
          printf("%*s  Piece %d,%d\n",2*depth,"",lend,rend);
#endif
#ifdef SHOW_PIECES
          Print_Seq(bseq+lend,rend-lend);
#endif
          pb += 3;
        }
      else // lend <= bbeg
        { rend = TRIM[gb+1];
          if (rend > bbeg)
            { if (rend > bend)
                rend = bend;
              memmove(target+tlen,bseq+bbeg,rend-bbeg);
              tlen += rend-bbeg;
#ifdef DEBUG_PATCHING
              printf("%*s  Piece %d,%d\n",2*depth,"",bbeg,rend);
#endif
#ifdef SHOW_PIECES
              Print_Seq(bseq+bbeg,rend-bbeg);
#endif
            }
          else
            pb += 3;
	}
      if (rend >= bend)
        break;
    }
  if (gb >= ge)
    fprintf(stderr,"Fatal: Should not happen\n");
 
  if (patch[0] < 0)
    Complement_Seq(target,tlen);

  return (tlen);
}

int main(int argc, char *argv[])
{ DAZZ_READ  *reads;
  int         nreads;
  int64       boff;

  int         nnew, nmax;
  int64       ntot;

  int        *segfate;
  int         nsegs;

  char       *pwd1, *root1;  //  inputs
  char       *pwd2, *root2;
  int         VERBOSE;
  int         CUTOFF;
  int         HGAP_MIN;

  int         nfiles;        //  contents of source .db file
  int         nblocks;
  int64       bsize;
  char      **flist = NULL;
  char      **plist = NULL;
  int        *findx = NULL;
  int        *bindx = NULL;

  FILE       *NB_FILE;       //  files for writing target
  FILE       *IDX_FILE;
  FILE       *BPS_FILE;
  FILE       *MP_AFILE;
  FILE       *MP_DFILE;

  //  Process arguments

  { int   i, j, k;
    int   flags[128];
    char *eptr;

    ARG_INIT("DASedit")

    CUTOFF = 0;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("v")
            break;
          case 'x':
            ARG_NON_NEGATIVE(CUTOFF,"Min read length cutoff")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];

    if (argc != 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: Verbose mode, output statistics as proceed.\n");
        fprintf(stderr,"      -x: minimum length for edited reads\n");
        exit (1);
      }
  }

  //  Open source and make target has a different pwd/name

  { int status;

    status = Open_DB(argv[1],DB);
    if (status < 0)
      exit (1);
    if (status == 1)
      { fprintf(stderr,"%s: Cannot be called on a .dam index: %s\n",Prog_Name,argv[1]);
        exit (1);
      }
    if (DB->part)
      { fprintf(stderr,"%s: Cannot be called on a block: %s\n",Prog_Name,argv[1]);
        exit (1);
      }

    pwd1  = PathTo(argv[1]);
    root1 = Root(argv[1],".db");

    pwd2  = PathTo(argv[2]);
    root2 = Root(argv[2],".db");

    if (strcmp(root1,root2) == 0 && strcmp(pwd1,pwd2) == 0)
      { fprintf(stderr,"%s: source and target are the same !\n",Prog_Name);
        exit (1);
      }
  }

  //  Load the file and block structure of the source database

  { int   i;
    int   nid, oid;
    int   cutoff, allflag, ufirst;
    FILE *dstub;

    dstub = Fopen(Catenate(pwd1,"/",root1,".db"),"r");
    if (dstub == NULL)
      exit (1);

    if (fscanf(dstub,DB_NFILE,&nfiles) != 1)
      { fprintf(stderr,"%s: %s.db is corrupted, read failed\n",Prog_Name,root1);
        exit (1);
      }

    flist = (char **) Malloc(sizeof(char *)*nfiles,"Allocating file list");
    plist = (char **) Malloc(sizeof(char *)*nfiles,"Allocating prolog list");
    findx = (int *) Malloc(sizeof(int *)*nfiles,"Allocating file index");
    if (flist == NULL || plist == NULL || findx == NULL)
      exit (1);

    for (i = 0; i < nfiles; i++)
      { char prolog[MAX_NAME], fname[MAX_NAME];

        if (fscanf(dstub,DB_FDATA,findx+i,fname,prolog) != 3)
          { fprintf(stderr,"%s: %s.db is corrupted, read failed\n",Prog_Name,root1);
            exit (1);
          }
        if ((flist[i] = Strdup(fname,"Adding to file list")) == NULL)
          exit (1);
        if ((plist[i] = Strdup(prolog,"Adding to prolog list")) == NULL)
          exit (1);
      }

    if (fscanf(dstub,DB_NBLOCK,&nblocks) != 1)
      { fprintf(stderr,"%s: %s.db is corrupted, read failed\n",Prog_Name,root1);
        exit (1);
      }
    if (fscanf(dstub,DB_PARAMS,&bsize,&cutoff,&allflag) != 3)
      { fprintf(stderr,"%s: %s.db is corrupted, read failed\n",Prog_Name,root1);
        exit (1);
      }

    bindx = (int *) Malloc(sizeof(int *)*(nblocks+1),"Allocating block indices");
    if (bindx == NULL)
      exit (1);

    for (i = 0; i <= nblocks; i++)
      if (fscanf(dstub,DB_BDATA,&ufirst,bindx+i) != 2)
        { fprintf(stderr,"%s: %s.db is corrupted, read failed\n",Prog_Name,root1);
          exit (1);
        }

    fclose(dstub);

  // map read counts for each file into trimmed read counts for each file

    if (allflag)
      allflag = 0;
    else
      allflag = DB_BEST;
    reads = DB->reads;

    nid = 0;
    oid = 0;
    for (i = 0; i < nfiles; i++)
      { while (oid < findx[i])
          { if ((reads[oid].flags & DB_BEST) >= allflag && reads[oid].rlen >= cutoff)
              nid++;
            oid += 1;
          }
        findx[i] = nid;
      }
  }

  //  Can now trim source DB.  Also load .trim and .patch tracks

  Trim_DB(DB);

  reads  = DB->reads;
  nreads = DB->nreads;

  { DAZZ_TRACK *track;
    int         i;

    track = Load_Track(DB,"trim");
    if (track != NULL)
      { FILE *afile;
        char *aname;
        int   extra, tracklen, size;
        DAZZ_EXTRA  ex_hgap;

        TRIM_IDX = (int64 *) track->anno;
        TRIM     = (int *) track->data;
        for (i = 0; i <= nreads; i++)
          TRIM_IDX[i] /= sizeof(int);

        //  Get HGAP minimum from .trim extras

        aname = Strdup(Catenate(DB->path,".","trim",".anno"),"Allocating anno file");
        if (aname == NULL)
          exit (1);
        afile  = fopen(aname,"r");

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
        fclose(afile);

        HGAP_MIN = (int) ((int64 *) (ex_hgap.value))[0];
      }
    else
      { fprintf(stderr,"%s: Must have a 'trim' track, run DAStrim\n",Prog_Name);
        exit (1);
      }

    track = Load_Track(DB,"patch");
    if (track != NULL)
      { PATCH_IDX = (int64 *) track->anno;
        PATCH     = (int *) track->data;
        for (i = 0; i <= nreads; i++)
          PATCH_IDX[i] /= sizeof(int);
      }
    else
      { fprintf(stderr,"%s: Must have a 'patch' track, run DASfix\n",Prog_Name);
        exit (1);
      }
  }

  //  Allocate segment status: will have value of defined constants below or if > 0
  //    the the segment is followed by a patch of given size.
  //    Also open new db files here, to be certain they can be so manipulated before
  //    doing anything that would need to be reversed.

#define SHORT_LAST -4
#define GOOD_LAST  -3
#define TRIMMED    -2
#define SHORT      -1
#define GOOD        0

  nsegs = nreads + PATCH_IDX[nreads]/3;

  segfate = Malloc(sizeof(int)*nsegs,"Allocating block status vector");

  NB_FILE  = Fopen(Catenate(pwd2,"/",root2,".db"),"w");
  IDX_FILE = Fopen(Catenate(pwd2,PATHSEP,root2,".idx"),"w");
  BPS_FILE = Fopen(Catenate(pwd2,PATHSEP,root2,".bps"),"w");
  MP_AFILE = Fopen(Catenate(pwd2,PATHSEP,root2,".map.anno"),"w");
  MP_DFILE = Fopen(Catenate(pwd2,PATHSEP,root2,".map.data"),"w");
  if (NB_FILE == NULL || IDX_FILE == NULL || BPS_FILE == NULL ||
      MP_AFILE == NULL || MP_DFILE == NULL)
    exit (1);

  //  Patch all the reads creating a new compressed .bps file of said for the new
  //     DB.  Further tally the total bytes and number of cuts, and also produce
  //     the .map track.

  { int   i, ni, bi;
    int   gb, ge, pb;
    char *target, *aseq;
    int64 MP_INDEX;

    int64 ntrim, nshort;
    int64 nttot, nstot;
    int64 htrim, httot;

    ni = 0;
    fwrite(&ni,sizeof(int),1,MP_AFILE);
    ni = 8;
    fwrite(&ni,sizeof(int),1,MP_AFILE);
    MP_INDEX = 0;
    fwrite(&MP_INDEX,sizeof(int64),1,MP_AFILE);

    for (i = 0; i < STACK_SIZE; i++)
      BSTACK[i] = NULL;

    { int ml = DB->maxlen;
      DB->maxlen = 1.5*ml + 10000;
      target = New_Read_Buffer(DB);
      DB->maxlen = ml;
    }

    aseq   = BSTACK[0] = New_Read_Buffer(DB);
    segfate[0] = GOOD_LAST;

    ntrim = nshort = 0;
    nstot = nttot = 0;
    htrim = httot = 0;
    ntot  = 0;
    nnew  = nmax = 0;

    boff = 0;
    ni = 0;
    bi = 0;
    pb = 0;
    ge = 0;
    for (i = 0; i < nreads; i++)
      { int  tlen, clen;
        int  lend, rend, blen;
        int  gb1, bi1;

        gb = ge;
        ge = TRIM_IDX[i+1];

#ifdef DEBUG_PATCHING
        printf("Doing %d\n",i);
#endif

        if (reads[i].rlen < HGAP_MIN)
          { segfate[bi++] = TRIMMED;
            htrim += 1;
            httot += reads[i].rlen;
            continue;
          }

        if (ge <= gb)
          { segfate[bi++] = TRIMMED;
            ntrim += 1;
            nttot += reads[i].rlen;
            continue;
          }

        Load_Read(DB,i,aseq,0);

        nttot += TRIM[gb];

        gb1  = gb;
        bi1  = bi;
        tlen = 0;
        for ( ; gb < ge; gb += 3)
          { lend = TRIM[gb];
            rend = TRIM[gb+1];
            blen = rend - lend;

            memmove(target+tlen,aseq+lend,blen);
            tlen += blen;
#ifdef DEBUG_PATCHING
            printf("  Piece %d,%d (%d)\n",lend,rend,bi);
#endif
#ifdef SHOW_PIECES
            Print_Seq(aseq+lend,blen);
#endif

            if (gb+3 < ge)
              { if (PATCH[pb] != 0)
                  clen = Load_Model(PATCH+pb,target+tlen,1);
                else
                  clen = -1;
                pb += 3;

                if (clen >= 0)
                  { tlen += clen;
                    segfate[bi++] = clen;
                    continue;
                  }
              }

            if (tlen >= CUTOFF)
              { Compress_Read(tlen,target);
                clen = COMPRESSED_LEN(tlen);
                fwrite(target,1,clen,BPS_FILE);
                boff += clen;

#ifdef DEBUG_PATCHING
                printf("  Output %d(%d)\n",ni,tlen);
#endif
                ni += 1;

                fwrite(&i,sizeof(int),1,MP_DFILE);
                fwrite(&(reads[i].rlen),sizeof(int),1,MP_DFILE);
                MP_INDEX += 2*sizeof(int);
                while (gb1 < gb)
                  { fwrite(TRIM+gb1,sizeof(int),1,MP_DFILE);
                    fwrite(TRIM+(gb1+1),sizeof(int),1,MP_DFILE);
                    fwrite(segfate+bi1,sizeof(int),1,MP_DFILE);
                    MP_INDEX += 3*sizeof(int);
                    gb1 += 3;
                    bi1 += 1;
                  }
                fwrite(TRIM+gb1,sizeof(int),1,MP_DFILE);
                fwrite(TRIM+(gb1+1),sizeof(int),1,MP_DFILE);
                MP_INDEX += 2*sizeof(int);
                fwrite(&MP_INDEX,sizeof(int64),1,MP_AFILE);

                gb1 = gb+3;
                if (gb1 <= ge)
                  segfate[bi++] = GOOD;
                else
                  segfate[bi++] = GOOD_LAST;
                bi1 = bi;
                nnew += 1;
                ntot += tlen;
                if (tlen > nmax)
                  nmax = tlen;
              }
            else
              { gb1 = gb+3;
                if (gb1 <= ge)
                  segfate[bi++] = SHORT;
                else
                  segfate[bi++] = SHORT_LAST;
                bi1 = bi;
                nshort += 1;
                nstot  += tlen;
#ifdef DEBUG_PATCHING
                printf("  Remove %d(%d)\n",ni,tlen);
#endif
              }
            tlen = 0;
            if (gb1 <= ge)
              { nttot += TRIM[gb1]-rend;
#ifdef DEBUG_PATCHING
                printf("  Cutting %d,%d\n",rend,TRIM[gb1]);
#endif
              }
            else
              nttot += reads[i].rlen-rend;
          }
      }

    nsegs = bi;
    for (i = 0; i < STACK_SIZE; i++)
      if (BSTACK[i] == NULL)
        break;
      else
        free(BSTACK[i]-1);
    free(target-1);

    rewind(MP_AFILE);
    fwrite(&ni,sizeof(int),1,MP_AFILE);

    if (VERBOSE)
      { printf("\n  ");

        if (htrim > 0)
          { Print_Number(htrim,0,stdout);
            printf(" reads and ");
            Print_Number(httot,0,stdout);
            printf(" bases in reads < H-length (%d)\n\n  ",HGAP_MIN);
          }

        Print_Number(DB->nreads-htrim,0,stdout);
        printf(" reads and ");
        Print_Number(DB->totlen-httot,0,stdout);
        if (htrim > 0)
          printf(" bases in reads >= H-length in source DB \n\n  ");
        else
          printf(" bases in source DB\n\n  ");

        Print_Number(ntrim,0,stdout);
        printf(" reads and ");
        Print_Number(nttot,0,stdout);
        printf(" bases were trimmmed by scrubbing (DAStrim)\n\n  ");

        if (CUTOFF > 0)
          { Print_Number(nshort,0,stdout);
            printf(" edited reads < %d bases, totaling ",CUTOFF);
            Print_Number(nstot,0,stdout);
            printf(" bases were cut (-x option)\n\n  ");
          }

        printf("The edited DB has ");
        Print_Number((int64) nnew,0,stdout);
        printf(" edited reads and ");
        Print_Number(ntot,0,stdout);
        printf(" bases\n");
     }

    fclose(BPS_FILE);
    fclose(MP_AFILE);
    fclose(MP_DFILE);
  }


  //  Output the file structure for the new database, adjusting the number
  //    of reads in each file according to how reads are split, and also adjust
  //    block indices similarly.  Further compress all read records that are trimmed
  //    or produce no edited reads.

  { int i, s;
    int oi, ni;
    int nf, nb;

#ifdef DEBUG_MAPPING
    printf("\nMAPPING\n");
#endif
    nf = 0;
    nb = 1;
    oi = 0;
    ni = 0;
    for (i = 0; i < nsegs; i++)
      { s = segfate[i];
        if (s <= 0)
          { if (s == GOOD || s == GOOD_LAST)
              ni += 1;
            if (s <= TRIMMED)
              { oi += 1;
                if (oi == findx[nf])
                  { findx[nf++] = ni;
#ifdef DEBUG_MAPPING
                    printf(" %2d: %d->%d\n",nf-1,oi,ni);
#endif
                  }
                if (oi == bindx[nb])
                  { bindx[nb++] = ni;
#ifdef DEBUG_MAPPING
                    printf(" %2d: %d->%d\n",nb-1,oi,ni);
#endif
                  }
              }
          }
      }

#ifdef DEBUG_MAPPING
    printf("  Total reads = %d  Trimmed Reads = %d\n",ni,oi);
    if (nf != nfiles)
      printf("  File count not correct %d %d\n",nf,nfiles);
    if (nb != nblocks+1)
      printf("  Block count not correct %d %d\n",nb,nblocks+1);
    fflush(stdout);
#endif

    fprintf(NB_FILE,DB_NFILE,nfiles);

    for (i = 0; i < nfiles; i++)
      fprintf(NB_FILE,DB_FDATA,findx[i],flist[i],plist[i]);

    fprintf(NB_FILE,DB_NBLOCK,nblocks);
    fprintf(NB_FILE,DB_PARAMS,bsize,CUTOFF,1);

    for (i = 0; i <= nblocks; i++)
      fprintf(NB_FILE,DB_BDATA,bindx[i],bindx[i]);

    fclose(NB_FILE);
  }

  { int       i, s;
    int       bi, gb;
    int       tlen, first;
    DAZZ_DB   NB;
    DAZZ_READ newrec;
#ifdef DEBUG_INDEX
    int       ni;
#endif

  //  Write an index of the patched reads

    NB = *DB;
    NB.ureads = nnew;
    NB.treads = nnew;
    NB.cutoff = CUTOFF;
    NB.allarr = DB->allarr | DB_ALL;
    NB.totlen = ntot;
    NB.maxlen = nmax;

    fwrite(&NB,sizeof(DAZZ_DB),1,IDX_FILE);

#ifdef DEBUG_INDEX
    printf("\nINDEXING\n");
    ni = 0;
#endif

    boff = 0;
    i    = 0;
    gb   = 0;
    bi   = 0;
    while (bi < nsegs)
      { tlen  = 0;
        first = TRIM[gb];
        while ((s = segfate[bi++]) > 0)
          { tlen += (TRIM[gb+1] - TRIM[gb]) + s; 
#ifdef DEBUG_INDEX
            printf(" [%d,%d] %d",TRIM[gb],TRIM[gb+1],s);
#endif
            gb += 3;
          }
        if (s == GOOD || s == GOOD_LAST)
          { tlen += (TRIM[gb+1] - TRIM[gb]); 
#ifdef DEBUG_INDEX
            printf(" [%d,%d]",TRIM[gb],TRIM[gb+1]);
            printf("  GOOD %d: (%d,%d)\n",i,ni++,bi);
#endif
            newrec.origin = reads[i].origin;
            newrec.fpulse = reads[i].fpulse + first;
            newrec.rlen   = tlen;
            newrec.boff   = boff;
            newrec.coff   = i;
            newrec.flags  = (reads[i].flags & DB_QV) | DB_BEST;
            if (segfate[bi] > TRIMMED) 
              newrec.flags |= DB_CSS;
            fwrite(&newrec,sizeof(DAZZ_READ),1,IDX_FILE);
            boff += COMPRESSED_LEN(tlen);

            if (s == GOOD)
              gb += 3;
            else
              { gb += 2;
                i  += 1;
              }
          }
        else if (s == TRIMMED)
          {
#ifdef DEBUG_INDEX
            printf("  TRIM %d: (%d)\n",i,bi);
#endif
            i += 1;
          }
        else // s == SHORT || s == SHORT_LAST
          {
#ifdef DEBUG_INDEX
            printf(" [%d,%d]",TRIM[gb],TRIM[gb+1]);
            printf("  SHRT %d: (%d)\n",i,bi);
#endif
            if (s == SHORT)
              gb += 3;
            else
              { gb += 2;
                i  += 1;
              }
          }
      }

#ifdef DEBUG_INDEX
    printf("Finish %d %d %lld\n",ni,i,boff);
#endif

    fclose(IDX_FILE);
  }

  free(pwd2);
  free(root2);
  free(pwd1);
  free(root1);
  free(Prog_Name);

  exit (0);
}
