/*******************************************************************************************
 *
 *  Using overlap pile for each read compute estimated intrinisic quality values
 *
 *  Author:  Gene Myers
 *  Date  :  September 2015
 *
 *******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "DB.h"
#include "align.h"

#ifdef HIDE_FILES
#define PATHSEP "/."
#else
#define PATHSEP "/"
#endif

#undef  QV_DEBUG

static char *Usage = "[-v] -c<int> <source:db> <overlaps:las> ...";

#define  MAXQV   50     //  Max QV score is 50
#define  MAXQV1  51
#define  MINCOV   2     //  To have a score must be covered >= MINCOV in each direction (must be >0)

#define  PARTIAL .20    //  Partial terminal segments covering this percentage are scored

static int     QV_DEEP; //  # of best diffs to average for QV score
static int     VERBOSE;

static int     TRACE_SPACING;  //  Trace spacing (from .las file)
static int     TBYTES;         //  Bytes per trace segment (from .las file)

static HITS_DB _DB, *DB  = &_DB;   //  Data base
static int     DB_FIRST;           //  First read of DB to process
static int     DB_LAST;            //  Last read of DB to process (+1)
static int     DB_PART;            //  0 if all, otherwise block #

static FILE   *QV_AFILE;   //  .qual.anno
static FILE   *QV_DFILE;   //  .qual.data
static int64   QV_INDEX;   //  Current index into .qual.data file

//  Statistics

static int64 nreads, totlen;
static int64 qgram[MAXQV1], sgram[MAXQV1];


//  For each pile, calculate QV scores of the aread at tick spacing TRACE_SPACING

static void CALCULATE_QVS(int aread, Overlap *ovls, int novl)
{ static int    nmax = 0;
  static int   *hist = NULL;
  static int   *cist = NULL;
  static uint8 *qvec = NULL;
  static int    partial;

  int  alen, atick;
  int *tick, *cick;
  int  i;

  alen  = DB->reads[aread].rlen;
  atick = (alen + (TRACE_SPACING-1))/TRACE_SPACING;

#if defined(QV_DEBUG)
  printf("AREAD %d",aread);
  if (novl == 0)
    printf(" EMPTY");
  printf("\n");
#endif

  //  QV SCORES
  //       Allocate or expand data structures for qv calculation as needed

  if (atick > nmax)
    { nmax = atick*1.2 + 100;

      hist = (int *) Realloc(hist,nmax*MAXQV1*sizeof(int),"Allocating histograms");
      cist = (int *) Realloc(cist,nmax*MAXQV1*sizeof(int),"Allocating histograms");
      qvec = (uint8 *) Realloc(qvec,nmax*sizeof(uint8),"Allocating QV vector");

      if (hist == NULL || cist == NULL || qvec == NULL)
        exit (1);

      for (i = MAXQV1*nmax-1; i >= 0; i--)
        hist[i] = cist[i] = 0;
      partial = PARTIAL*TRACE_SPACING;
    }

  //  For every segment, fill histogram of match diffs for every one of the
  //   atick intervals, building separate histograms, hist & cist, for forward
  //   and reverse B-hits

  for (i = 0; i < novl; i++)
    { Path   *path;
      uint16 *trace;
      int    *ht;
      int     tlen, abit;
      int     a, b, x;

      path  = &(ovls[i].path);
      trace = (uint16 *) path->trace;
      tlen  = path->tlen;

      if (COMP(ovls[i].flags))
        ht = cist;
      else
        ht = hist;

      b = 0;
      a = (path->abpos/TRACE_SPACING)*MAXQV1;

      abit = (path->abpos % TRACE_SPACING);
      if (abit != 0)
        { a += MAXQV1;
          b += 2;
        }

      abit = (path->aepos % TRACE_SPACING);
      if (abit != 0)
        tlen -= 2;

      while (b < tlen)
        { x = (int) ((200.*trace[b]) / (TRACE_SPACING + trace[b+1]));
          if (x > MAXQV)
            x = MAXQV;
          ht[a + x] += 1;
          a += MAXQV1;
          b += 2;
        }

      if (path->aepos == alen && abit >= partial)
        { x = (int) ((200.*trace[tlen]) / (abit + trace[tlen+1]));
          if (x > MAXQV)
            x = MAXQV;
          ht[a + x] += 1;
        }
    }

  //  For every segment, qv score is the maximum of the averages of the QV_DEEP lowest
  //    in the forward and reverse directions (if each is QV_DEEP), or the average
  //    of overlap scores (if between MINCOV and QV_DEEP-1), or MAXQV if no overlaps at all.
  //    Reset histogram for segment to zeros.

  tick = hist;
  cick = cist;
  for (i = 0; i < atick; i++)
    { int  v, y;
      int  qvn, qvc;
      int  cntn, cntc;
      int  sumn, sumc;

#ifdef QV_DEBUG
      { int min, max;

        printf("   [%5d,%5d]:",i*TRACE_SPACING,(i+1)*TRACE_SPACING);

        for (v = 0; v <= MAXQV; v++)
          if (tick[v] > 0)
            break;
        min = v;

        for (v = MAXQV; v >= 0; v--)
          if (tick[v] > 0)
            break;
        max = v;

        for (v = min; v <= max; v++)
          if (tick[v] == 1)
            printf(" %2d",v);
          else if (tick[v] > 1)
            printf(" %2d(%d)",v,tick[v]);

        printf("\n                :");

        for (v = 0; v <= MAXQV; v++)
          if (cick[v] > 0)
            break;
        min = v;

        for (v = MAXQV; v >= 0; v--)
          if (cick[v] > 0)
            break;
        max = v;

        for (v = min; v <= max; v++)
          if (cick[v] == 1)
            printf(" %2d",v);
          else if (cick[v] > 1)
            printf(" %2d(%d)",v,cick[v]);
       }
#endif

      if (VERBOSE)
        for (v = 0; v <= MAXQV; v++)
          sgram[v] += tick[v] + cick[v];
     
      cntn = sumn = 0;
      for (v = 0; v <= MAXQV; v++)
        { y = tick[v];
          tick[v] = 0;
          cntn += y;
          sumn += y*v;
          if (cntn >= QV_DEEP)
            { sumn -= (cntn-QV_DEEP)*v;
              cntn  = QV_DEEP;
              break;
            }
        }
      for (v++; v <= MAXQV; v++)
        tick[v] = 0;
     
      cntc = sumc = 0;
      for (v = 0; v <= MAXQV; v++)
        { y = cick[v];
          cick[v] = 0;
          cntc += y;
          sumc += y*v;
          if (cntc >= QV_DEEP)
            { sumc -= (cntc-QV_DEEP)*v;
              cntc  = QV_DEEP;
              break;
            }
        }
      for (v++; v <= MAXQV; v++)
        cick[v] = 0;

      if (cntn >= MINCOV)
        qvn = sumn/cntn;
      else
        qvn = MAXQV;

      if (cntc >= MINCOV)
        qvc = sumc/cntc;
      else
        qvc = MAXQV;

      if (qvn > qvc)
        qvec[i] = (uint8) qvn;
      else
        qvec[i] = (uint8) qvc;

      tick += MAXQV1;
      cick += MAXQV1;

#ifdef QV_DEBUG
      printf(" >> %2d %2d = %2d <<\n",qvn,qvc,qvec[i]);
#endif
    }

  //  Accumulate qv histogram (if VERBOSE) and append qv's to .qual file

  if (VERBOSE)
    { for (i = 0; i < atick; i++)
        qgram[qvec[i]] += 1;
      nreads += 1;
      totlen += alen;
    }

  fwrite(qvec,sizeof(uint8),atick,QV_DFILE);
  QV_INDEX += atick;
  fwrite(&QV_INDEX,sizeof(int64),1,QV_AFILE);
}

  //  Read in each successive pile and call ACTION on it.  Read in the traces only if
  //   "trace" is nonzero

static int make_a_pass(FILE *input, void (*ACTION)(int, Overlap *, int), int trace)
{ static Overlap *ovls = NULL;
  static int      omax = 500;
  static uint16  *paths = NULL;
  static int      pmax = 100000;

  int64 i, j, novl;
  int   n, a;
  int   pcur;
  int   max;

  if (ovls == NULL)
    { ovls = (Overlap *) Malloc(sizeof(Overlap)*omax,"Allocating overlap buffer");
      if (ovls == NULL)
        exit (1);
    }
  if (trace && paths == NULL)
    { paths = (uint16 *) Malloc(sizeof(uint16)*pmax,"Allocating path buffer");
      if (paths == NULL)
        exit (1);
    }

  rewind(input);
  fread(&novl,sizeof(int64),1,input);
  fread(&TRACE_SPACING,sizeof(int),1,input);
  if (TRACE_SPACING <= TRACE_XOVR)
    TBYTES = sizeof(uint8);
  else
    TBYTES = sizeof(uint16);

  Read_Overlap(input,ovls);
  if (trace)
    { if (ovls[0].path.tlen > pmax)
        { pmax  = 1.2*(ovls[0].path.tlen)+10000;
          paths = (uint16 *) Realloc(paths,sizeof(uint16)*pmax,"Expanding path buffer");
          if (paths == NULL) exit (1);
        }
      fread(paths,TBYTES,ovls[0].path.tlen,input);
      if (TBYTES == 1)
        { ovls[0].path.trace = paths;
          Decompress_TraceTo16(ovls);
        }
    }
  else
    fseek(input,TBYTES*ovls[0].path.tlen,SEEK_CUR);

  if (ovls[0].aread < DB_FIRST)
    { fprintf(stderr,"%s: .las file overlaps don't correspond to reads in block %d of DB\n",
                     Prog_Name,DB_PART);
      exit (1);
    }

  pcur = 0;
  n = max = 0;
  for (j = DB_FIRST; j < DB_LAST; j++)
    { ovls[0] = ovls[n];
      a = ovls[0].aread;
      if (a != j)
        n = 0;
      else
        { if (trace)
            memmove(paths,paths+pcur,sizeof(uint16)*ovls[0].path.tlen);
          n = 1;
          pcur = ovls[0].path.tlen;
          while (1)
            { if (Read_Overlap(input,ovls+n) != 0)
                { ovls[n].aread = INT32_MAX;
                  break;
                }
              if (trace)
                { if (pcur + ovls[n].path.tlen > pmax)
                    { pmax = 1.2*(pcur+ovls[n].path.tlen)+10000;
                      paths = (uint16 *) Realloc(paths,sizeof(uint16)*pmax,"Expanding path buffer");
                      if (paths == NULL) exit (1);
                    }
                  fread(paths+pcur,TBYTES,ovls[n].path.tlen,input);
                  if (TBYTES == 1)
                    { ovls[n].path.trace = paths+pcur;
                      Decompress_TraceTo16(ovls+n);
                    }
                }
              else
                fseek(input,TBYTES*ovls[n].path.tlen,SEEK_CUR);
              if (ovls[n].aread != a)
                break;
              pcur += ovls[n].path.tlen;
              n    += 1;
              if (n >= omax)
                { omax = 1.2*n + 100;
                  ovls = (Overlap *) Realloc(ovls,sizeof(Overlap)*omax,"Expanding overlap buffer");
                  if (ovls == NULL) exit (1);
                }
            }

          if (n >= max)
            max = n;
          pcur = 0;
          for (i = 0; i < n; i++)
            { ovls[i].path.trace = paths+pcur;
              pcur += ovls[i].path.tlen;
            }
        }
      ACTION(j,ovls,n);
    }

  return (max);
}

int main(int argc, char *argv[])
{ FILE  *input;
  char  *root, *dpwd;
  char  *las, *lpwd;
  int64  novl;
  int    c, COVERAGE;

  //  Process arguments

  { int  i, j, k;
    int  flags[128];
    char *eptr;

    ARG_INIT("DASqv")

    COVERAGE = -1;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("v")
            break;
          case 'c':
            ARG_POSITIVE(COVERAGE,"Voting depth")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];

    if (argc < 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }

    if (COVERAGE < 0)
      { fprintf(stderr,"%s: Must supply -c parameter\n",Prog_Name);
        exit (1);
      }
    else
      { if (COVERAGE >= 40)
          QV_DEEP = COVERAGE/8;
        else if (COVERAGE >= 20)
          QV_DEEP = 5;
        else if (COVERAGE >= 4)
          QV_DEEP = COVERAGE/4;
        else
          { fprintf(stderr,"%s: Average coverage is too low (< 4X), cannot infer qv's\n",Prog_Name);
            exit (1);
          }
      }
  }

  //  Open trimmed DB

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
    Trim_DB(DB);
  }

  //  Initialize statistics gathering

  if (VERBOSE)
    { int i;

      nreads = 0;
      totlen = 0;
      for (i = 0; i <= MAXQV; i++)
        qgram[i] = sgram[i] = 0;

      printf("\nDASqv -c%d %s",COVERAGE,argv[1]);
      for (i = 2; i < argc; i++)
        printf(" %s",argv[i]);
      printf("\n");
    }

  //  Determine if overlap block is being processed and if so get first and last read
  //    from .db file

  dpwd = PathTo(argv[1]);
  root = Root(argv[1],".db");

  for (c = 2; c < argc; c++)
    { las  = Root(argv[c],".las");

      { FILE *dbfile;
        char  buffer[2*MAX_NAME+100];
        char *p, *eptr;
        int   i, part, nfiles, nblocks, cutoff, all, oindx;
        int64 size;

        DB_PART  = 0;
        DB_FIRST = 0;
        DB_LAST  = DB->nreads;

        p = rindex(las,'.');
        if (p != NULL)
          { part = strtol(p+1,&eptr,10);
            if (*eptr == '\0' && eptr != p+1)
              { dbfile = Fopen(Catenate(dpwd,"/",root,".db"),"r");
                if (dbfile == NULL)
                  exit (1);
                if (fscanf(dbfile,DB_NFILE,&nfiles) != 1)
                  SYSTEM_ERROR
                for (i = 0; i < nfiles; i++)
                  if (fgets(buffer,2*MAX_NAME+100,dbfile) == NULL)
                    SYSTEM_ERROR
                if (fscanf(dbfile,DB_NBLOCK,&nblocks) != 1)
                  SYSTEM_ERROR
                if (fscanf(dbfile,DB_PARAMS,&size,&cutoff,&all) != 3)
                  SYSTEM_ERROR
                for (i = 1; i <= part; i++)
                  if (fscanf(dbfile,DB_BDATA,&oindx,&DB_FIRST) != 2)
                    SYSTEM_ERROR
                if (fscanf(dbfile,DB_BDATA,&oindx,&DB_LAST) != 2)
                  SYSTEM_ERROR
                fclose(dbfile);
                DB_PART = part;
                *p = '\0';
              }
          }
      }

      //   Set up preliminary trimming track

      if (DB_PART > 0)
        { QV_AFILE = Fopen(Catenate(dpwd,PATHSEP,root,
                                    Numbered_Suffix(".",DB_PART,".qual.anno")),"w");
          QV_DFILE = Fopen(Catenate(dpwd,PATHSEP,root,
                                    Numbered_Suffix(".",DB_PART,".qual.data")),"w");
        }
      else
        { QV_AFILE = Fopen(Catenate(dpwd,PATHSEP,root,".qual.anno"),"w");
          QV_DFILE = Fopen(Catenate(dpwd,PATHSEP,root,".qual.data"),"w");
        }
      if (QV_AFILE == NULL || QV_DFILE == NULL)
        exit (1);

      { int size, nreads;

        nreads = DB_LAST - DB_FIRST;
        size   = sizeof(int64);
        fwrite(&nreads,sizeof(int),1,QV_AFILE);
        fwrite(&size,sizeof(int),1,QV_AFILE);
        QV_INDEX = 0;
        fwrite(&QV_INDEX,sizeof(int64),1,QV_AFILE);
      }

      //  Open overlap file

      lpwd = PathTo(argv[c]);
      if (DB_PART > 0)
        input = Fopen(Catenate(lpwd,"/",las,Numbered_Suffix(".",DB_PART,".las")),"r");
      else
        input = Fopen(Catenate(lpwd,"/",las,".las"),"r");
      if (input == NULL)
        exit (1);

      free(lpwd);
      free(las);

      //  Get trace point spacing information

      fread(&novl,sizeof(int64),1,input);
      fread(&TRACE_SPACING,sizeof(int),1,input);

      //  Process each read pile

      make_a_pass(input,CALCULATE_QVS,1);

      fclose(QV_AFILE);
      fclose(QV_DFILE);
    }

  //  If verbose output statistics summary to stdout

  if (VERBOSE)
    { int   i;
      int64 ssum, qsum;
      int64 stotal, qtotal;

      printf("\nInput:  ");
      Print_Number(nreads,7,stdout);
      printf("reads,  ");
      Print_Number(totlen,12,stdout);
      printf(" bases\n");

      stotal = qtotal = 0;
      for (i = 0; i <= MAXQV; i++)
        { stotal += sgram[i];
          qtotal += qgram[i];
        }

      printf("\nHistogram of q-values (average %d best)\n",2*QV_DEEP);
      printf("\n                 Input                 QV\n");
      qsum = qgram[MAXQV];
      ssum = sgram[MAXQV];
      printf("\n    %2d:  %9lld  %5.1f%%    %9lld  %5.1f%%\n\n",
             MAXQV,sgram[MAXQV],(100.*ssum)/stotal,qgram[MAXQV],(100.*qsum)/qtotal);

      qtotal -= qsum;
      stotal -= ssum;
      ssum = qsum = 0;
      for (i = MAXQV-1; i >= 0; i--) 
        if (qgram[i] > 0)
          { ssum += sgram[i];
            qsum += qgram[i];
            printf("    %2d:  %9lld  %5.1f%%    %9lld  %5.1f%%\n",
                   i,sgram[i],(100.*ssum)/stotal,
                     qgram[i],(100.*qsum)/qtotal);
          }
    }

  //  Clean up

  free(dpwd);
  free(root);

  Close_DB(DB);
  free(Prog_Name);

  exit (0);
}
