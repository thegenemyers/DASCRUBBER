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

static char *Usage = "[-v] [-c<int>] <source:db> <overlaps:las> ...";

#define  MAXQV   50     //  Max QV score is 50
#define  MAXQV1  51
#define  MINCOV   2     //  To have a score must be covered >= MINCOV in each direction (must be >0)

#define  PARTIAL .20    //  Partial terminal segments covering this percentage are scored

static int     VERBOSE;
static int     COVERAGE;  //  Estimated coverage of genome
static int       QV_DEEP;   //  # of best diffs to average for QV score
static int     HGAP_MIN;  //  Under this length do not process for HGAP

static int     TRACE_SPACING;  //  Trace spacing (from .las file)
static int     TBYTES;         //  Bytes per trace segment (from .las file)

static DAZZ_DB _DB, *DB  = &_DB;   //  Data base
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

  if (alen < HGAP_MIN)
    { fwrite(&QV_INDEX,sizeof(int64),1,QV_AFILE);
      return;
    }

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
      printf(" >> %2d(%d) %2d(%d) = %2d <<\n",qvn,cntn,qvc,cntc,qvec[i]);
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

  if (Read_Overlap(input,ovls) != 0)
    ovls[0].aread = INT32_MAX;
  else if (trace)
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

  if (ovls[n].aread < INT32_MAX)
    { fprintf(stderr,"%s: .las file overlaps don't correspond to reads in block %d of DB\n",
                     Prog_Name,DB_PART);
      exit (1);
    }

  return (max);
}

int main(int argc, char *argv[])
{ char  *root, *dpwd;
  int64  novl;
  int    c;

  DAZZ_EXTRA ex_hgap, ex_covr;
  DAZZ_EXTRA ex_cest, ex_qvs, ex_dif;

  char *cest_name = "Coverage Estimate";
  char *qvs_name  = "Histogram of QVs";
  char *dif_name  = "Histogram of Tile Differences";
  int64 cover64;

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
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: Verbose mode, output statistics as proceed.\n");
        fprintf(stderr,"      -c: Use this as the average coverage (not DAScover estimate).\n");
        exit (1);
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

  //  Get .covr track information

  { FILE      *afile;
    char      *aname;
    int        extra, cmax;
    int64     *cgram;

    aname = Strdup(Catenate(DB->path,".","covr",".anno"),"Allocating anno file");
    if (aname == NULL)
      exit (1);
    afile  = fopen(aname,"r");
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

    HGAP_MIN = (int) ((int64 *) (ex_hgap.value))[0];
    cgram    = (int64 *) (ex_covr.value);
    cmax     = ex_covr.nelem - 1;

    if (COVERAGE < 0)
      { int i;

        i = 0;
        while (cgram[i+1] < cgram[i])
          i += 1;
        for (COVERAGE = i++; i < cmax; i++)
          if (cgram[COVERAGE] < cgram[i])
            COVERAGE = i;
      }

    if (COVERAGE >= 40)
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

  //  Setup extras

  ex_cest.vtype = DB_INT;      // Estimated coverage (same for every .las)
  ex_cest.nelem = 1;
  ex_cest.accum = DB_EXACT;
  ex_cest.name  = cest_name;
  cover64 = COVERAGE;
  ex_cest.value = &cover64;

  ex_qvs.vtype = DB_INT;      //  Histogram of MAXQV1 trace-point diff counts
  ex_qvs.nelem = MAXQV1;
  ex_qvs.accum = DB_SUM;
  ex_qvs.name  = qvs_name;
  ex_qvs.value = &qgram;

  ex_dif.vtype = DB_INT;     //  Histogram of MAXQV1 intrinisic qv counts
  ex_dif.nelem = MAXQV1;
  ex_dif.accum = DB_SUM;
  ex_dif.name  = dif_name;
  ex_dif.value = &sgram;

  //  For each .las file do

  dpwd = PathTo(argv[1]);
  root = Root(argv[1],".db");

  for (c = 2; c < argc; c++)
    { Block_Looper *parse;
      FILE         *input;

      parse = Parse_Block_Arg(argv[c]);

      while ((input = Next_Block_Arg(parse)) != NULL)
        { DB_PART  = 0;
          DB_FIRST = 0;
          DB_LAST  = DB->nreads;

          //  Determine if overlap block is being processed and if so get first and last read
          //    from .db file

          { FILE *dbfile;
            char  buffer[2*MAX_NAME+100];
            char *p, *eptr;
            int   i, part, nfiles, nblocks, cutoff, all, oindx;
            int64 size;

            p = rindex(Block_Arg_Root(parse),'.');
            if (p != NULL)
              { part = strtol(p+1,&eptr,10);
                if (*eptr == '\0' && eptr != p+1)
                  { dbfile = Fopen(Catenate(dpwd,"/",root,".db"),"r");
                    if (dbfile == NULL)
                      exit (1);
                    if (fscanf(dbfile,DB_NFILE,&nfiles) != 1)
                      SYSTEM_READ_ERROR
                    for (i = 0; i < nfiles; i++)
                      if (fgets(buffer,2*MAX_NAME+100,dbfile) == NULL)
                        SYSTEM_READ_ERROR
                    if (fscanf(dbfile,DB_NBLOCK,&nblocks) != 1)
                      SYSTEM_READ_ERROR
                    if (fscanf(dbfile,DB_PARAMS,&size,&cutoff,&all) != 3)
                      SYSTEM_READ_ERROR
                    for (i = 1; i <= part; i++)
                      if (fscanf(dbfile,DB_BDATA,&oindx,&DB_FIRST) != 2)
                        SYSTEM_READ_ERROR
                    if (fscanf(dbfile,DB_BDATA,&oindx,&DB_LAST) != 2)
                      SYSTEM_READ_ERROR
                    fclose(dbfile);
                    DB_PART = part;
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

          { int size, length;

            length = DB_LAST - DB_FIRST;
            size   = sizeof(int64);
            fwrite(&length,sizeof(int),1,QV_AFILE);
            fwrite(&size,sizeof(int),1,QV_AFILE);
            QV_INDEX = 0;
            fwrite(&QV_INDEX,sizeof(int64),1,QV_AFILE);
          }

          //  Get trace point spacing information

          fread(&novl,sizeof(int64),1,input);
          fread(&TRACE_SPACING,sizeof(int),1,input);

          //  Initialize statistics gathering

          { int i;

            nreads = 0;
            totlen = 0;
            for (i = 0; i <= MAXQV; i++)
              qgram[i] = sgram[i] = 0;
          }

          if (VERBOSE)
            { printf("\n\nDASqv");
              if (HGAP_MIN > 0)
                printf(" -H%d",HGAP_MIN);
              printf(" -c%d %s %s\n\n",COVERAGE,argv[1],argv[c]);
              fflush(stdout);
            }

          //  Process each read pile

          make_a_pass(input,CALCULATE_QVS,1);

          //  Write out extras and close .qual track

          Write_Extra(QV_AFILE,&ex_hgap);
          Write_Extra(QV_AFILE,&ex_cest);
          Write_Extra(QV_AFILE,&ex_qvs);
          Write_Extra(QV_AFILE,&ex_dif);

          fclose(QV_AFILE);
          fclose(QV_DFILE);
          fclose(input);

          //  If verbose output statistics summary to stdout

          if (VERBOSE)
            { int   i;
              int64 ssum, qsum;
              int64 stotal, qtotal;
              int   gval, bval;

              printf("\n  Input:  ");
              Print_Number(nreads,7,stdout);
              printf("reads,  ");
              Print_Number(totlen,12,stdout);
              printf(" bases");
              if (HGAP_MIN > 0)
                { printf(" (another ");
                  Print_Number((DB_LAST-DB_FIRST) - nreads,0,stdout);
		  printf(" were < H-length)");
		}
	      printf("\n");

	      stotal = qtotal = 0;
	      for (i = 0; i <= MAXQV; i++)
	        { stotal += sgram[i];
                  qtotal += qgram[i];
                }

              printf("\n  Histogram of q-values (average %d best)\n",2*QV_DEEP);
              printf("\n                   Input                 QV\n");
              qsum = qgram[MAXQV];
              ssum = sgram[MAXQV];
              printf("\n      %2d:  %9lld  %5.1f%%    %9lld  %5.1f%%\n\n",
                     MAXQV,sgram[MAXQV],(100.*ssum)/stotal,qgram[MAXQV],(100.*qsum)/qtotal);

              bval    = gval = -1;
              qtotal -= qsum;
              stotal -= ssum;
              ssum = qsum = 0;
              for (i = MAXQV-1; i >= 0; i--) 
                if (qgram[i] > 0)
                  { ssum += sgram[i];
                    qsum += qgram[i];
                    printf("      %2d:  %9lld  %5.1f%%    %9lld  %5.1f%%\n",
                           i,sgram[i],(100.*ssum)/stotal,
                             qgram[i],(100.*qsum)/qtotal);
                    if ((100.*qsum)/qtotal > 7. && bval < 0)
                      bval = i;
                    if ((100.*qsum)/qtotal > 20. && gval < 0)
                      gval = i;
                  }

              printf("\n    Recommend \'DAStrim -g%d -b%d'\n\n",gval,bval);
            }
        }

      Free_Block_Arg(parse);
    }

  //  Clean up

  free(dpwd);
  free(root);

  Close_DB(DB);
  free(Prog_Name);

  exit (0);
    }
