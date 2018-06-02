/*******************************************************************************************
 *
 *  Using overlap pile for each read compute estimated coverage of the underlying genome
 *    generating a .covr track containing a histogram of the coverage of every unmasked
 *    trace-point tile.
 *
 *  Author:  Gene Myers
 *  Date  :  January 2018
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

#undef COVER_DEBUG

static char *Usage = "[-v] [-H<int>] [-m<track>]+ <source:db> <overlaps:las> ...";

#define  MAX_COVER  1000

static int    VERBOSE;
static int    HGAP_MIN;    //  Under this length do not process for HGAP

typedef struct
  { char  *name;
    int64 *idx;
    int   *data;
  } Mask;

static int   NUM_MASK;
static Mask *MASKS;

static int64  Cov_Hist[MAX_COVER+1];   //  [0..MAX_COVER] counts of trace-interval coverages
static char  *Cov_Name  = "Coverage Histogram";
static char  *Hgap_Name = "Hgap threshold";

static int     TRACE_SPACING;  //  Trace spacing (from .las file)
static int     TBYTES;         //  Bytes per trace segment (from .las file)

static DAZZ_DB   _DB, *DB  = &_DB;    //  Data base
static int        DB_FIRST;           //  First read of DB to process
static int        DB_LAST;            //  Last read of DB to process (+1)
static int        DB_PART;            //  0 if all, otherwise block #
static DAZZ_READ *Reads;	      //  Data base reads array

static FILE   *CV_AFILE;   //  .covr.anno
static char   *CV_ANAME;   //  ".covr.anno"


//  Statistics

static int64 nreads, totlen;

//  For each pile, calculate QV scores of the aread at tick spacing TRACE_SPACING

static void HISTOGRAM_COVER(int aread, Overlap *ovls, int novl)
{ static int    nmax = 0;
  static int   *local;
  static int   *cover;
  static int   *maskd;

  int  DBreads;
  int  alen, atick;
  int  bread, cssr;
  int  alow, ahigh;
  int  t2, t1;
  int  i, j, a, e;

  alen  = DB->reads[aread].rlen;
  atick = (alen + (TRACE_SPACING-1))/TRACE_SPACING;

  if (alen < HGAP_MIN)
    return;

#if defined(COVER_DEBUG)
  printf("AREAD %d",aread);
  if (novl == 0)
    printf(" EMPTY");
  printf("\n");
#endif

  //  COVERAGE
  //       Allocate or expand data structures for cover calculation as needed

  if (nmax == 0)
    { nmax = (DB->maxlen + (TRACE_SPACING-1))/TRACE_SPACING;

      local = (int *) Malloc(nmax*sizeof(int),"Allocating bread cover");
      cover = (int *) Malloc(nmax*sizeof(int),"Allocating aread cover");
      maskd = (int *) Malloc(nmax*sizeof(int),"Allocating aread cover");

      if (local == NULL || cover == NULL || maskd == NULL)
        exit (1);

      for (i = nmax-1; i >= 0; i--)
        local[i] = cover[i] = maskd[i] = 0;
    }

  //  For every segment, fill histogram of match diffs for every one of the
  //   atick intervals, building separate histograms, hist & cist, for forward
  //   and reverse B-hits

  t2 = TRACE_SPACING/2;
  t1 = t2-1;

  DBreads = DB->nreads;
  for (j = aread+1; j < DBreads; j++)
    if ((Reads[j].flags & DB_CSS) == 0)
      break;
  ahigh = j;
  for (j = aread; j >= 0; j--)
    if ((Reads[j].flags & DB_CSS) == 0)
      break;
  alow = j;

  for (i = 0; i < novl; i = j)
    { bread = ovls[i].bread;

      if (alow <= bread && bread < ahigh)
        { j += 1;
          continue;
        }

      for (j = bread+1; j < DBreads; j++)
        if ((Reads[j].flags & DB_CSS) == 0)
          break;
      cssr = j;

      for (j = i; j < novl; j++)
        if (ovls[j].bread < cssr)
          { e = (ovls[j].path.aepos + t2) / TRACE_SPACING;
            for (a = (ovls[j].path.abpos + t1) / TRACE_SPACING; a < e; a++)
              local[a] = 1;
          }
        else
          break;

      for (j = i; j < novl; j++)
        if (ovls[j].bread < cssr)
          { e = (ovls[j].path.aepos + t2) / TRACE_SPACING;
            for (a = (ovls[j].path.abpos + t1) / TRACE_SPACING; a < e; a++)
              if (local[a] > 0)
                { local[a] = 0;
                  cover[a] += 1;
                }
          }
        else
          break;
    }

  for (i = 0; i < NUM_MASK; i++)
    { Mask *m = (Mask *) (MASKS+i);
      int   k, u, e;

#ifdef COVER_DEBUG
      printf("  %s:",m->name);
#endif
      for (k = m->idx[aread]; k < m->idx[aread+1]; k += 2)
        { e = m->data[k+1];
#ifdef COVER_DEBUG
          printf(" [%d..%d]",m->data[k],e);
#endif
          e = (e + (TRACE_SPACING-1))/TRACE_SPACING;
          for (u = m->data[k]/TRACE_SPACING; u < e; u++) 
            maskd[u] = 1;
        }
#ifdef COVER_DEBUG
      printf("\n");
#endif
    }

#ifdef COVER_DEBUG
  printf("Mask: ");
  for (a = 0; a < atick; a++)
    printf("%d",maskd[a]);
  printf("\n");

  printf("Cover: ");
  for (a = 0; a < atick; a++)
    printf(" %d",cover[a]);
  printf("\n");
#endif

  for (a = 0; a < atick; a++)
    { if (maskd[a] == 0)
        { e = cover[a];
          if (e >= MAX_COVER)
            Cov_Hist[MAX_COVER] += 1;
          else
            Cov_Hist[e] += 1;
        }
      cover[a] = 0;
    }

  for (i = 0; i < NUM_MASK; i++)
    { Mask *m = (Mask *) (MASKS+i);
      int   k, u, e;

      for (k = m->idx[aread]; k < m->idx[aread+1]; k += 2)
        { e = m->data[k+1];
          e = (e + (TRACE_SPACING-1))/TRACE_SPACING;
          for (u = m->data[k]/TRACE_SPACING; u < e; u++) 
            maskd[u] = 0;
        }
    }

  if (VERBOSE)
    { nreads += 1;
      totlen += alen;
    }
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
  int64  novl, hgap64;
  int    c;

  DAZZ_EXTRA ex_covr, ex_hgap;

  //  Process arguments

  { int  i, j, k;
    int  flags[128];
    char *eptr;
    int   mmax;

    ARG_INIT("DAScover")

    HGAP_MIN = 0;
    NUM_MASK = 0;

    mmax  = 10;
    MASKS = (Mask *) Malloc(mmax*sizeof(Mask),"Allocating mask array");
    if (MASKS == NULL)
      exit (1);

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("v")
            break;
          case 'm':
            if (NUM_MASK >= mmax)
              { mmax  = 1.2*NUM_MASK + 10;
                MASKS = (Mask *) Realloc(MASKS,mmax*sizeof(Mask),"Reallocating mask array");
                if (MASKS == NULL)
                  exit (1);
              }
            MASKS[NUM_MASK++].name = argv[i]+2;
            break;
          case 'H':
            ARG_POSITIVE(HGAP_MIN,"HGAP threshold (in bp.s)")
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
        fprintf(stderr,"      -H: HGAP minimum length threshold.\n");
        fprintf(stderr,"      -m: repeat masks, stats not collected over these intervals\n");
        exit (1);
      }
  }

  //  Open trimmed DB and any mask tracks

  { DAZZ_TRACK *track;
    int64      *anno;
    int         status, kind;
    int         i, j, k;

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
    Reads = DB->reads;

    k = 0;
    for (i = 0; i < NUM_MASK; i++)
      { status = Check_Track(DB,MASKS[i].name,&kind);
        switch (status)
        { case 0:
            fprintf(stderr,"%s: [WARNING] %s track is for the *un*trimmed DB?\n",
                           Prog_Name,MASKS[i].name);
            continue;
          case -1:
            fprintf(stderr,"%s: [WARNING] %s track size not correct for trimmed DB.\n",
                           Prog_Name,MASKS[i].name);
            continue;
          case -2:
            fprintf(stderr,"%s: [WARNING] -m%s option given but no track found.\n",
                           Prog_Name,MASKS[i].name);
            continue;
          default:
            if (kind != MASK_TRACK)
              { fprintf(stderr,"%s: [WARNING] %s track is not a mask track.\n",
                               Prog_Name,MASKS[i].name);
                continue;
              }
            break;
        }

        track = Load_Track(DB,MASKS[i].name);

        anno = (int64 *) (track->anno);
        for (j = 0; j <= DB->nreads; j++)
          anno[j] /= sizeof(int);

        MASKS[k].name = MASKS[i].name;
        MASKS[k].idx  = anno;
        MASKS[k].data = (int *) (track->data);
        k += 1;
      }
    NUM_MASK = k;
  }

  //  Determine if overlap block is being processed and if so get first and last read
  //    from .db file

  ex_covr.vtype = DB_INT;
  ex_covr.nelem = MAX_COVER+1;
  ex_covr.accum = DB_SUM;
  ex_covr.name  = Cov_Name;
  ex_covr.value = Cov_Hist; 

  ex_hgap.vtype = DB_INT;
  ex_hgap.nelem = 1;
  ex_hgap.accum = DB_EXACT;
  ex_hgap.name  = Hgap_Name;
  hgap64 = HGAP_MIN;
  ex_hgap.value = &hgap64;

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

          //  Set up cover extra's track

          if (DB_PART > 0)
           CV_ANAME = Strdup(Catenate(dpwd,PATHSEP,root,
                           Numbered_Suffix(".",DB_PART,".covr.anno")),"Allocating cover anno name");
          else
           CV_ANAME = Strdup(Catenate(dpwd,PATHSEP,root,".covr.anno"),"Allocating cover anno name");
          CV_AFILE = Fopen(CV_ANAME,"w");
          if (CV_ANAME == NULL || CV_AFILE == NULL)
            exit (1);

          { int size, length;

            length = 0;
            size   = 1;
            fwrite(&length,sizeof(int),1,CV_AFILE);
            fwrite(&size,sizeof(int),1,CV_AFILE);
          }

          //  Get trace point spacing information

          fread(&novl,sizeof(int64),1,input);
          fread(&TRACE_SPACING,sizeof(int),1,input);

          //  Initialize statistics gathering

          if (VERBOSE)
            { nreads = 0;
              totlen = 0;
              printf("\nDAScover %s %s\n",argv[1],argv[c]);
            }

          { int i;

            for (i = 0; i <= MAX_COVER; i++)
              Cov_Hist[i] = 0;
          }

          //  Process each read pile

          make_a_pass(input,HISTOGRAM_COVER,0);

          //  If verbose output statistics summary to stdout

          if (VERBOSE)
            { int   i, cover;
              int64 ssum, stotal;

              printf("\nInput:  ");
              Print_Number(nreads,7,stdout);
              printf(" reads,  ");
              Print_Number(totlen,12,stdout);
              printf(" bases");
              if (HGAP_MIN > 0)
                { printf(" (another ");
                  Print_Number((DB_LAST-DB_FIRST) - nreads,0,stdout);
                  printf(" were < H-length)");
                }
              printf("\n");
 
              stotal = 0;
              for (i = 0; i <= MAX_COVER; i++)
                stotal += Cov_Hist[i];

              printf("\nCoverage Histogram\n\n");
              ssum = Cov_Hist[MAX_COVER];
              if (ssum > 0)
                printf("    %4d:  %9lld  %5.1f%%\n\n",
                       MAX_COVER,Cov_Hist[MAX_COVER],(100.*ssum)/stotal);
              stotal -= ssum;
              ssum    = 0;
              for (i = MAX_COVER-1; i >= 0; i--) 
                if (Cov_Hist[i] > 0)
                  { ssum += Cov_Hist[i];
                    printf("    %4d:  %9lld  %5.1f%%\n",
                           i,Cov_Hist[i],(100.*ssum)/stotal);
                  }

              i = 0;
              while (Cov_Hist[i+1] < Cov_Hist[i])
                i += 1;
              for (cover = i++; i < MAX_COVER; i++)
                if (Cov_Hist[cover] < Cov_Hist[i])
                  cover = i;

              printf("\n  Coverage is estimated at %d\n\n",cover);
            }

          //  Output coverage histogram

          Write_Extra(CV_AFILE,&ex_covr);
          Write_Extra(CV_AFILE,&ex_hgap);

          fclose(CV_AFILE);
          free(CV_ANAME);
          fclose(input);
        }

      Free_Block_Arg(parse);
    }

  free(dpwd);
  free(root);

  Close_DB(DB);
  free(Prog_Name);

  exit (0);
}
