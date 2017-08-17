/*******************************************************************************************
 *
 *  Using overlap pile for each read,intrinisic quality values, and trimmed hq-intervals
 *    for each read, determine the B-read and segment thereof to use to patch each low
 *    quality segment between hq-intervals.
 *
 *  Author:  Gene Myers
 *  Date  :  June 2016
 *
 *******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include "DB.h"
#include "align.h"

#ifdef HIDE_FILES
#define PATHSEP "/."
#else
#define PATHSEP "/"
#endif

#undef   DEBUG_GAP_FILL
#undef     SHOW_PAIRS
#undef   DEBUG_SUMMARY

//  Command format and global parameter variables

static char *Usage = " [-v] <source:db> <overlaps:las> ...";

static int     BAD_QV;     //  qv >= and you are "bad"
static int     GOOD_QV;    //  qv <= and you are "good"
static int     VERBOSE;

//  Gap states

#define LOWQ  0   //  Gap is spanned by many LAs and patchable
#define SPAN  1   //  Gap has many paired LAs and patchable
#define SPLIT 2   //  Gap is a chimer or an unpatchable gap
#define ADPRE 3   //  Gap is an adatper break and prefix should be removed
#define ADSUF 4   //  Gap is an adatper break and suffix should be removed
#define NOPAT 3   //  Gap could not be patched (internal only)

#define  COVER_LEN     400  //  An overlap covers a point if it extends COVER_LEN to either side.
#define  ANCHOR_MATCH  .25  //  Delta in trace interval at both ends of patch must be < this %.

static int ANCHOR_THRESH;


//  Global Variables (must exist across the processing of each pile)

  //  Read-only

static int     TRACE_SPACING;  //  Trace spacing (from .las file)

static HITS_DB _DB, *DB  = &_DB;   //  Data base
static int     DB_FIRST;           //  First read of DB to process
static int     DB_LAST;            //  Last read of DB to process (+1)
static int     DB_PART;            //  0 if all, otherwise block #

static int64  *QV_IDX;     //  qual track index
static uint8  *QV;         //  qual track values

static int64  *TRIM_IDX;   //  trim track index
static int    *TRIM;       //  trim track values

  //  Read & Write

static FILE   *PR_AFILE;   //  .trim.anno
static FILE   *PR_DFILE;   //  .trim.data
static int64   PR_INDEX;   //  Current index into .trim.data file as it is being written

  // Statistics

static int fpatch, npatch;

//  Data Structures

typedef struct   //  General read interval [beg..end]
  { int beg;
    int end;
  } Interval;


/*******************************************************************************************
 *
 *  FIND ANY UNREMOVED ADAPTER (OR POLYMERASE SWITCHES) AND TRIM SMALLER PARTS
 *
 ********************************************************************************************/

typedef struct
  { int bread;   //  bread^comp[beg..end] is the patch sequence
    int comp;
    int beg;
    int end;

    int anc;     //  maximum anchor interval match
    int bad;     //  number of segments that are bad
    int avg;     //  average QV of the patch
  } Patch;

//  Evaluate the quality of lov->bread = rov->bread spaning [lcv,rcv] as a patch

static Patch *compute_patch(int lft, int rgt, Overlap *lov, Overlap *rov)
{ static Patch ans;

  uint16  *tr;
  int      bread, bcomp, blen;
  int      bb, be;
  int      t, te;
  int      bl, br;
  uint8   *qb;
  int      avg, anc, bad;

  if (lft > lov->path.aepos || rgt < rov->path.abpos)   // Cannot anchor
    return (NULL);
  if (lov->path.abpos > lft-TRACE_SPACING || rgt+TRACE_SPACING > rov->path.aepos)
    return (NULL);

  //  Get max of left and right anchors as anchor score

  tr = (uint16 *) lov->path.trace;
  te = 2 * (((lft + (TRACE_SPACING-1)) - lov->path.abpos)/TRACE_SPACING);
  if (te == 0)
    return (NULL);
  anc = tr[te-2];

  bb = lov->path.bbpos;
  for (t = 1; t < te; t += 2)
    bb += tr[t];

  tr = (uint16 *) rov->path.trace;
  te = 2 * (((rgt + (TRACE_SPACING-1)) - rov->path.abpos)/TRACE_SPACING);
  if (te >= rov->path.tlen)
    return (NULL);
  if (tr[te] > anc)
    anc = tr[te];

  be = rov->path.bepos;
  for (t = rov->path.tlen-1; t > te; t -= 2)
    be -= tr[t];

  if (bb >= be)
    return (NULL);

  bread = lov->bread;
  bcomp = COMP(lov->flags);

  //  Compute metrics for b-read patch

  if (bcomp)
    { blen = DB->reads[bread].rlen;
      t  = blen - be;
      be = blen - bb;
      bb = t;
    }

  bl = bb/TRACE_SPACING;
  br = (be+(TRACE_SPACING-1))/TRACE_SPACING;
  qb = QV + QV_IDX[bread];
  if (bl >= br)
    { avg = qb[bl];
      if (avg >= BAD_QV)
        bad = 1;
      else
        bad = 0;
    }
  else
    { avg = 0;
      bad = 0;
      for (t = bl; t < br; t++)
        { avg += qb[t];
          if (qb[t] >= BAD_QV)
            bad += 1;
        }
      avg /= (br-bl);
    }

  ans.bread = bread;
  ans.comp  = bcomp;
  ans.beg   = bb;
  ans.end   = be;
  ans.anc   = anc;
  ans.bad   = bad;
  ans.avg   = avg;

  return (&ans);
}

static int unsuitable(int bread, int lft, int rgt)
{ int tb, te;

  tb = TRIM_IDX[bread];
  te = TRIM_IDX[bread+1];
  for ( ; tb < te; tb += 3)
    if (TRIM[tb+1] >= lft)
      break;
  if (tb >= te || TRIM[tb] > lft)
    return (1);
  for ( ; tb < te ; tb += 3)
    { if (TRIM[tb+1] >= rgt)
        break;
      if (TRIM[tb+2] == SPLIT)
        return (1);
    }
  if (tb >= te || TRIM[tb] > rgt)
    return (1);
  return (0);
}

//  Categorize each gap and if appropriate return the best patch for each

static Patch *lowq_patch(Overlap *ovls, int novl, Interval *lblock, Interval *rblock)
{ static Patch patch;

  int j;
  int lft, rgt; 
  int lcv, rcv; 

  lft = lblock->end;
  rgt = rblock->beg;
  lcv = lft - COVER_LEN;
  rcv = rgt + COVER_LEN;
  if (lcv < lblock->beg)
    lcv = lblock->beg;
  if (rcv > rblock->end)
    rcv = rblock->end;

  patch.bread = -1;
  patch.anc   = TRACE_SPACING;
  patch.avg   = 100;
  for (j = 0; j < novl; j++)
    if (ovls[j].path.abpos <= lcv && ovls[j].path.aepos >= rcv)
      { Patch *can;

        can = compute_patch(lft,rgt,ovls+j,ovls+j);
            
        if (can == NULL) continue;

        if (unsuitable(can->bread,can->beg,can->end))
          continue;

        if (can->anc <= ANCHOR_THRESH && can->avg <= GOOD_QV && can->bad == 0 &&
            can->avg + can->anc < patch.anc + patch.avg)
          patch = *can;
      }

#ifdef DEBUG_GAP_FILL
  if (patch.bread >= 0)
    printf("    LOWQ  PATCH = %d%c[%d..%d]  %d (%d)\n",
           patch.bread,patch.comp?'c':'n',patch.beg,patch.end,patch.anc,patch.avg);
  else
    printf("    LOWQ PATCH FAIL\n");
#endif

  return (&patch);
}

static Patch *span_patch(Overlap *ovls, int novl, Interval *lblock, Interval *rblock)
{ static Patch patch;

  int j, k;
  int lft, rgt; 
  int lcv, rcv; 
  int bread, bcomp, blen;
  int ab, ae;
  int lidx, ridx, sidx, cidx;
  Patch *can;

  lft = lblock->end;
  rgt = rblock->beg;
  lcv = lft - COVER_LEN;
  rcv = rgt + COVER_LEN;
  if (lcv < lblock->beg)
    lcv = lblock->beg;
  if (rcv > rblock->end)
    rcv = rblock->end;

  //  Find LA pairs or LAs spanning the gap flank [lcv,rcv]

  patch.bread = -1;
  patch.bad   = DB->maxlen;
  patch.avg   = 100;
  for (j = 0; j < novl; j = k)
    { bread = ovls[j].bread;
      blen  = DB->reads[bread].rlen;
      bcomp = COMP(ovls[j].flags);
      if (bcomp)
        cidx = j;

      lidx = ridx = sidx = -1;    //  For all LA's with same b-read
      for (k = j; k < novl; k++)
        { if (ovls[k].bread != bread)
            break;
          if (COMP(ovls[k].flags) != (uint32) bcomp)   //  Note when b switches orientation
            { cidx  = k;
              bcomp = COMP(ovls[k].flags);
            }
          ab = ovls[k].path.abpos;
          ae = ovls[k].path.aepos;

#ifdef SHOW_PAIRS
          printf("\n %5d [%5d,%5d] %c",bread,ab,ae,COMP(ovls[k].flags)?'c':'n');
          if (ab <= lcv && ae >= rcv)
            printf("s");
          else
            printf(" ");
#endif

          //  Is LA a spanner, left-partner, or right partner

          if (ab <= lcv && ae >= rcv)
            { sidx = k;
              lidx = ridx = -1;
              continue;
            }

#ifdef SHOW_PAIRS
         if (ae >= rcv && ab <= rcv && ab - ovls[k].path.bbpos <= lft - COVER_LEN)
           printf("r");
         else
           printf(" ");
         if (ab <= lcv && ae >= lcv && ae + (blen-ovls[j].path.bepos) >= rgt + COVER_LEN)
           printf("l");
         else
           printf(" ");
#endif

          if (ae >= rcv && ab <= rcv && ab - ovls[k].path.bbpos <= lft - COVER_LEN)
            ridx = k;
          if (ab <= lcv && ae >= lcv && ae + (blen-ovls[j].path.bepos) >= rgt + COVER_LEN)
            lidx = k;
        }
      if (! bcomp)
        cidx = k;

#ifdef SHOW_PAIRS
      printf(" =");
      if (sidx >= 0)
        printf(" S");
      if (lidx >= 0)
        printf(" L");
      if (ridx >= 0)
        printf(" R");
      if (0 <= lidx && lidx < ridx && (ridx < cidx || lidx >= cidx))
        printf(" G");
      if ((0<=ridx && ridx<cidx && cidx<=lidx) || (0<=lidx && lidx<cidx && cidx<=ridx))
        printf(" A");
#endif

      //  Check for spanning pair

      if (sidx >= 0)
        lidx = ridx = sidx;
      else if (0 > lidx || lidx >= ridx || (ridx >= cidx && cidx > lidx))
        continue;

      //  Otherwise consider the gap linkable and try to patch it, declaring a split
      //    iff all patch attemtps fail

#ifdef DEBUG_GAP_FILL
      if (lidx != ridx)
        printf("       %5d [%5d,%5d] [%5d,%5d]",
               ovls[lidx].bread,ovls[lidx].path.abpos,ovls[lidx].path.aepos,
               ovls[ridx].path.abpos,ovls[ridx].path.aepos);
      else
        printf("       %5d [%5d,%5d] SSS",
               ovls[lidx].bread,ovls[lidx].path.abpos,ovls[lidx].path.aepos);
#endif

      can = compute_patch(lft,rgt,ovls+lidx,ovls+ridx);

      if (can != NULL)
        {
#ifdef DEBUG_GAP_FILL
          printf(" %d",can->end - can->beg);
#endif
          if ( ! unsuitable(can->bread,can->beg,can->end) && can->anc <= ANCHOR_THRESH)
            { if (can->bad < patch.bad)
                patch = *can;
              else if (can->bad == patch.bad)
                { if (can->avg < patch.avg)
                    patch = *can;
                }
#ifdef DEBUG_GAP_FILL
              printf("  AA %d %d(%d)",can->anc,can->bad,can->avg);
#endif
            }
        }
#ifdef DEBUG_GAP_FILL
      printf("\n");
#endif
    } 

#ifdef DEBUG_GAP_FILL
  if (patch.bread >= 0)
    printf("    SPAN %5d:  PATCH = %d%c[%d..%d]  %d %d(%d)\n",rgt-lft,
           patch.bread,patch.comp?'c':'n',patch.beg,
           patch.end,patch.anc,patch.bad,patch.avg);
  else
    printf("    SPAN PATCH FAIL\n");
#endif

  return (&patch);
}

/*******************************************************************************************
 *
 *  SCRUB EACH PILE:
 *     Trim low-quality tips of reads and patch low quality intervals within a sequence
 *     Trim adapter (and associated redundant prefix or suffix)
 *     Break chimers or all unscaffoldable no-coverage gaps of reads
 *
 ********************************************************************************************/

//  Analyze all the gaps between the good patches found in the first pass.
//  Consider a hole between two good intervals [lb,le] and [rb,re].  An overlap
  //    is anchored to the left of the whole if abpos <= le-COVER_LEN and aepos >= rb+COVER_LEN

static void PATCH_GAPS(int aread, Overlap *ovls, int novl)
{ static Patch dummy = { 0, 0, 0, 0, 0, 0, 0 };

#ifdef DEBUG_SUMMARY
  static char *status_string[4] = { "LOWQ", "SPAN", "SPLIT", "NOPAT" };
#endif

  Interval  lblock, rblock;
  Patch    *patch = NULL;
  int       status;
  int       tb, te;
  int       val;

#if defined(DEBUG_GAP_FILL) || defined(DEBUG_SUMMARY)
  printf("\n");
  printf("AREAD %d\n",aread);
#endif

  //  Determine patch for every LOWQ and SPAN gap and output dummy 0-patch
  //    for all SPLIT decisions

  tb = TRIM_IDX[aread];
  te = TRIM_IDX[aread+1];
  if (tb+2 < te)
    { lblock.beg = TRIM[tb];
      lblock.end = TRIM[tb+1];
      for (tb += 3; tb < te; tb += 3)
        { status     = TRIM[tb-1];
          rblock.beg = TRIM[tb];
          rblock.end = TRIM[tb+1];

          if (status == LOWQ)
            { patch = lowq_patch(ovls,novl,&lblock,&rblock);
              if (patch->bread < 0)
                status = SPAN;
            }
          if (status == SPAN)
            patch = span_patch(ovls,novl,&lblock,&rblock);

          if (status == SPLIT)
            { val = 0;
              patch = &dummy;
            }
          else
            { if (patch->bread < 0)
                { val = 0;
                  fpatch += 1;
#ifdef DEBUG_SUMMARY
                  TRIM[tb-1] = NOPAT;
#endif
                }
              else if (patch->comp)
                val = -(patch->bread+1);
              else
                val = patch->bread+1;
              npatch += 1;
            }
          fwrite(&val,sizeof(int),1,PR_DFILE);
          fwrite(&(patch->beg),sizeof(int),1,PR_DFILE);
          fwrite(&(patch->end),sizeof(int),1,PR_DFILE);
          PR_INDEX += 3*sizeof(int);

          lblock = rblock;
        }
    }
  fwrite(&PR_INDEX,sizeof(int64),1,PR_AFILE);

#ifdef DEBUG_SUMMARY
  tb = TRIM_IDX[aread];
  te = TRIM_IDX[aread+1];
#ifdef DEBUG_GAP_FILL
  if (tb+2 < te)
    printf("  FINAL:\n");
#endif
  if (tb < te)
    { printf("    [%d,%d]",TRIM[tb],TRIM[tb+1]);
      for (tb += 3; tb < te; tb += 3)
        printf(" %s [%d,%d]",status_string[TRIM[tb-1]],TRIM[tb],TRIM[tb+1]);
      printf("\n");
    }
#endif
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
  int   tbytes;

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
    tbytes = sizeof(uint8);
  else
    tbytes = sizeof(uint16);

  if (novl <= 0)
    return (0);

  Read_Overlap(input,ovls);
  if (trace)
    { if (ovls[0].path.tlen > pmax)
        { pmax  = 1.2*(ovls[0].path.tlen)+10000;
          paths = (uint16 *) Realloc(paths,sizeof(uint16)*pmax,"Expanding path buffer");
          if (paths == NULL) exit (1);
        }
      fread(paths,tbytes,ovls[0].path.tlen,input);
      if (tbytes == 1)
        { ovls[0].path.trace = paths;
          Decompress_TraceTo16(ovls);
        }
    }
  else
    fseek(input,tbytes*ovls[0].path.tlen,SEEK_CUR);

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
                  fread(paths+pcur,tbytes,ovls[n].path.tlen,input);
                  if (tbytes == 1)
                    { ovls[n].path.trace = paths+pcur;
                      Decompress_TraceTo16(ovls+n);
                    }
                }
              else
                fseek(input,tbytes*ovls[n].path.tlen,SEEK_CUR);
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
{ FILE       *input;
  char       *root, *dpwd;
  char       *las, *lpwd;
  int64       novl;
  HITS_TRACK *track;
  int         c;

  //  Process arguments

  { int  i, j, k;
    int  flags[128];

    ARG_INIT("DASpatch")

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("v")
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
  }

  //  Open trimmed DB and .qual and .trim tracks

  { int i, status;

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
        fread(&GOOD_QV,sizeof(int),1,afile);
        fread(&BAD_QV,sizeof(int),1,afile);
        fclose(afile);

        { int a, t, x;
          int tb, te;

          x = 0;
          for (a = 0; a < DB->nreads; a++)
            { tb = TRIM_IDX[a];
              te = TRIM_IDX[a+1];
              if (tb+2 < te)
                { if (TRIM[tb+2] == ADPRE)
                    tb += 3;
                  if (TRIM[te-3] == ADSUF)
                    te -= 3; 
                }
              TRIM_IDX[a] = x;
              for (t = tb; t < te; t++)
                TRIM[x++] = TRIM[t]; 
            }
          TRIM_IDX[DB->nreads] = x;
        }
      }
    else
      { fprintf(stderr,"%s: Must have a 'trim' track, run DAStrim\n",Prog_Name);
        exit (1);
      }
  }

  //  Initialize statistics gathering

  if (VERBOSE)
    { npatch = 0;
      fpatch = 0;

      printf("\nDASpatch -g%d -b%d %s",GOOD_QV,BAD_QV,argv[1]);
      for (c = 2; c < argc; c++)
        printf(" %s",argv[c]);
      printf("\n");
    }

  //  For each .las block/file

  dpwd = PathTo(argv[1]);
  root = Root(argv[1],".db");

  for (c = 2; c < argc; c++)
    { las  = Root(argv[c],".las");

      //  Determine if a .las block is being processed and if so get first and last read
      //    from .db file

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

      //  Set up QV trimming track

      { int len, size;
										      
        if (DB_PART > 0)
          { PR_AFILE = Fopen(Catenate(dpwd,PATHSEP,root,
                                      Numbered_Suffix(".",DB_PART,".patch.anno")),"w");
            PR_DFILE = Fopen(Catenate(dpwd,PATHSEP,root,
                                      Numbered_Suffix(".",DB_PART,".patch.data")),"w");
          }
        else
          { PR_AFILE = Fopen(Catenate(dpwd,PATHSEP,root,".patch.anno"),"w");
            PR_DFILE = Fopen(Catenate(dpwd,PATHSEP,root,".patch.data"),"w");
          }
        if (PR_AFILE == NULL || PR_DFILE == NULL)
          exit (1);
										      
        len  = DB_LAST - DB_FIRST;
        size = 8;
        fwrite(&len,sizeof(int),1,PR_AFILE);
        fwrite(&size,sizeof(int),1,PR_AFILE);
        PR_INDEX = 0;
	fwrite(&PR_INDEX,sizeof(int64),1,PR_AFILE);
      }

      //  Open overlap file

      lpwd = PathTo(argv[2]);
      if (DB_PART)
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
      ANCHOR_THRESH = ANCHOR_MATCH * TRACE_SPACING;

      make_a_pass(input,PATCH_GAPS,1);

      //  Clean up

      fclose(PR_AFILE);
      fclose(PR_DFILE);
    }

  //  If verbose output statistics summary to stdout

  if (VERBOSE)
    { if (fpatch == 0)
        printf("%s: All %d patches were successful\n",Prog_Name,npatch);
      else
        printf("%s: %d out of %d total patches failed\n",Prog_Name,fpatch,npatch);
    }

  free(dpwd);
  free(root);

  Close_DB(DB);
  free(Prog_Name);

  exit (0);
}
