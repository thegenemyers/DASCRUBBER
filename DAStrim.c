/*******************************************************************************************
 *
 *  Using overlap pile for each read and intrinisic quality values, determine the
 *    high quality segments with interspersed gaps.  Any unremoved
 *    adaptemer sequences are dectected and the shorter side trimmed.
 *    Every gap is analyzed and either patched or splits the read.
 *
 *  Author:  Gene Myers
 *  Date  :  June 2016
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

#undef    DEBUG_HQ_BLOCKS     //  Various DEBUG flags (normally all off)
#undef     SHOW_EVENTS
#undef    DEBUG_HOLE_FINDER
#undef    DEBUG_GAP_STATUS
#undef     SHOW_PAIRS
#undef    DEBUG_PATCHING
#undef    DEBUG_SUMMARY

#undef    DEBUG_CLASS

#define  ANNOTATE   //  Output annotation tracks for DaViewer


//  Command format and global parameter variables

static char *Usage = " [-v] [-g<int>] [-b<int>] <source:db> <overlaps:las> ...";

static int     COVERAGE;   //  estimated coverage
static int     BAD_QV;     //  qv >= and you are "bad"
static int     GOOD_QV;    //  qv <= and you are "good"
static int     HGAP_MIN;   //  less than this length do not process for HGAP
static int     VERBOSE;

//  Gap states

#define LOWQ  0   //  Gap is spanned by many LAs and patchable
#define SPAN  1   //  Gap has many paired LAs and patchable
#define SPLIT 2   //  Gap is a chimer or an unpatchable gap
#define ADAPT 3   //  Gap is due to adaptemer (internal only)

//  Good patch constants

#define MIN_BLOCK  500    //  Minimum length of a good patch

//  Gap constants

#define  MIN_COVER       3  //  A coverage gap occurs at or below this level
#define  COVER_LEN     400  //  An overlap covers a point if it extends COVER_LEN to either side.
#define  ANCHOR_MATCH  .25  //  Delta in trace interval at both ends of patch must be < this %.

#define  MIN_OVERLAP   900  //  Was COVER_LEN, too small?

//  Wall Constants

#define MIN_PNT 5        //  Minimum # of events in a wall
#define MAX_SEP 25       //  Maximum separation between two events in a wall
#define AVE_SEP 5.       //  Maximum average separation between two events in a wall


//  Global Variables (must exist across the processing of each pile)

  //  Input

static int     TRACE_SPACING;  //  Trace spacing (from .las file)

static DAZZ_DB _DB, *DB  = &_DB;   //  Data base
static int     DB_FIRST;           //  First read of DB to process
static int     DB_LAST;            //  Last read of DB to process (+1)
static int     DB_PART;            //  0 if all, otherwise block #

static int64  *QV_IDX;     //  qual track index
static uint8  *QV;         //  qual track values

  //  Output

static FILE   *TR_AFILE;   //  .trim.anno
static FILE   *TR_DFILE;   //  .trim.data
static int64   TR_INDEX;   //  Current index into .trim.data file as it is being written

#ifdef ANNOTATE

static FILE   *HQ_AFILE;   //  .hq.anno
static FILE   *HQ_DFILE;   //  .hq.data
static int64   HQ_INDEX;   //  Current index into .hq.data file as it is being written

static FILE   *HL_AFILE;   //  .hole.anno
static FILE   *HL_DFILE;   //  .hole.data
static int64   HL_INDEX;   //  Current index into .hole.data file as it is being written

static FILE   *SN_AFILE;   //  .span.anno
static FILE   *SN_DFILE;   //  .span.data
static int64   SN_INDEX;   //  Current index into .span.data file as it is being written

static FILE   *SP_AFILE;   //  .split.anno
static FILE   *SP_DFILE;   //  .split.data
static int64   SP_INDEX;   //  Current index into .split.data file as it is being written

static FILE   *AD_AFILE;   //  .adapt.anno
static FILE   *AD_DFILE;   //  .adapt.data
static int64   AD_INDEX;   //  Current index into .adapt.data file as it is being written

static FILE   *KP_AFILE;   //  .keep.anno
static FILE   *KP_DFILE;   //  .keep.data
static int64   KP_INDEX;   //  Current index into .keep.data file as it is being written

#endif

  //  Statistics

static int64 nreads, totlen;
static int64 nelim, nelimbp;
static int64 n5trm, n5trmbp;
static int64 n3trm, n3trmbp;
static int64 natrm, natrmbp;
static int64 ngaps, ngapsbp;
  static int64 nlowq, nlowqbp;
  static int64 nspan, nspanbp;
  static int64 nchim, nchimbp;


//  Data Structures

typedef struct   //  General read interval [beg..end]
  { int beg;
    int end;
  } Interval;

  //  Coverage events, type (one of 7 below) and position

#define ADD  0  //  leftmost A-position of LA
#define LFT  1  //  ADD position + COVER_LEN of LA (>= 2*COVER_LEN long)
#define LGP  2  //  left end of an HQ-block
#define CTR  3  //  A-center of LA < 2*COVER_LEN long
#define RGP  4  //  right end of an HQ-block
#define RGT  5  //  DEL position - COVER_LEN of LA
#define DEL  6  //  rightmost A-position of LA

#ifdef SHOW_EVENTS

static char   Symbol[7] = { 'A', 'L', '[', 'C', ']', 'R', 'D' };

#endif

typedef struct
  { int  type;   
    int  pos;
  } Event;

  //  Wall: there are cnt LFT/RGT events ending in the interval [beg,end] going
  //        from coverage depth cov up to cov+cnt

typedef struct
  { int beg;
    int end;
    int cnt;
    int cov;
  } Wall;

/*******************************************************************************************
 *
 *  FIND ALL HIGH_QV BLOCKS OF EACH READ
 *
 ********************************************************************************************/

  //  Find "good" blocks of trace point intervals:
  //    0. A good block must begin and end with an interval <= GOOD_QV
  //    1. Any stretch all < BAD_QV at least MIN_BLOCK long
  //    2. Any stretch all <= GOOD_QV at least MIN_BLOCK-TRACE_SPACING long
  //    3. Any stretch all <= GOOD_QV only 1 interval away from another good patch

  //  Global Inputs:  QV, QV_IDX, GOOD_QV, BAD_QV 
  //  HQ_BLOCKS[0..*nblk-1] contain the good patches in increase sequencing across aread.
  //  Parameter aread is input-only, and p_nblk is output-only.

static Interval *HQ_BLOCKS(int aread, int *p_nblk)
{ int              nblk;
  static int      *alive = NULL;
  static Interval *block = NULL;

  int    alen, atick;
  uint8 *qvec;

  alen  = DB->reads[aread].rlen;
  atick = (alen + (TRACE_SPACING-1))/TRACE_SPACING;

  if (alive == NULL)
    { int max = DB->maxlen/TRACE_SPACING+2;
      alive = (int *) Malloc(max*sizeof(int),"Allocating alive vector");
      block = (Interval *) Malloc(max*sizeof(Interval),"Allocating block vector");
      if (alive == NULL || block == NULL)
        exit (1);
    }

  qvec = QV + QV_IDX[aread];
  nblk = 0;

  //  Find all blocks < BAD_QV with either len >= MIN_BLOCK or all <= GOOD_QV in block[0..nblk)
  //    Mark those satisfying 1. or 2. as "alive" (.alv)

  { int lmost = 0, rmost = 0, thr;
    int i, in;

    thr = (MIN_BLOCK-1)/TRACE_SPACING;
    in  = 0;
    for (i = 0; i <= atick; i++)
      { int q, alv;

        if (i < atick)
          q = qvec[i];
        else
          q = BAD_QV;
        if (in)
          { if (q >= BAD_QV)
              { alv = (lmost-rmost >= thr);
                if (alv)
                  { block[nblk].beg = rmost; 
                    block[nblk].end = lmost + 1;
                    alive[nblk] = alv;
                    nblk += 1;
                  }
                else
                  { int j, k;

                    for (j = rmost; j <= lmost; j = k)
                      { for (k = j+1; k <= lmost; k++)
                          if (qvec[k] > GOOD_QV)
                            break;
                        block[nblk].beg = j; 
                        block[nblk].end = k;
                        alive[nblk] = (k-j >= thr);
                        nblk += 1;
                        for ( ; k <= lmost; k++)
                          if (qvec[k] <= GOOD_QV)
                            break;
                      }
                  }
                in = 0;
              }
            else if (q <= GOOD_QV)
              lmost = i;
          }
        else
          { if (q <= GOOD_QV)
              { rmost = lmost = i;
                in = 1;
              }
          }
      }
  }

  //  Mark as alive all short, all-good blocks that satisfy 3.

  { int i, j;

    for (i = 0; i < nblk; i++)
      if (alive[i])
        { for (j = i-1; j >= 0 && ! alive[j]; j--)
            if (block[j+1].beg - block[j].end == 1)
              alive[j] = 1;  
            else
              break;
          for (j = i+1; j < nblk && ! alive[j]; j++)
            if (block[j].beg - block[j-1].end == 1)
              alive[j] = 1;  
            else
              break;
        }
  }

  //  Remove all blocks that are not alive

  { int i, j;

    j = 0;
    for (i = 0; i < nblk; i++)
      if (alive[i])
        { block[j].beg = block[i].beg * TRACE_SPACING;
          block[j].end = block[i].end * TRACE_SPACING;
          j += 1;
        }
    nblk = j;
    if (nblk > 0 && block[nblk-1].end > alen)
      block[nblk-1].end = alen;
  }

#ifdef DEBUG_HQ_BLOCKS
  { int i;

    printf("         %3d:",nblk);
    for (i = 0; i < nblk; i++)
      printf(" [%5d,%5d]",block[i].beg,block[i].end);
    printf("\n");
  }
#endif

  *p_nblk = nblk;
  return (block);
}


/*******************************************************************************************
 *
 *  WALL ANALYZER TO HELP AVOID REPEAT BOUNDARIES
 *
 ********************************************************************************************/

   //  Find intervals of LFT/RGT events where no two events are separated by more than
   //    MAX_SEP, the average arrival rate is AVE_SEP, and there are at least MIN_PNT
   //    events in the interval.

static Wall *wall_detector(int *ev, int b, int e, Wall *next)
{ int idx;

  { int    i, n, max;
    double ave;

    n = e-b;

    if (n < MIN_PNT) return (next);  //  Too small: done

    idx = b;
    max = -1;                        //  Find the position of the largest separation between
    for (i = b+1; i < e; i++)        //     two tips in ev[b..e)
      if (ev[i] - ev[i-1] > max)
        { max = ev[i] - ev[i-1];
          idx = i;
        }

    ave = (ev[e-1] - ev[b]) / (n-1.);        //  Check if the current interval is a wall
    if (ave <= AVE_SEP && max <= MAX_SEP)
      { if (max <= 4.*(ave+1.))              //  Max separation < 4*average separation ?
          { next->beg = b;
            next->end = e;
            next->cnt = n;
            return (next+1);
          }
      }
  }

  next = wall_detector(ev,b,idx,next);    //  If not then split on the largest separation
  next = wall_detector(ev,idx,e,next);    //    and recurse on the two parts
  return (next);
}

   //  Find LFT/RGT event walls

static Wall *find_walls(int novl, Event *queue, int *anum, int *dnum)
{ static int   nmax = 0;

  Wall        *aptr, *dptr;
  static Wall *wall = NULL;

  int          ntip;
  static int  *adds = NULL;
  static int  *dels;

  if (novl == 0)
    return (NULL);

  if (novl > nmax)
    { nmax = novl*1.2 + 1000; 
      wall = (Wall *) Realloc(wall,sizeof(Wall)*2*(nmax/MIN_PNT),"Reallocating wall vector");
      adds = (int *)  Realloc(adds,sizeof(int)*2*nmax,"Reallocating add+del vectors");
      if (wall == NULL || adds == NULL)
        exit (1);
      dels = adds + nmax;
    }

  //  Make separate arrays of add and del tips (LFT and RGT events) in sorted order in
  //    which to seek "walls".

  { int i, j, x;

    i = x = 0;                      //  A bit tricky: less than novl tips due to CTR events
    for (j = 0; x < novl; j++)      //   that don't generate tips, so analyze events until
      if (queue[j].type == CTR)     //   have counted all LA's.  Furthermore adds and dels
        x += 1;                     //   are sorted because queue is sorted.
      else if (queue[j].type == LFT)
        { x += 1;
          adds[i++] = queue[j].pos;
        }
    ntip = i;

    i = 0;
    for (j = 0; i < ntip; j++)
      if (queue[j].type == RGT)
        dels[i++] = queue[j].pos;
  }

  //  Find LFT walls and RGT walls in [walls,aptr) and [aptr,dptr)

  aptr = wall_detector(adds,0,ntip,wall);
  dptr = wall_detector(dels,0,ntip,aptr);

  //  For each wall, determine the coverage of its base with a merged traversal
  //    of the adds and dels arrays

  { Wall *a, *d;
    int   i, j, x;

    x = 0;
    a = wall;
    d = aptr;;
    i = j = 0;
    while (j < ntip)
      if (i < ntip && adds[i] < dels[j])
        { if (a->beg == i)
            a->cov = x;
          else if (a->end == i+1)
            { a += 1;
              if (a >= aptr)
                a -= 1;
            }
          x += 1;
          i += 1;
        }
      else
        { if (d->beg == j)
            d->cov = x - d->cnt;
          else if (d->end == j+1)
            { d += 1;
              if (d >= dptr)
                d -= 1;
            }
          x -= 1;
          j += 1;
        }
  }
  
  //  Sneaky, switch beg/end from an index into the adds or dels array, to the actually
  //    coordinate of the event.

  { Wall *a;

    for (a = wall; a < aptr; a++)
      { a->beg = adds[a->beg];
        a->end = adds[a->end-1];
      }
    for (a = aptr; a < dptr; a++)
      { a->beg = dels[a->beg];
        a->end = dels[a->end-1];
      }
  }

  *anum = aptr-wall;
  *dnum = dptr-aptr;
  return (wall);
}


/*******************************************************************************************
 *
 *  COVERAGE ANALYSIS TO FIND ALL HOLES (regions of very low coverage/support)
 *
 ********************************************************************************************/

  //  Find intervals for which there are MIN_COVER or fewer LAs that project at least COVER_LEN
  //    bases to the left and right of the interval.  These are called holes.
  //  Holes are usually found between HQ-blocks.  However occasionally they intersect one or
  //    more blocks and this requires the HQ-blocks be refined as follows:
  //      a. Hole spans an HQ-block:
  //           The block needs to be removed as HQ *if* it is not based on 5 or more LA's
  //           (this usually never happens, 10^-5 or less)
  //      b. Hole is contained in an HQ-block:
  //           The block needs to be split around the hole because one needs to verify that
  //           the left and right regions on each side of a hole actually belong together
  //           (this happens occasionaly, ~ 10^-3)
  //      c. Hole overlaps an HQ-block:
  //           If this happens, then the overlap is very small and the block is left unperturbed.
  //           (this worries me a bit, but in all testing it (very small overlap) remains so)
  //  Given the above possibilities, the list of HQ-blocks can be modified by FIND_HOLES.

static int ESORT(const void *l, const void *r)
{ Event *x = (Event *) l;
  Event *y = (Event *) r;

  if (x->pos == y->pos)
    return (x->type - y->type);
  return (x->pos - y->pos);
}

static int FIND_HOLES(int aread, Overlap *ovls, int novl, Interval *block, int nblk)
{ static int       nmax = 0;

  int              nev;
  static Event    *queue = NULL;  //  Event queue[0..nev)

  int              nhole;
  static Interval *holes = NULL;  //  Detected holes[0..nhole)

  static int       pmax;
  static Interval *cover = NULL;  //  Coverage at block ends [0..nblk)
  static Interval *nwblk;         //  Modified block list [0..nblk')

  int   anum = 0, dnum = 0;       //  LFT and RGT walls, awall[0..anum) & dwall[0..dnum)
  Wall *awall,   *dwall;

  if (cover == NULL)
    { pmax = DB->maxlen/TRACE_SPACING + 2;
      cover = (Interval *) Malloc(2*pmax*sizeof(Interval),"Allocating patch vector");
      nwblk = cover + pmax;
    }

  if (4*novl + pmax > nmax)
    { nmax = 4.8*novl + pmax + 100;
      queue = (Event *) Realloc(queue,(nmax+1)*sizeof(Event),"Allocating event queue");
      holes = (Interval *) Realloc(holes,(nmax/4)*sizeof(Interval),"Allocating hole vector");
      if (queue == NULL || holes == NULL)
        exit (1);
    }

  { int i;

    //  For each trimmed overlap: add its events to the queue

    nev = 0;
    for (i = 0; i < novl; i++)
      { queue[nev].type = ADD;
        queue[nev].pos  = ovls[i].path.abpos;
        nev += 1;

        queue[nev].type = DEL;
        queue[nev].pos  = ovls[i].path.aepos;
        nev += 1;

        if (ovls[i].path.abpos + 2*COVER_LEN + 10 > ovls[i].path.aepos)
          { queue[nev].type = CTR;
            queue[nev].pos  = (ovls[i].path.abpos + ovls[i].path.aepos) / 2;
            nev += 1;
          }
        else
          { queue[nev].type = LFT;
            queue[nev].pos  = ovls[i].path.abpos + COVER_LEN;
            nev += 1;
            queue[nev].type = RGT;
            queue[nev].pos  = ovls[i].path.aepos - COVER_LEN;
            nev += 1;
          }
      }

    //  For each HQ-block: add its events to the queue

    for (i = 0; i < nblk; i++)
      { queue[nev].type = LGP;
        queue[nev].pos  = block[i].beg;
        nev += 1;
        queue[nev].type = RGP;
        queue[nev].pos  = block[i].end;
        nev += 1;
      }

    queue[nev].pos = DB->reads[aread].rlen;
  }

  //  Sort the events

  qsort(queue,nev,sizeof(Event),ESORT);

  //  Find all LFT and RGT walls

  awall = find_walls(novl,queue,&anum,&dnum);
  dwall = awall + anum;

#ifdef DEBUG_HOLE_FINDER
  { int i;

    printf("\n");
    for (i = 0; i < anum; i++)
      printf("  Add [%5d,%5d] %d %d\n",awall[i].beg,awall[i].end,awall[i].cnt,awall[i].cov);
    for (i = 0; i < dnum; i++)
      printf("  Del [%5d,%5d] %d %d\n",dwall[i].beg,dwall[i].end,dwall[i].cnt,dwall[i].cov);
    printf("\n");
  }
#endif

  //  Move through events in order keeping track of inc, dec, & cnf so that the
  //    invariant stated below holds

  { int  cnf, inc, dec;
    int  cblk;

    int   in;
    int   nbeg, nend = 0;
    int   first, last;

    int  i;

    in    = 1;
    first = -1;
    cblk  = 0;
    nhole = 0;
    inc = dec = cnf = 0;
    for (i = 0; i < nev; i++)
      { switch (queue[i].type)
        { case ADD:
            inc += 1;
            break;
          case LFT:
            inc  -= 1;
            cnf  += 1;
            break;
          case LGP:
            cover[cblk].beg = cnf + inc + dec;  // = coverage depth at block[cblk].beg
            continue;
          case CTR:
            inc -= 1;
            dec += 1;
            continue;
          case RGP:
            cover[cblk].end = cnf + inc + dec;  // = coverage depth at block[cblk].end
            cblk += 1;
            continue;
          case RGT:
            cnf  -= 1;
            dec  += 1;
            break;
          case DEL:
            dec -= 1;
            break;
        }

        //  For position x = queue[i].pos:
        //    inc = # of LA's between (ADD,LFT] positions
        //    dec = # of LA's between (RGT,DEL] positions
        //    cnf = # of LA's between (LFT,RGT] positions (= # of LAs tat project at least
        //             COVER_LEN bases to the right and left of x!

#ifdef SHOW_EVENTS
        printf(" %5d %c: %3d<  %3d  >%3d  %3d\n",
               queue[i].pos,Symbol[queue[i].type],inc,cnf,dec,dec-inc);
#endif

        //  When truncated coverage, cnf, transitions below MIN_COVER(3), note the fact (in = 1)
        //    and record the index first of the event (must be a RGT) and the number of LA's
        //    currently in their (RGT,DEL] interval

        if (cnf <= MIN_COVER)
          { if ( ! in)
              { in = 1;
                nend  = dec;
                first = i;
              }
          }

        //  When truncated coverage transitions above MIN_COVER, we declare it a hole
        //    if interval below MIN_COVER is at least COVER_LEN long, there are at least
        //    4 LA's that are "ending" at the left (i.e. in (RGT,DEL] interval, and
        //    at least 4 LA's ending at the right.

        else
          { if (in && first >= 0 && queue[i].pos - queue[first].pos >= COVER_LEN &&
                  nend >= 4 && inc >= 4)
              { int lflank, rflank;
                int dpos, apos;

                nbeg = inc;
                last = i;

                //  Need to find the boundaries of the hole.  In principle, this is
                //    [dpos + COVER_LEN, apos - COVER_LEN] where apos = queue[first].pos
                //    and dpos = queue[last].pos, i.e. the entry and exit into the low
                //    truncated cover interval.  However, walls induced by repeat boundaries
                //    and/or uneveness in the end-points of LA's can cause the above to be
                //    quite far off.  So ...


                //  First try the average of the 2nd and 3rd quartile of the nend RGT events
                //    before dpos.   The requisite number of events must exist by the definition
                //    of nend.  While one is at it determine the index of the first of the
                //    nend RGT events in lflank.

                { int64 sum;
                  int   q1, q3, n;
                  int   a, d, k;
                  int   acov, dcov;

                  q1  = nend/4;
                  q3  = (3*nend)/4;
                  sum = 0;
                  n   = 0;
                  for (lflank = first; n < nend; lflank--)
                    if (queue[lflank].type == RGT || queue[lflank].type == CTR)
                      { if (n >= q1 && n < q3)
                          sum += queue[lflank].pos;
                        n += 1;
                      }
                  dpos = sum/(q3-q1);
                  lflank += 1;

#ifdef DEBUG_HOLE_FINDER
                  printf("  Dev %5d-%3d-%5d -> %5d",queue[lflank].pos,nend,queue[first].pos,dpos);
#endif

                //  Second, look for the leftmost RGT-(LFT-)wall that overlaps the left (right)
                //    flank, i.e. queue[lflank,first].pos (queue[last,rflank].pos), and if found
                //    take the average position of the flank position within the wall.

                  for (d = dnum-1; d >= 0; d--)
                    if (dwall[d].beg <= queue[first].pos)
                      break;
                  if (d >= 0 && dwall[d].end >= queue[lflank].pos)
                    { sum = 0;
                      n   = 0;
                      for (k = first; k >= lflank; k--)
                        if (queue[k].type == RGT || queue[k].type == CTR)
                          { if (queue[k].pos < dwall[d].beg)
                              break;
                            if (queue[k].pos <= dwall[d].end)
                              { sum += queue[k].pos;
                                n += 1;
                              }
                          }
                      dpos = sum/n;
#ifdef DEBUG_HOLE_FINDER
                      printf(" Map [%5d,%5d] -> %4d\n",dwall[d].beg,dwall[d].end,dpos);
#endif
                      dcov = dwall[d].cov + dwall[d].cnt;
		      d -= 1;
		    }
                  else
                    { dcov = nend + MIN_COVER;
#ifdef DEBUG_HOLE_FINDER
                      printf(" No wall mapping\n");
#endif
                    }

                // First try on LFT events (replace nend with nbeg, RGT with LFT, before
                //   with after, and dpos with apos, first with last, and lflank with rflank.

                  q1  = nbeg/4;
                  q3  = (3*nbeg)/4;
                  sum = 0;
                  n   = 0;
                  for (rflank = last; n < nbeg; rflank++)
                    if (queue[rflank].type == LFT || queue[rflank].type == CTR)
                      { if (n >= q1 && n < q3)
                          sum += queue[rflank].pos;
                        n += 1;
                      }
                  apos = sum/(q3-q1);
                  rflank -= 1;

#ifdef DEBUG_HOLE_FINDER
                  printf("  Aev %5d-%3d-%5d -> %5d",queue[i].pos,nbeg,queue[rflank].pos,apos);
#endif

                // Second look at LFT events.

                  for (a = 0; a < anum; a++)
                    if (awall[a].end >= queue[i].pos)
                      break;
                  if (a < anum && awall[a].beg <= queue[rflank].pos) 
                    { sum = 0;
                      n   = 0;
                      for (k = i; k <= rflank; k++)
                        if (queue[k].type == LFT || queue[k].type == CTR)
                          { if (queue[k].pos > awall[a].end)
                              break;
                            if (queue[k].pos >= awall[a].beg)
                              { sum += queue[k].pos;
                                n += 1;
                              }
                          }
                      apos = sum/n;
#ifdef DEBUG_HOLE_FINDER
                      printf(" Map [%5d,%5d] -> %4d\n",awall[a].beg,awall[a].end,apos);
#endif
                      acov = awall[a].cov + awall[a].cnt;
                      a += 1;
                    }
                  else
                    { acov = nbeg + MIN_COVER;
#ifdef DEBUG_HOLE_FINDER
                      printf(" No wall mapping\n");
#endif
                    }

                // If apos and dpos are still so close that the implied hole boundaries
                //   are out of order by 50 or more bases, then walk back through ascending
                //   walls (if present) until this is no longer true or there are no more
                //   more walls left.  If both left and right options exist, always take
                //   the wall starting at the lower current height.

                  while (apos - dpos < 2*COVER_LEN - 50)
                    { if (d >= 0 && dwall[d].cov >= dcov)
                        if (a < anum && awall[a].cov >= acov)
                          { if (dcov < acov)
                              { dcov = dwall[d].cov + dwall[d].cnt;
                                dpos = dwall[d--].beg;
#ifdef DEBUG_HOLE_FINDER
                                printf("  Push <- %d\n",dpos);
#endif
                              }
                            else
                              { acov = awall[a].cov + awall[a].cnt;
                                apos = awall[a++].end;
#ifdef DEBUG_HOLE_FINDER
                                printf("  Push -> %d\n",apos);
#endif
                              }
                          }
                        else
                          { dcov = dwall[d].cov + dwall[d].cnt;
                            dpos = dwall[d--].beg;
#ifdef DEBUG_HOLE_FINDER
                            printf("  Push <- %d\n",dpos);
#endif
                          }
                      else
                        if (a < anum && awall[a].cov >= acov)
                          { acov = awall[a].cov + awall[a].cnt;
                            apos = awall[a++].end;
#ifdef DEBUG_HOLE_FINDER
                            printf("  Push -> %d\n",apos);
#endif
                          }
                        else
                          {
#ifdef DEBUG_HOLE_FINDER
                            printf("  FAULT\n");
#endif
                            break;
                          }
                    }
                }

                // Finalize and record the hole boundaries.

                holes[nhole].beg = dpos + COVER_LEN;
                holes[nhole].end = apos - COVER_LEN;
                nhole += 1;
              }
            in   = 0;
          }
      }
  }

  //  See if the holes remove or split any HQ-blocks and build the revised list
  //    in newblk[0..q).

  { int i, p, q, x;
    int lhang, rhang;
#ifdef DEBUG_HOLE_FINDER
    int reverse;
#endif

    //  For each hole in left-to-right order

    p = q = 0;
    for (i = 0; i < nhole; i++)
      { if (holes[i].beg > holes[i].end)
          { x = holes[i].beg;
            holes[i].beg = holes[i].end;
            holes[i].end = x;
#ifdef DEBUG_HOLE_FINDER
            reverse = 1;
#endif
          }
#ifdef DEBUG_HOLE_FINDER
        else
          reverse = 0;
#endif

        //  Advance to the next block p that intersects with or is to the right of hole
        //    moving blocks being skipped over to the new block list

        while (p < nblk && block[p].end <= holes[i].beg)
          nwblk[q++] = block[p++];

#ifdef DEBUG_HOLE_FINDER
        printf("  HOLE: %5d [%5d,%5d]\n",
               aread+1,holes[i].beg,holes[i].end);
#endif
      
        //  While the current block intersects the current hole

        while (p < nblk && block[p].beg < holes[i].end)
          { lhang = (holes[i].beg < block[p].beg);
            rhang = (holes[i].end > block[p].end);
            if (lhang)
              { if (rhang)

                  //  Hole i contains block p: remove it if coverage <= 4 at both ends

                  { if (block[p].end - block[p].beg >= MIN_BLOCK &&
                          (cover[p].beg > 4 || cover[p].end > 4))
                      nwblk[q++] = block[p];
                    p += 1;
#ifdef DEBUG_HOLE_FINDER
                    printf("      INTERSECT %5d S [%5d,%5d] %3d %3d",
                           aread+1,block[p-1].beg,block[p-1].end,cover[p-1].beg,cover[p-1].end);
                    if (reverse)
                      printf(" REV");
                    printf("\n");
#endif
                  }

                 //  Hole i intersect the left tip of block p: nothing to do

                else
                  {
#ifdef DEBUG_HOLE_FINDER
                    printf("      INTERSECT %5d Z %5d [..,%5d] %3d",
                           aread+1,holes[i].end-block[p].beg,holes[i].end,cover[p].beg);
                    if (reverse)
                      printf(" REV");
                    printf("\n");
#endif
                    break;
                  }
              }
            else
              if (rhang)

                //  Hole i intersect the right tip of block p: move p to new block list

                { nwblk[q++] = block[p++];
#ifdef DEBUG_HOLE_FINDER
                  printf("      INTERSECT %5d Z %5d [%5d,..] %3d",
                         aread+1,block[p-1].end-holes[i].beg,holes[i].beg,cover[p-1].end);
                  if (reverse)
                    printf(" REV");
                  printf("\n");
#endif
                }
              else

                //  Hole i is contained within block p:  Break block into two parts at
                //    TRACE_SPACING ticks left and right of hole, and keep each piece
                //    if they are greater than MIN_BLOCK long.

                { int beg, end;

#ifdef DEBUG_HOLE_FINDER
                  printf("      INTERSECT %5d C %5d [%5d,%5d]",
                         aread+1,holes[i].end-holes[i].beg,block[p].beg,block[p].end);
                  if (reverse)
                    printf(" REV");
		  printf("\n");
#endif
                  beg = (holes[i].beg/TRACE_SPACING);
                  end = (holes[i].end-1)/TRACE_SPACING+1;
                  if (beg == end)
                    { beg -= 1; end += 1; }
                  beg *= TRACE_SPACING;
                  end *= TRACE_SPACING;
                  if (beg - block[p].beg >= MIN_BLOCK)
                    { nwblk[q].beg = block[p].beg;
                      nwblk[q++].end = beg;
                    }
                  if (block[p].end - end >= MIN_BLOCK)
                    block[p].beg = end;
                  else
                    p += 1;
                  break;
                }
          }
      }

    //  Remove any remaining blocks to the new list

    while (p < nblk)
      nwblk[q++] = block[p++];
    nblk = q;

    //  Transfer new blocks to original block vector

    for (i = 0; i < nblk; i++)
      block[i] = nwblk[i];
  }

#ifdef ANNOTATE
  { int i;

    for (i = 0; i < nhole; i++)
      if (holes[i].end - holes[i].beg < 75)
        { holes[i].end += 50;
          holes[i].beg -= 50;
          fwrite(&(holes[i].beg),sizeof(int),1,HL_DFILE);
          fwrite(&(holes[i].end),sizeof(int),1,HL_DFILE);
          holes[i].end -= 50;
          holes[i].beg += 50;
        }
      else
        { fwrite(&(holes[i].beg),sizeof(int),1,HL_DFILE);
          fwrite(&(holes[i].end),sizeof(int),1,HL_DFILE);
        }
    HL_INDEX += 2*nhole*sizeof(int);
    fwrite(&HL_INDEX,sizeof(int64),1,HL_AFILE);
  }
#endif

  // Return the list of holes holes[0..nhole) and the new list of blocks, nwblk[0..nblk)

  return (nblk);
}


/*******************************************************************************************
 *
 *  FIND ANY UNREMOVED ADAPTER (OR POLYMERASE SWITCHES) AND TRIM SMALLER PARTS
 *
 ********************************************************************************************/

typedef struct
  { int lidx;    //  left LA index
    int ridx;    //  right LA index
    int delta;   //  Difference between A-gap and B-gap
    int soft;    //  0 (soft) = pair is not close to gap border or adjacent gap border on both sides
                 //  1 (anciliary) = pair is not close to gap border but is to adjacent gap border
                 //  2 (hard) = pair is close to gap border on both sides
  } Spanner;

typedef struct
  { int bread;   //  bread^comp[beg..end] is the patch sequence
    int comp;
    int beg;
    int end;
    int anc;     //  maximum anchor interval match
    int bad;     //  number of segments that are bad
    int avg;     //  average QV of the patch
  } Patch;

static int GSORT(const void *l, const void *r)
{ Spanner *x = (Spanner *) l;
  Spanner *y = (Spanner *) r;

  return (x->delta - y->delta);
}

//  Return match score of lov->bread with "anchor" lov->aread[lft-TRACE_SPACING,lft]

static int eval_lft_anchor(int lft, Overlap *lov)
{ uint16  *tr;
  int      te;

  if (lft > lov->path.aepos)
    return (50);
  tr = (uint16 *) lov->path.trace;
  te = 2 * (((lft + (TRACE_SPACING-1)) - lov->path.abpos)/TRACE_SPACING);
  if (te <= 0)
    return (50);
  return (tr[te-2]);
}

//  Return match score of lov->bread with "anchor" lov->aread[rgt,rgt+TRACE_SPACING]

static int eval_rgt_anchor(int rgt, Overlap *rov)
{ uint16  *tr;
  int      te;

  if (rgt < rov->path.abpos)
    return (50);
  tr = (uint16 *) rov->path.trace;
  te = 2 * (((rgt + (TRACE_SPACING-1)) - rov->path.abpos)/TRACE_SPACING);
  if (te >= rov->path.tlen)
    return (50);
  return (tr[te]);
}

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

  bread = lov->bread;
  bcomp = COMP(lov->flags);
  blen  = DB->reads[bread].rlen;

  if (blen < HGAP_MIN)
    return (NULL);

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

  ans.bread = bread;
  ans.comp  = bcomp;
  ans.beg   = bb;
  ans.end   = be;
  ans.anc   = anc;

  //  Compute metrics for b-read patch

  if (bcomp)
    { t  = blen - be;
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

  ans.bad   = bad;
  ans.avg   = avg;

  return (&ans);
}

//  Examine the spanning pairs for a gap.  Group those with sufficient density
//    i.e. with 20 + 10% of the last one.  If test == 1, keep groups that have at
//    least 4 members, and are either 60% hard or at least 8 hard pairs, but trim
//    away any extremal non-hard pairs.  If test == 0, keep the largest group that has
//    at least 4 members, and is either 60% not soft or at least 8 non-soft pairs,
//    but trim awya ny extermal soft pairs.  If "move" is non-zero then compress gsort
//    accordingly, return the size of gsort after trimming in all instances.

static int analyze_gap_pairs(int gsize, Spanner *gsort, Overlap *ovls,
                             int gcnt, int scnt, int test, int move)
{ int j, l, c, x, w;
  int bord, soft, keeper;
  int ncnt, biggest;

  (void) ovls;

  biggest = 0;
  ncnt = 0;
#ifdef SHOW_PAIRS
  if (move)
    printf("  Gsort: %d\n",gsize);
#endif
  if (gcnt > 0)
    c = gsize - gsort[0].delta;
  else
    c = gsize;
  w = 0;
  for (j = 0; j <= gcnt; j++)
    { l = c;
      if (j >= gcnt)
        bord = 1;
      else
        { c = gsize - gsort[j].delta;
          if (l < 0)
            bord = (l-c >= 20-.1*c);
          else
            bord = (l-c >= 20+.1*l);
        }
      if (bord)
        { soft = 0;
          for (x = w; x < j; x++)
            soft += (gsort[x].soft <= test);
          keeper = (j-w >= 4 && (soft < .4*(j-w) || (j-w)-soft >= 8));
          if (test == 0)
            { if (keeper && j-w > biggest)
                { biggest = j-w;
                  ncnt = 0;
                }
              else
                keeper = 0;
            }
#ifdef SHOW_PAIRS
          if (move)
            { printf("----\n");
              for (x = w; x < j; x++)
                { printf("  %3d: %5d %5d",x,gsort[x].delta,gsize-gsort[x].delta);
                  printf(" %5d",ovls[gsort[x].lidx].bread);
                  if (gsort[x].soft == 0)
                    printf(" @");
                  else if (gsort[x].soft == 1)
                    printf(" #");
                  else
                    printf("  ");
                  if (!keeper)
                    printf(" X");
                  printf("\n");
                }
            }
#endif
          if (keeper)
            { for (x = w; x < j; x++)
                if (gsort[x].soft > test)
                  break;
              for (w = j; gsort[w-1].soft <= test; w--)
                ;
              if (move)
                while (x < w)
                  gsort[ncnt++] = gsort[x++];
              else
                ncnt += w-x;
            }
          w = j;
        }
    }
  if (move)
    for (x = gcnt; x < gcnt+scnt; x++)
      gsort[ncnt++] = gsort[x];
  else
    ncnt += scnt;

  return (ncnt);
}

//  Categorize each gap and if appropriate return the best patch for each

static int gap_status(Overlap *ovls, int novl, Interval *lblock, Interval *rblock,
                      int *p_lft, int *p_rgt)
{ static int  nmax = 0;

  static Spanner *gsort = NULL;   //  A-B delta and idx-pair for all B-reads spanning a gap

  static int      ANCHOR_THRESH;

  static Interval *FirstB;
  static Interval *LastB;

  int j;
  int lft, rgt;
  int lcv, rcv; 
  int cnt;

  if (p_lft == NULL)
    { if (novl > nmax)
        { nmax = 1.2*novl + 500;
          gsort = (Spanner *) Realloc(gsort,nmax*sizeof(Spanner),"Allocating gap vector");
          if (gsort == NULL)
            exit (1);
          ANCHOR_THRESH = ANCHOR_MATCH * TRACE_SPACING;
        }

      FirstB = lblock;
      LastB  = rblock-1;
      return (0);
    }

  lft = lblock->end;
  rgt = rblock->beg;
  lcv = lft - COVER_LEN;
  rcv = rgt + COVER_LEN;
  if (lcv < lblock->beg)
    lcv = lblock->beg;
  if (rcv > rblock->end)
    rcv = rblock->end;

#ifdef DEBUG_GAP_STATUS
  printf("  GAP [%5d,%5d]  <%5d,%5d>\n",lft,rgt,lcv,rcv);
#endif

  //  If the gap flank [lcv,rcv] is covered by 10 or more LAs, then a LOWQ gap

  cnt = 0;
  for (j = 0; j < novl; j++)
    if (ovls[j].path.abpos <= lcv && ovls[j].path.aepos >= rcv)
      { cnt += 1;
        if (cnt >= 10)
          break;
      }

  //  If so and it is patchable then report LOWQ

  if (cnt >= 10)
    { for (j = 0; j < novl; j++)
        if (ovls[j].path.abpos <= lcv && ovls[j].path.aepos >= rcv)
          { Patch *can;

            can = compute_patch(lft,rgt,ovls+j,ovls+j);
            
            if (can == NULL) continue;

            if (can->anc <= ANCHOR_THRESH && can->avg <= GOOD_QV && can->bad == 0)
              {
#ifdef DEBUG_GAP_STATUS
                printf("    LOWQ  PATCHABLE = %d%c[%d..%d]  %d (%d)\n",
                       can->bread,can->comp?'c':'n',can->beg,
                       can->end,can->anc,can->avg);
#endif
                return (LOWQ);
              }
          }
#ifdef DEBUG_GAP_STATUS
      printf("    FAILING TO PATCH_LOWQ\n");
#endif
    }

  { int bread, bcomp, blen, blast;
    int ab, ae;
    int lstack[10], ltop;
    int rstack[10], rtop;
    int lcnt, rcnt, scnt, gcnt, acnt;
    int lidx, ridx, sidx, Lidx, Ridx;
    int k;

    //  Find LA pairs or LAs spanning the gap flank [lcv,rcv]

    bread = -1;
    lidx = ridx = -1;
    lcnt = rcnt = scnt = gcnt = acnt = 0;
    for (j = 0; j < novl; j = k)
      { blast = bread;
        bread = ovls[j].bread;
        blen  = DB->reads[bread].rlen;
        bcomp = COMP(ovls[j].flags);
        Lidx  = lidx;
        Ridx  = ridx;

        ltop = rtop = 0;
        lidx = ridx = sidx = -1;    //  For all LA's with same b-read
        for (k = j; k < novl; k++)
          { if (ovls[k].bread != bread)
              break;
            if (COMP(ovls[k].flags) != (uint32) bcomp)   //  Note when b switches orientation
              break;
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
            //    A partner is hard if end=point is within COVER_LEN of the gap boundary
            //    Record rigthmost/leftmost left/right hard partners (if any)

            if (ab <= lcv && ae >= rcv)
              { sidx = k;
                lidx = ridx = -1;
                ltop = rtop = 0;
                continue;
              }

            if (sidx >= 0)
              continue;
            
            if (ae >= rcv && ab > lft)
              { if (rtop < 10)
                  rstack[rtop++] = k;
#ifdef SHOW_PAIRS
                printf("r");
#endif
                if (ab <= rcv && ridx < 0)
                  { ridx = k;
#ifdef SHOW_PAIRS
                    printf("+");
#endif
                  }
              }

            if (ab <= lcv && ae < rgt)
              { if (ltop < 10)
                  lstack[ltop++] = k;
#ifdef SHOW_PAIRS
                printf("l");
#endif
                if (ae >= lcv)
                  { lidx = k;
#ifdef SHOW_PAIRS
                    printf("+");
#endif
                  }
	      }
          }

        //  Check for a hard contra pair and if found add
        //  Then check for a spanner and if so then add to gsort list.
        //  Then check for a spanning pair: use hard pair if available, otherwise
        //    use tightest pair and term it anciliary if endpoints are within an adjacent
        //    gap boundary, or soft otherwise.
        //  Finally, if left or right hard (but unpaired) then record as a conflict if
        //    projection extends MIN_OVERLAP past the other side.

        if (blast == bread)
          { if (ridx < 0)
              { if (lidx >= 0 && Ridx >= 0 && Lidx < 0)
                  { acnt += 1;
                    if (Ridx >= 0 && ovls[Ridx].path.abpos - ovls[Ridx].path.bbpos
                                             <= lft - MIN_OVERLAP)
                      rcnt -= 1;
#ifdef SHOW_PAIRS
                    printf(" = A");
#endif
                    continue;
                  }
              }
            else
              { if (lidx < 0 && Ridx < 0 && Lidx >= 0)
                  { acnt += 1;
                    if (Lidx >= 0 && ovls[Lidx].path.aepos + (blen-ovls[Lidx].path.bepos)
                                            >= rgt + MIN_OVERLAP)
                      lcnt -= 1;
#ifdef SHOW_PAIRS
                    printf(" = A");
#endif
                    continue;
                  }
              }
          }

        if (sidx >= 0)
          { gsort[gcnt].delta = DB->maxlen;
            gsort[gcnt].lidx  = sidx;
            gsort[gcnt].ridx  = sidx;
            gcnt += 1;
            scnt += 1;
#ifdef SHOW_PAIRS
            printf(" = S");
#endif
            continue;
          }

        if (ltop > 0 && rtop > 0)
          { int lok, rok, x;

            if (lidx < 0 || ridx < 0)
              { int dif, bst;
                int x, y;

                bst = 0x7fffffff;
                for (ltop--; ltop >= 0; ltop--)
                  { x = lstack[ltop];
                    for (rtop--; rtop >= 0; rtop--)
                      { y = rstack[rtop];
                        dif = (ovls[y].path.abpos - ovls[x].path.aepos)
                           - (ovls[y].path.bbpos - ovls[x].path.bepos);
                        dif = abs(dif);
                        if (dif < bst)
                          { bst  = dif;
                            lidx = x;
                            ridx = y;
                            dif = (ovls[ridx].path.abpos - ovls[lidx].path.aepos)
                                - (ovls[ridx].path.bbpos - ovls[lidx].path.bepos);
#ifdef SHOW_PAIRS
                            printf(" C(%d,%d = %d)",x,y,dif);
#endif
                          }
                      }
                  }
              }

            lok = 2;
            if (ovls[lidx].path.aepos < lcv)
              { x = ovls[lidx].path.aepos;
                lok = (lblock > FirstB && x <= lblock->beg && x >= lblock[-1].end - COVER_LEN);
              }
            rok = 2;
            if (ovls[ridx].path.abpos > rcv)
              { x = ovls[ridx].path.abpos;
                rok = (rblock < LastB && x >= rblock->end && x <= rblock[1].beg + COVER_LEN);
              }
  
            if (lok >= 2 && rok >= 2)
              gsort[gcnt].soft = 2;
            else if (lok >= 1 && rok >= 1)
              gsort[gcnt].soft = 1;
            else
              gsort[gcnt].soft = 0;
            gsort[gcnt].delta = (ovls[ridx].path.abpos - ovls[lidx].path.aepos)
                              - (ovls[ridx].path.bbpos - ovls[lidx].path.bepos);
            gsort[gcnt].lidx  = lidx;
            gsort[gcnt].ridx  = ridx;
            gcnt += 1;
#ifdef SHOW_PAIRS
            printf(" = G%d",gsort[gcnt-1].delta);
#endif
            continue;
          }

        if (ridx >= 0 && ovls[ridx].path.abpos - ovls[ridx].path.bbpos <= lft - MIN_OVERLAP)
          { rcnt += 1;
#ifdef SHOW_PAIRS
            printf(" = R");
#endif
          }
        if (lidx >= 0 && ovls[lidx].path.aepos + (blen-ovls[lidx].path.bepos) >= rgt + MIN_OVERLAP)
          { lcnt += 1;
#ifdef SHOW_PAIRS
            printf(" = L");
#endif
          }
      }

#ifdef SHOW_PAIRS
    printf("\n");
#endif

    { int ccnt, ocnt;

      if (lcnt < rcnt)
        ccnt = lcnt;
      else
        ccnt = rcnt;

      //  Analyze pair list gsort: if standard analysis (only hard pairs count) does not yield
      //    a span, then consider anciliary pair spanners (rarely makes a difference but does
      //    save a few.
  
      qsort(gsort,gcnt,sizeof(Spanner),GSORT);
      gcnt -= scnt;

      ocnt = gcnt;
      gcnt = analyze_gap_pairs(rgt-lft,gsort,ovls,gcnt,scnt,1,0);
      if (scnt < 4 && gcnt < 10 && gcnt < ccnt)
        {
#ifdef SHOW_PAIRS
          printf("  SPECIAL\n");
#endif
          gcnt = analyze_gap_pairs(rgt-lft,gsort,ovls,ocnt,scnt,0,1);
#ifdef SHOW_PAIRS
          if (gcnt >= 10 || gcnt >= ccnt)
            printf("  SWITCH\n");
#endif
        }
      else
        analyze_gap_pairs(rgt-lft,gsort,ovls,ocnt,scnt,1,1);

#ifdef DEBUG_GAP_STATUS
      printf("    lcnt = %d  gcnt = %d  scnt = %d  rcnt = %d acnt = %d\n",
             lcnt,gcnt-scnt,scnt,rcnt,acnt);
#endif

      //  Lots of contra pairs and less spanning support, call it an adaptamer gap.
  
      if (acnt >= .3*ccnt && gcnt < acnt) 
        {
#ifdef DEBUG_GAP_STATUS
          printf("    ADAPT %3d (vs %d)\n",acnt,gcnt);
#endif
          return (ADAPT);
        }

      //  If there is insufficient evidence for a span, then split.

      if (scnt < 4 && gcnt < 10 && gcnt <= ccnt)
        {
#ifdef DEBUG_GAP_STATUS
          if (ccnt >= 20)
            printf("    STRONG SPLIT\n");
          else
            printf("    WEAK SPLIT\n");
          if (gcnt >= 10)
            printf("  UNCERTAIN %5.1f %5d %3d\n",gcnt/(1.*ccnt),rgt-lft,gcnt);
#endif
          return (SPLIT);
        }

      //  Otherwise consider the gap spannable and try to find a viable patch, declaring a split
      //    iff all patch attemtps fail
  
      else
        { Patch    *can;
          int       ncand;
          uint8    *qa;
          Interval *clb, *crb;

          qa   = QV + QV_IDX[ovls[0].aread];
          clb  = lblock;
          crb  = rblock;

          //  First make sure enough partners provide anchors, and if not
          //     shift them back to the next good segment of A-read

          { int    nshort;

            nshort = 0;
            for (j = 0; j < gcnt; j++)
              { if (lft > ovls[gsort[j].lidx].path.aepos)
                  nshort += 1;
              }
            if (nshort > .2*gcnt)
              do
                { lft -= TRACE_SPACING;
                  if (lft <= clb->beg)
                    { if (clb <= FirstB)
                        break;
                      clb -= 1;
                      lft = clb->end;
                    }
                }
              while (qa[lft/TRACE_SPACING-1] > GOOD_QV);

            nshort = 0;
            for (j = 0; j < gcnt; j++)
              { if (rgt < ovls[gsort[j].ridx].path.abpos)
                  nshort += 1;
              }
            if (nshort > .2*gcnt)
              do
                { rgt += TRACE_SPACING;
                  if (rgt >= crb->end)
                    { if (crb >= LastB)
                        break;
                      crb += 1;
                      rgt = crb->beg;
                    }
                }
              while (qa[rgt/TRACE_SPACING] > GOOD_QV);

            //  Could not find primary anchor pair, then declare a SPLIT

            if (clb < FirstB || crb > LastB)
              {
#ifdef DEBUG_GAP_STATUS
                printf("    ANCHOR FAIL (BOUNDS)\n");
#endif
                return (SPLIT);
              }
          }

          //  Count all patch candidates that have a good anchor pair

          ncand  = 0;
          for (j = 0; j < gcnt; j++)
            { lidx = gsort[j].lidx;
              ridx = gsort[j].ridx;

#ifdef DEBUG_PATCHING
              if (lidx != ridx)
                printf("       %5d [%5d,%5d] [%5d,%5d] %4d",
                       ovls[lidx].bread,ovls[lidx].path.abpos,ovls[lidx].path.aepos,
                       ovls[ridx].path.abpos,ovls[ridx].path.aepos,gsort[j].delta);
              else
                printf("       %5d [%5d,%5d] SSS",
                       ovls[lidx].bread,ovls[lidx].path.abpos,ovls[lidx].path.aepos);
#endif

              can = compute_patch(lft,rgt,ovls+lidx,ovls+ridx);

              if (can != NULL)
                {
#ifdef DEBUG_PATCHING
                  printf(" %d",can->end-can->beg);
#endif
                  if (can->anc <= ANCHOR_THRESH)
                    { ncand += 1;
#ifdef DEBUG_PATCHING
                      printf("  AA %d %d(%d)",can->anc,can->bad,can->avg);
#endif
                    }
                }
#ifdef DEBUG_PATCHING
              printf("\n");
#endif
            } 

          //  If there are less than 5 of them, then seek better anchor points a bit
          //    further back

          if (ncand < 5)
            { int x, best, nlft, nrgt;
              int nanchor, ntry;

#ifdef DEBUG_PATCHING
              printf("     NOT ENOUGH\n");
#endif

              //  Try 4 additional anchor spots located at good intervals of A (if available)
              //     One can cross other gaps in the search.  Try the one with the most
              //     partners having match scores below the anchor threshold.  Do this to the
              //     left and right.  (A better search could be arranged (i.e. find smallest
              //     spanning pair of adjusted anchors, but this situation happens 1 in 5000
              //     times, so felt it was not worth it).

              ntry = 0;
              nlft = lft;
              best = -1;
              for (x = lft; ntry < 5; x -= TRACE_SPACING)
                if (x <= clb->beg)
                  { if (clb <= FirstB)
                      break;
                    clb -= 1;
                    x = clb->end + TRACE_SPACING;
                  }
                else if (qa[x/TRACE_SPACING-1] <= GOOD_QV)
                  { ntry += 1;
                    nanchor = 0;
                    for (j = 0; j < gcnt; j++)
                      if (eval_lft_anchor(x,ovls+gsort[j].lidx) <= ANCHOR_THRESH)
                        nanchor += 1;
#ifdef DEBUG_PATCHING
                    printf("       %5d: %3d\n",x,nanchor);
#endif
                    if (nanchor > best)
                      { best = nanchor;
                        nlft = x;
                      }
                  }
#ifdef DEBUG_PATCHING
              printf("          %5d->%5d\n",lft,nlft);
#endif

              ntry = 0;
              nrgt = rgt;
              best = -1;
              for (x = rgt; ntry < 5; x += TRACE_SPACING)
                if (x >= crb->end)
                  { if (crb >= LastB)
                      break;
                    crb += 1;
                    x = crb->beg - TRACE_SPACING;
                  }
                else if (qa[x/TRACE_SPACING] <= GOOD_QV)
                  { ntry += 1;
                    nanchor = 0;
                    for (j = 0; j < gcnt; j++)
                      if (eval_rgt_anchor(x,ovls+gsort[j].ridx) <= ANCHOR_THRESH)
                        nanchor += 1;
#ifdef DEBUG_PATCHING
                    printf("       %5d: %3d\n",x,nanchor);
#endif
                    if (nanchor > best)
                      { best = nanchor;
                        nrgt = x;
                      }
                  }
#ifdef DEBUG_PATCHING
              printf("          %5d->%5d\n",rgt,nrgt);
#endif

              //  If a better candidate pair of anchor points does not exist, then split.

              if (lft == nlft && rgt == nrgt)
                {
#ifdef DEBUG_GAP_STATUS
                  printf("    ANCHOR FAIL (ONCE) %d\n",ncand);
#endif
                  return (SPLIT);
                }
              lft = nlft;
              rgt = nrgt;

              //  Check out if the new anchor pair has 5 or more candidate patches

              ncand  = 0;
              for (j = 0; j < gcnt; j++)
                { lidx = gsort[j].lidx;
                  ridx = gsort[j].ridx;

#ifdef DEBUG_PATCHING
                  if (lidx != ridx)
                    printf("       %5d [%5d,%5d] [%5d,%5d] %4d",
                           ovls[lidx].bread,ovls[lidx].path.abpos,ovls[lidx].path.aepos,
                           ovls[ridx].path.abpos,ovls[ridx].path.aepos,gsort[j].delta);
                  else
                    printf("       %5d [%5d,%5d] SSS",
                           ovls[lidx].bread,ovls[lidx].path.abpos,ovls[lidx].path.aepos);
#endif

                  if (lft <= ovls[lidx].path.aepos && rgt >= ovls[ridx].path.abpos)
                    { can = compute_patch(lft,rgt,ovls+lidx,ovls+ridx);
    
                      if (can != NULL)
                        {
#ifdef DEBUG_PATCHING
                          printf(" %d",can->end-can->beg);
#endif
                          if (can->anc <= ANCHOR_THRESH)
                            { ncand += 1;
#ifdef DEBUG_PATCHING
                              printf("  AA %d %d(%d)",can->anc,can->bad,can->avg);
#endif
                            }
                        }
                    }
#ifdef DEBUG_PATCHING
                  printf("\n");
#endif
                } 

              //  Could not arrange 5 patch candidates, give up and split.

              if (ncand < 5)
                {
#ifdef DEBUG_GAP_STATUS
                  printf("    ANCHOR FAIL (TWICE) %d\n",ncand);
#endif
                  return (SPLIT);
                }
            }

          *p_lft = lft;
          *p_rgt = rgt;

#ifdef DEBUG_GAP_STATUS
          printf("    SPAN %5d:  PATCHABLE\n",rgt-lft);
#endif
          return (SPAN);
        }
    }
  }
}

static int *GAP_ANALYSIS(int aread, Overlap *ovls, int novl, Interval *block, int *p_nblk)
{ static int  bmax = 0;
  static int *status = NULL;      //  Status of gaps between HQ_blocks
 
#if defined(DEBUG_SUMMARY) || defined(DEBUG_CLASS)
  static char *status_string[4] = { "LOWQ", "SPAN", "SPLIT", "ADAPT" };
#endif

  int nblk;
  int i, j;
  int slft = 0, srgt = 0;

  (void) aread;

  nblk = *p_nblk;
  if (nblk > bmax)
    { bmax = 1.2*nblk + 100;
      status = (int *) Realloc(status,bmax*sizeof(int),"Allocating status vector");
      if (status == NULL)
        exit (1);
    }

  gap_status(ovls,novl,block,block+nblk,NULL,NULL);  //  Initialization call
  j = 0;
  for (i = 1; i < nblk; i++)
    { status[i] = gap_status(ovls,novl,block+j,block+i,&slft,&srgt);
      if (status[i] == SPAN)
        { while (slft < block[j].beg)
            j -= 1;
          block[j].end = slft;
          j += 1;
          status[j] = status[i];
          while (srgt > block[i].end)
            i += 1;
          block[j] = block[i];
          block[j].beg = srgt;
        }
      else
        { j += 1;
          status[j] = status[i];
          block[j]  = block[i];
        }
    }
  nblk = j+1;

#ifdef DEBUG_SUMMARY
#ifdef DEBUG_GAP_STATUS
  printf("  FINAL:\n");
#endif
  printf("    [%d,%d]",block[0].beg,block[0].end);
  for (i = 1; i < nblk; i++)
    printf(" %s [%d,%d]",status_string[status[i]],block[i].beg,block[i].end);
#endif

#ifdef DEBUG_CLASS
  for (i = 1; i < nblk; i++)
    printf("AREAD %d  %s  [%d,%d]\n",aread,status_string[status[i]],block[i-1].end,block[i].beg);
#endif
  
  *p_nblk  = nblk;
  return (status);
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

static void GAPS(int aread, Overlap *ovls, int novl)
{ int       alen;

  int       nblk;
  Interval *block;
  int      *status;

#if defined(DEBUG_HQ_BLOCKS) || defined(DEBUG_HOLE_FINDER) || defined(DEBUG_GAP_STATUS) || defined(DEBUG_SUMMARY)
  printf("\n");
  printf("AREAD %d\n",aread);
#endif

  alen = DB->reads[aread].rlen;
  if (alen < HGAP_MIN)
    {
#ifdef ANNOTATE
      fwrite(&HQ_INDEX,sizeof(int64),1,HQ_AFILE);
      fwrite(&SN_INDEX,sizeof(int64),1,SN_AFILE);
      fwrite(&SP_INDEX,sizeof(int64),1,SP_AFILE);
      fwrite(&AD_INDEX,sizeof(int64),1,AD_AFILE);

      fwrite(&HL_INDEX,sizeof(int64),1,HL_AFILE);

      fwrite(&KP_INDEX,sizeof(int64),1,KP_AFILE);
#endif

      fwrite(&TR_INDEX,sizeof(int64),1,TR_AFILE);
      return;
    }

  nreads += 1;
  totlen += alen;

  //  Partition into HQ-blocks

  block = HQ_BLOCKS(aread,&nblk);

  //  Find holes and modify HQ-blocks if necessary

  if (nblk > 0)
    nblk = FIND_HOLES(aread,ovls,novl,block,nblk);

  //  Determine the status of each gap between a pair of blocks

  if (nblk > 0)
    status = GAP_ANALYSIS(aread,ovls,novl,block,&nblk);

  //  No blocks? ==> nothing to do

  if (nblk <= 0)
    { nelim += 1;
      nelimbp += alen;
#ifdef ANNOTATE
      fwrite(&HQ_INDEX,sizeof(int64),1,HQ_AFILE);
      fwrite(&SN_INDEX,sizeof(int64),1,SN_AFILE);
      fwrite(&SP_INDEX,sizeof(int64),1,SP_AFILE);
      fwrite(&AD_INDEX,sizeof(int64),1,AD_AFILE);

      fwrite(&HL_INDEX,sizeof(int64),1,HL_AFILE);

      fwrite(&KP_INDEX,sizeof(int64),1,KP_AFILE);
#endif

      fwrite(&TR_INDEX,sizeof(int64),1,TR_AFILE);
      return;
    }

#ifdef ANNOTATE
  { int i;

    for (i = 0; i < nblk; i++)
      { fwrite(&(block[i].beg),sizeof(int),1,HQ_DFILE);
        fwrite(&(block[i].end),sizeof(int),1,HQ_DFILE);
        if (i > 0)
          { if (status[i] == SPAN || status[i] == LOWQ)
              { fwrite(&(block[i-1].end),sizeof(int),1,SN_DFILE);
                fwrite(&(block[i].beg),sizeof(int),1,SN_DFILE);
                SN_INDEX += 2*sizeof(int);
              }
            else if (status[i] == SPLIT)
              { fwrite(&(block[i-1].end),sizeof(int),1,SP_DFILE);
                fwrite(&(block[i].beg),sizeof(int),1,SP_DFILE);
                SP_INDEX += 2*sizeof(int);
              }
            else //  status[i] == ADAPT 
              { fwrite(&(block[i-1].end),sizeof(int),1,AD_DFILE);
                fwrite(&(block[i].beg),sizeof(int),1,AD_DFILE);
                AD_INDEX += 2*sizeof(int);
              }
          }
      }
    HQ_INDEX += 2*sizeof(int)*nblk;
    fwrite(&HQ_INDEX,sizeof(int64),1,HQ_AFILE);
    fwrite(&SN_INDEX,sizeof(int64),1,SN_AFILE);
    fwrite(&SP_INDEX,sizeof(int64),1,SP_AFILE);
    fwrite(&AD_INDEX,sizeof(int64),1,AD_AFILE);
  }
#endif

  //   Find largest non-adaptemer/subread range: block[abeg..aend)

  { int cmax, amax, abeg = 0, aend = 0;
    int p, i;

    amax = 0;
    p = 0;
    cmax = block[0].end-block[0].beg;
    for (i = 1; i < nblk; i++)
      if (status[i] == ADAPT)
        { if (cmax > amax)
            { amax = cmax;
              abeg = p;
              aend = i;
            }
          p = i;
          cmax = block[i].end - block[i].beg; 
        }
      else if (status[i] != SPLIT)
        cmax += block[i].end - block[i-1].end;
      else
        cmax += block[i].end - block[i].beg;
    if (cmax > amax)
      { amax = cmax;
        abeg = p;
        aend = nblk;
      }

    if (block[aend-1].end - block[aend-1].beg < TRACE_SPACING)
      { aend -= 1;
        nblk -= 1;  // assert: aend == nblk && status[aend-1] = SPLIT
      }
    if (block[aend-1].end == alen)
      block[aend-1].end = (alen/TRACE_SPACING)*TRACE_SPACING;
 
#ifdef DEBUG_SUMMARY
    printf("  :::  Keeping [%d,%d]\n",block[abeg].beg,block[aend-1].end);
#endif

  //  Accummulate statistics

    if (block[0].beg > 0)
      { n5trm += 1;
        n5trmbp += block[0].beg;
      }
    if (block[nblk-1].end < alen)
      { n3trm += 1;
        n3trmbp += alen - block[nblk-1].end;
      }
    if (abeg > 0)
      { natrm   += 1;
        natrmbp += block[abeg].beg - block[0].beg;
      }
    if (aend < nblk)
      { natrm   += 1;
        natrmbp += (block[nblk-1].end - block[aend-1].end);
      }
    for (i = abeg+1; i < aend; i++)
      { ngaps += 1;
        ngapsbp += block[i].beg - block[i-1].end;
        if (status[i] == LOWQ)
          { nlowq   += 1;
            nlowqbp += block[i].beg - block[i-1].end;
          }
        else if (status[i] == SPAN)
          { nspan   += 1;
            nspanbp += block[i].beg - block[i-1].end;
          }
        else // status[i] == SPLIT
          { nchim   += 1;
            nchimbp += block[i].beg - block[i-1].end;
          }
      }

#ifdef ANNOTATE
    fwrite(&(block[abeg].beg),sizeof(int),1,KP_DFILE);
    for (i = abeg+1; i < aend; i++)
      if (status[i] == SPLIT)
        { fwrite(&(block[i-1].end),sizeof(int),1,KP_DFILE);
          fwrite(&(block[i].beg),sizeof(int),1,KP_DFILE);
          KP_INDEX += 2*sizeof(int);
        }
    fwrite(&(block[aend-1].end),sizeof(int),1,KP_DFILE);
    KP_INDEX += 2*sizeof(int);
    fwrite(&KP_INDEX,sizeof(int64),1,KP_AFILE);
#endif

  //  Output .trim track for this read

    fwrite(&(block[abeg].beg),sizeof(int),1,TR_DFILE);
    fwrite(&(block[abeg].end),sizeof(int),1,TR_DFILE);
    TR_INDEX += 2*sizeof(int);
    for (i = abeg+1; i < aend; i++)
      { fwrite(status+i,sizeof(int),1,TR_DFILE);
        fwrite(&(block[i].beg),sizeof(int),1,TR_DFILE);
        fwrite(&(block[i].end),sizeof(int),1,TR_DFILE);
        TR_INDEX += 3*sizeof(int);
      }
    fwrite(&TR_INDEX,sizeof(int64),1,TR_AFILE);
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

  if (Read_Overlap(input,ovls) != 0)
    ovls[0].aread = INT32_MAX;
  else if (trace)
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

  if (ovls[n].aread < INT32_MAX)
    { fprintf(stderr,"%s: .las file overlaps don't correspond to reads in block %d of DB\n",
                     Prog_Name,DB_PART);
      exit (1);
    }

  return (max);
}

int main(int argc, char *argv[])
{ char       *root, *dpwd;
  int64       novl;
  DAZZ_TRACK *track;
  int         c;

  DAZZ_EXTRA  ex_hgap, ex_cest;
  DAZZ_EXTRA  ex_good, ex_bad, ex_trim;

  char *good_name = "Good QV threshold";
  char *bad_name  = "Bad QV threshold";
  char *trim_name = "Trimming statistics";
  int64 good64, bad64, tstats[18];

  //  Process arguments

  { int  i, j, k;
    int  flags[128];
    char *eptr;

    ARG_INIT("DAStrim")

    BAD_QV    = -1;
    GOOD_QV   = -1;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("v")
            break;
          case 'b':
            ARG_NON_NEGATIVE(BAD_QV,"Minimum QV score for being considered bad")
            break;
          case 'g':
            ARG_NON_NEGATIVE(GOOD_QV,"Maximum QV score for being considered good")
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
        fprintf(stderr,"      -g: Use as good qv threshold (and not auto-estimate).\n");
        fprintf(stderr,"      -b: Use as bad qv threshold (and not auto-estimate).\n");
        exit (1);
      }
  }

  //  Open trimmed DB and the qual-track

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

  //  Get .qual track and extras

  track = Load_Track(DB,"qual");
  if (track != NULL)
    { FILE       *afile;
      char       *aname;
      int         extra, tracklen, size;
      DAZZ_EXTRA  ex_qvs, ex_dif;

      QV_IDX = (int64 *) track->anno;
      QV     = (uint8 *) track->data;

      aname = Strdup(Catenate(DB->path,".","qual",".anno"),"Allocating anno file");
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

      COVERAGE = (int) ((int64 *) (ex_cest.value))[0];
      HGAP_MIN = (int) ((int64 *) (ex_hgap.value))[0];

      //  Compute -g and -b parameters

      { int64  qsum, qtotal;
        int64 *qgram;
        int    i, maxqv;
        int    gv, bv;

        qgram = (int64 *) (ex_qvs.value);
        maxqv = ex_qvs.nelem - 1;

        qtotal = 0;
        for (i = 0; i < maxqv; i++)
          qtotal += qgram[i];

        bv = gv = -1;
        qsum = 0;
        for (i = maxqv-1; i >= 0; i--)
          if (qgram[i] > 0)
            { qsum += qgram[i];
              if ((100.*qsum)/qtotal > 7. && bv < 0)
                bv = i+1;
              if ((100.*qsum)/qtotal > 20. && gv < 0)
                gv = i+1;
            }

        if (GOOD_QV < 0)
          GOOD_QV = gv;
        if (BAD_QV < 0)
          BAD_QV = bv;
      }
    }
  else
    { fprintf(stderr,"%s: Must have a 'qual' track, run DASqv\n",Prog_Name);
      exit (1);
    }

  if (GOOD_QV > BAD_QV)
    { fprintf(stderr,"%s: Good QV threshold (%d) > Bad QV threshold (%d) ?\n",
                     Prog_Name,GOOD_QV,BAD_QV);
      exit (1);
    }

  //  Setup extras

  ex_good.vtype = DB_INT;     // Good QV threshold
  ex_good.nelem = 1;
  ex_good.accum = DB_EXACT;
  ex_good.name  = good_name;
  good64 = GOOD_QV;
  ex_good.value = &good64;

  ex_bad.vtype = DB_INT;      // Bad QV threshold
  ex_bad.nelem = 1;
  ex_bad.accum = DB_EXACT;
  ex_bad.name  = bad_name;
  bad64 = BAD_QV;
  ex_bad.value = &bad64;

  ex_trim.vtype = DB_INT;     // Trim statistics
  ex_trim.nelem = 16;
  ex_trim.accum = DB_SUM;
  ex_trim.name  = trim_name;
  ex_trim.value = &tstats;

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

          //  Set up QV trimming track

#define SETUP(AFILE,DFILE,INDEX,anno,data,S)					\
{ int len, size;								\
										\
  if (DB_PART > 0)								\
    { AFILE = Fopen(Catenate(dpwd,PATHSEP,root,					\
                                Numbered_Suffix(".",DB_PART,anno)),"w");	\
      DFILE = Fopen(Catenate(dpwd,PATHSEP,root,					\
                                Numbered_Suffix(".",DB_PART,data)),"w");	\
    }										\
  else										\
    { AFILE = Fopen(Catenate(dpwd,PATHSEP,root,anno),"w");			\
      DFILE = Fopen(Catenate(dpwd,PATHSEP,root,data),"w");			\
    }										\
  if (AFILE == NULL || DFILE == NULL)						\
    exit (1);									\
										\
  len  = DB_LAST - DB_FIRST;							\
  size = S;									\
  fwrite(&len,sizeof(int),1,AFILE);						\
  fwrite(&size,sizeof(int),1,AFILE);						\
  INDEX = 0;									\
  fwrite(&INDEX,sizeof(int64),1,AFILE);						\
}

          SETUP(TR_AFILE,TR_DFILE,TR_INDEX,".trim.anno",".trim.data",8)

#ifdef ANNOTATE
          SETUP(HQ_AFILE,HQ_DFILE,HQ_INDEX,".hq.anno",".hq.data",0)
          SETUP(SN_AFILE,SN_DFILE,SN_INDEX,".span.anno",".span.data",0)
          SETUP(SP_AFILE,SP_DFILE,SP_INDEX,".split.anno",".split.data",0)
          SETUP(AD_AFILE,AD_DFILE,AD_INDEX,".adapt.anno",".adapt.data",0)

          SETUP(HL_AFILE,HL_DFILE,HL_INDEX,".hole.anno",".hole.data",0)

          SETUP(KP_AFILE,KP_DFILE,KP_INDEX,".keep.anno",".keep.data",0)
#endif

          //  Get trace point spacing information

          fread(&novl,sizeof(int64),1,input);
          fread(&TRACE_SPACING,sizeof(int),1,input);

          //  Initialize statistics gathering

          nreads  = 0;
          totlen  = 0;
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

          if (VERBOSE)
            { printf("\nDAStrim");
              if (HGAP_MIN  > 0)
                printf(" -H%d",HGAP_MIN);
              printf(" -c%d -g%d -b%d %s %s\n",
                     COVERAGE,GOOD_QV,BAD_QV,argv[1],Block_Arg_Root(parse));
            }

          //  Process each read pile

          make_a_pass(input,GAPS,1);

          //  Write out extras and close .trim track

          tstats[ 0] = nelim;
          tstats[ 1] = n5trm;
          tstats[ 2] = n3trm;
          tstats[ 3] = natrm;
          tstats[ 4] = nelimbp;
          tstats[ 5] = n5trmbp;
          tstats[ 6] = n3trmbp;
          tstats[ 7] = natrmbp;
          tstats[ 8] = ngaps;
          tstats[ 9] = nlowq;
          tstats[10] = nspan;
          tstats[11] = nchim;
          tstats[12] = ngapsbp;
          tstats[13] = nlowqbp;
          tstats[14] = nspanbp;
          tstats[15] = nchimbp;

          Write_Extra(TR_AFILE,&ex_hgap);
          Write_Extra(TR_AFILE,&ex_cest);
          Write_Extra(TR_AFILE,&ex_good);
          Write_Extra(TR_AFILE,&ex_bad);
          Write_Extra(TR_AFILE,&ex_trim);

          fclose(TR_AFILE);
          fclose(TR_DFILE);

#ifdef ANNOTATE
          fclose(HQ_AFILE);
          fclose(HQ_DFILE);

          fclose(SN_AFILE);
          fclose(SN_DFILE);

          fclose(SP_AFILE);
          fclose(SP_DFILE);

          fclose(AD_AFILE);
          fclose(AD_DFILE);

          fclose(HL_AFILE);
          fclose(HL_DFILE);

          fclose(KP_AFILE);
          fclose(KP_DFILE);
#endif

          fclose(input);

          //  If verbose output statistics summary to stdout

          if (VERBOSE)
            { printf("\n  Input:    ");
              Print_Number((int64) nreads,7,stdout);
              printf(" (100.0%%) reads     ");
              Print_Number(totlen,12,stdout);
              printf(" (100.0%%) bases");
              if (HGAP_MIN > 0)
                { printf(" (another ");
                  Print_Number((int64) ((DB_LAST-DB_FIRST)-nreads),0,stdout);
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

              printf("  Adapter:  ");
              Print_Number(natrm,7,stdout);
              printf(" (%5.1f%%) reads     ",(100.*natrm)/nreads);
              Print_Number(natrmbp,12,stdout);
              printf(" (%5.1f%%) bases\n",(100.*natrmbp)/totlen);

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
              Print_Number(n5trm+n3trm+nelim+natrm+nchim,7,stdout);
              printf(" clips              ");
              Print_Number(n5trmbp+n3trmbp+nelimbp+natrmbp,12,stdout);
              printf(" (%5.1f%%) bases\n",(100.*(n5trmbp+n3trmbp+nelimbp+natrmbp+nchimbp))/totlen);

              printf("  Patched:  ");
              Print_Number(nlowq+nspan,7,stdout);
              printf(" patches            ");
              Print_Number(nlowqbp+nspanbp,12,stdout);
              printf(" (%5.1f%%) bases\n",(100.*(nlowqbp+nspanbp))/totlen);
            }
        }

      Free_Block_Arg(parse);
    }

  free(dpwd);
  free(root);

  Close_DB(DB);
  free(Prog_Name);

  exit (0);
}
