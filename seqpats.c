/*******************************************************************************************
 *
 *  Library of bit-vector based sequence anlyzers for finding teriminal palindromes
 *    (unremoved adapter signal), palindromes (alignments across the main anti-diagonal),
 *    and satellites (alignments just off the main diagonal).
 *
 *  Author:  Gene Myers
 *  Date  :  March 2015
 *
 *******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <sys/file.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "DB.h"
#include "align.h"
#include "seqpats.h"

#undef   DEBUG_DP
#undef   DEBUG_CONFIRM
#undef   DEBUG_SEEDS
#undef   DEBUG_EXPLORE
#undef   SHOW
#undef   SHOW_UNITS

#define MAX_ERATE .35

#define SIGMA 4
#define WORD  uint64
#define W     64
#define ALL   0xffffffffffffffffll
#define WBIT  0x8000000000000000ll

#define FLENG      4
#define PLENG      6
#define SLENG      4

#define FLIPMIN  200    //   Minimum seed match length (must be <= FlipMax)
#define FGAPMAX 1000    //   Maximum permitted gap in flip palindrome

static int FlipMax = FLENG * W;
static int FlipThr = FLENG * W * MAX_ERATE;

static int PalOff = 32;
static int PalMax = PLENG * W;
static int PalThr = PLENG * (W/2) * MAX_ERATE;

#define MICROMAX  20   //  Anything <= 20bp is a micro, versus mini (use different detectors)
#define MICRORNG 256   //  Check 1..MICRORNG for micro spacing

static int SatMax = SLENG * W;
static int SatThr = W * MAX_ERATE;

static int Sat_Width[MICROMAX+1] =
    { -1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5 };

#define SAT_BORD  .25



/****************************************************************************************\
*                                                                                        *
*  BASIC BLOCK DP MACROS                                                                 *
*                                                                                        *
\****************************************************************************************/

#ifdef DEBUG

static int Letter[4] = { 'a', 'c', 'g', 't' };

#endif

#define ADVANCE_BLOCK(p,m,b,c)		\
{ WORD P, M, U, X, Y;			\
					\
  U = b;				\
  P = p;				\
  M = m;				\
					\
  if (c < 0)				\
    Y  = U | 1;				\
  else					\
    Y = U;				\
  X  = (((Y & P) + P) ^ P) | Y;		\
  U |= M;				\
					\
  Y = P;				\
  P = M | ~ (X | Y);			\
  M = Y & X;				\
    					\
  switch (c)				\
  { case 1:				\
      Y = (P << 1) | 1;			\
      p = (M << 1) | ~ (U | Y);		\
      break;				\
    case 0:				\
      Y = (P << 1);			\
      p = (M << 1) | ~ (U | Y);		\
      break;				\
    case -1:				\
       Y = (P << 1);			\
       p = (M << 1) | ~ (U | Y) | 1;	\
       break;				\
    }					\
  m = Y & U;				\
    					\
  if (P & WBIT)				\
    c = 1;				\
  else if (M & WBIT)			\
    c = -1;				\
  else					\
    c = 0;				\
}

#define BLOCK_DELTA(p,m,delta)		\
{ WORD P, M, U;				\
  int  o, b;				\
					\
  delta = 0;				\
  o = 0;				\
  P = p;				\
  M = m;				\
  U = WBIT;				\
  for (b = 0; b < W; b++)		\
    { if (P & U)			\
        { o -= 1;			\
          if (o < delta) delta = o;	\
        }				\
      else if (M & U)			\
        o += 1;				\
      U >>= 1;				\
    }					\
}

static void complement_read(char *read2, char *read1, int rlen)
{ int i, k;

  k = rlen-1;
  for (i = 0; i < rlen; i++)
    read2[i] = 3-read1[k--];
  read2[rlen] = read2[-1] = 4;
}

//   Setup the match vectors for first len symbols of read

void setup_search(char *read, int len, WORD **tran, WORD *_tran)
{ WORD *b, bvc, one;
  int   a, p, i, k;

  b = _tran;
  for (a = 0; a < SIGMA; a++)
    { tran[a] = b;
      for (p = 0; p < len; p += W)
        { bvc = 0;
          one = 1;
          k   = p+W;
          for (i = p; i < k; i++)
            { if (a == read[i])
                bvc |= one;
              one <<= 1;
            }
          *b++ = bvc;
        }
    }

#ifdef DEBUG_DP
  for (a = 0; a < SIGMA; a++)
    { printf("%3c ",Letter[a]);
      b = tran[a];
      for (p = 0; p < len; p += W)
        { for (i = 0; i < W; i++)
            if (*b & (0x1ll<<i))
              printf("1");
            else
              printf("0");
          printf(" ");
          b += 1;
        }
      printf("\n");
    }
#endif
}


/****************************************************************************************\
*                                                                                        *
*  POLYMERASE FLIP FINDER                                                                *
*                                                                                        *
\****************************************************************************************/

static int flip_confirm(Alignment *align, Work_Data *work, Align_Spec *spec,
                        int x, int y, Flip_Hit *hit)
{ Path *path;
  int   ab, bb, ae, be;
  int   qvm, qvg, gap, rlen;

  Local_Alignment(align,work,spec,x-y,x-y,x+y,-1,-1);

  path = align->path;
  rlen = align->alen;
  ab = path->abpos;
  bb = path->bbpos;
  ae = path->aepos;
  be = path->bepos;

  if (bb == 0 && ab != ae)

    { qvm = (100*path->diffs)/(ae-ab);
      if (qvm >= 100.*MAX_ERATE)
        return (0);

#ifdef DEBUG_CONFIRM
      printf("    (%d,%d) -> (%d,%d)",ab,bb,ae,be);
#endif

      gap = rlen-(ae+be);
      if (gap <= 0)
        { uint16 *trace = (uint16 *) path->trace;
          int     pt; 

          pt = 1;
          be = bb + trace[pt];
          ae = ((ab/100)+1)*100;
          while (ae+be < rlen)
            { pt += 2;
              bb  = be;
              be += trace[pt];
              ab  = ae;
              ae += 100;
              if (ae > path->aepos)
                ae = path->aepos;
            }
          hit->divpt = (rlen+(ae-be))/2;
          hit->qvmat = qvm;
          hit->gap   = 0;
#ifdef SHOW
          printf(" %5d: %5d %3d%% +++\n",hit->divpt,(path->aepos-path->abpos)/2,qvm);
#endif
          return (1);
        }
      else if (gap <= FGAPMAX && gap < ae-ab)
        { path->abpos = ae;
          path->bbpos = be;
          path->aepos = rlen-be;
          path->bepos = rlen-ae;
          path->diffs = gap;
          Compute_Trace_ALL(align,work);
          qvg = (100*path->diffs)/gap;

          hit->divpt = (rlen+(ae-be))/2;
          hit->qvmat = qvm;
          hit->gap   = gap;
          hit->qvgap = qvg;
#ifdef SHOW
          printf(" %5d: %5d %3d%%  %5d %3d%% +\n",hit->divpt,ae-ab,qvm,gap/2,qvg);
#endif

#ifdef DEBUG_CONFIRM
          printf(" [%d] (%d,%d) -> (%d,%d)",gap,rlen-be,rlen-ae,rlen-bb,rlen-ab);
          printf("\n");
          Print_Alignment(stdout,align,work,4,100,10,0,5);
          printf("\n");
#endif
          return (1);
        }
#ifdef SHOW
      else
        printf(" %5d %3d%%  %5d --\n",ae-ab,qvm,gap/2);
#endif
    }
#ifdef SHOW
  else
    printf("  No match ---\n");
#endif

  return (0);
}


//   Find potential flip middle + gap length, return non-zero if find one

typedef struct
  { WORD P;
    WORD M;
    int  V;
  } Scell;

static int flip_filter(Alignment *align, Work_Data *work, Align_Spec *spec, Flip_Hit *hit)
{ WORD  *e, *tran[SIGMA], _tran[SIGMA*FLENG];
  Scell *s, S[FLENG], *SE = S+(FLENG-1);
  int    besti, bestv;
  int    i, c;

  int    rlen = align->alen;
  char  *aseq = align->aseq;
  char  *bseq = align->bseq;

  setup_search(bseq,FlipMax,tran,_tran);

  for (s = S; s <= SE; s++)
    { s->P = ALL;
      s->M = 0;
      s->V = ((s-S)+1)*W;
    }

  // for each column of the dp matrix

  besti = 0;
  bestv = FlipThr;
  for (i = 0; i < rlen; i++)
    { e  = tran[(int) (aseq[i])];

#ifdef DEBUG_DP
      printf("  %c: ",Letter[(int) aseq[i]]);
      printf("Col =");
#endif

      // compute the first FLENG blocks

      c = 0;
      for (s = S; s <= SE; s++)
        { ADVANCE_BLOCK(s->P,s->M,*e++,c)
          s->V += c;
#ifdef DEBUG_DP
          { int k;
            for (k = 0; k < W; k++)
              if (s->P & (((long) 1) << k))
                printf("+");
              else if (s->M & (((long) 1) << k))
                printf("-");
              else
                printf("0");
            printf(" [%d]",s->V);
          }
#endif
        }
#ifdef DEBUG_DP
      printf("\n");
#endif

      // for each interval of columns scoring less than MAX_ERATE error
      //   on row FlipMax, call flip_confirm on the best one, and if get
      //   a hit return with said.

      if (SE->V < FlipThr)
        { if (SE->V < bestv)
            { besti = i+1;
              bestv = SE->V;
            }
        }
      else
        { if (bestv < FlipThr)
            { if (flip_confirm(align,work,spec,besti,FlipMax,hit))
                return (1);
              bestv = FlipThr;
            }
        }
    }

  // for the last column, walk down the boundary from row FlipMax
  //   calling flip_confim on any potential seeds.  Tricky part is
  //   the below-threshold intervals along the row: need to wrap
  //   around the boundary to follow the last column.

  { WORD   p, m, u;
    double f, bestf;
    int    j, bestj;

    if (bestv >= FlipThr)
      bestf = MAX_ERATE;
    else
      { besti = rlen;
        bestj = FlipMax;
        bestf = (1.*bestv)/FlipMax;
      } 

    c = SE->V;
    j = FlipMax;
    for (s = SE; s > S && j >= FLIPMIN; s--)
      { p = s->P;
        m = s->M;
        u = WBIT;
        for (i = 0; i < W && j >= FLIPMIN; i++)

          { f = (1.*c)/j;
            if (f < MAX_ERATE)
              { if (f < bestf)
                  { besti = rlen;
                    bestj = j;
                    bestf = f;
                  }
              }
            else
              { if (bestf < MAX_ERATE)
                  { if (flip_confirm(align,work,spec,besti,bestj,hit))
                      return (1);
                    bestf = MAX_ERATE;
                  }
              }

            if (p & u)
              c -= 1;
            else if (m & u)
              c += 1;
            u >>= 1;
            j -= 1;
          }
      }

    if (bestf < MAX_ERATE)
      { if (flip_confirm(align,work,spec,besti,bestj,hit))
          return (1);
      }
  }

  return (0);
}

Flip_Hit *Flip_Finder(HITS_DB *db, int i, int *nhit, int beg, int end)
{ static Flip_Hit   *hit = NULL;
  static int         hmax = 0;
  static Work_Data  *work;
  static Align_Spec *spec;
  static char       *read1, *read2;
  static int         firstime = 1;

  Alignment  _align, *align = &_align;
  Path       _path,  *path  = &_path;

  int rlen, clen, div, off;
  int phit, shit;

  if (firstime)
    { firstime = 0;

      spec = New_Align_Spec( .70, 100, db->freq);
      work = New_Work_Data();

      read1 = New_Read_Buffer(db);
      read2 = New_Read_Buffer(db);
    }

  align->path = path;
  clen = db->reads[i].rlen;           //  Fetch read

#ifdef SHOW
  printf("Analyzing %d(%x) %d\n",i,(db->reads[i].flags & DB_BEST) != 0,clen);
  fflush(stdout);
#endif

  if (clen < FlipMax)
    { *nhit = 0;
      return (hit);
    }

  if (Load_Read(db,i,read1,0))
    exit (1);
  complement_read(read2,read1,clen);

  if (hmax == 0)
    { hmax = 10;
      hit  = Malloc(sizeof(Flip_Hit)*hmax,"Allocating flip vector");
      if (hit == NULL)
        exit (1);
    }

  align->aseq = read2 + (clen-end);
  align->bseq = read1 + beg;
  align->alen = rlen = (end-beg);
  align->blen = rlen;
  off = clen-end;
  align->aseq[rlen] = align->bseq[rlen] = 4;
  align->aseq[-1]   = align->bseq[-1]   = 4;

  phit = 0;
  while (flip_filter(align,work,spec,hit+phit))
    { div   = hit[phit].divpt;
      off  += div;
      rlen -= div;
      hit[phit].divpt = clen - off;
      align->alen = align->blen = rlen;
      align->aseq = read2 + off;
      align->aseq[rlen] = align->bseq[rlen] = 4;
      align->aseq[-1]   = align->bseq[-1]   = 4;

      phit += 1;
      if (phit >= hmax)
        { hmax = 1.2*phit + 10;
          hit  = Realloc(hit,sizeof(Flip_Hit)*hmax,"Allocating flip vector");
          if (hit == NULL)
            exit (1);
        }
    }

  shit = phit;

  clen = db->reads[i].rlen;           //  Start over
  if (Load_Read(db,i,read1,0))
    exit (1);
  complement_read(read2,read1,clen);

  if (phit > 0)
    beg = hit[0].divpt;

  align->aseq = read1 + beg;
  align->bseq = read2 + (clen-end);
  align->alen = rlen = (end-beg);
  align->blen = rlen;
  off = beg;
  align->aseq[rlen] = align->bseq[rlen] = 4;
  align->aseq[-1]   = align->bseq[-1]   = 4;

  while (flip_filter(align,work,spec,hit+shit))
    { div   = hit[shit].divpt;
      off  += div;
      rlen -= div;
      hit[shit].divpt = off;
      align->alen = align->blen = rlen;
      align->aseq = read1 + off;
      align->aseq[rlen] = align->bseq[rlen] = 4;
      align->aseq[-1]   = align->bseq[-1]   = 4;

      shit += 1;
      if (shit >= hmax)
        { hmax = 1.2*shit + 10;
          hit  = Realloc(hit,sizeof(Flip_Hit)*hmax,"Allocating flip vector");
          if (hit == NULL)
            exit (1);
        }
    }

  { int      i, j;
    Flip_Hit x;

    j = phit-1;
    for (i = 0; i < j; i++, j--)
      { x = hit[i];
        hit[i] = hit[j];
        hit[j] = x;
      }
  }

  *nhit = shit;
  return (hit);
}


/****************************************************************************************\
*                                                                                        *
*  PAPLINDROME FINDER                                                                    *
*                                                                                        *
\****************************************************************************************/

static int pal_confirm(Alignment *align, Work_Data *work, Align_Spec *spec, int x, Pal_Hit *hit)
{ Path *path;
  int   ab, ae, bb, be, rlen;
  int   y;

#ifdef DEBUG_CONFIRM
  printf("  SEED %d %d\n",x,align->alen-(x+PalOff));
#endif

  rlen = align->alen;

  y = rlen-(x+PalOff);
  Local_Alignment(align,work,spec,x-y,x-y,x+y,-1,-1);

  path = align->path;
  ab = path->abpos;
  ae = path->aepos;
  bb = path->bbpos;
  be = path->bepos;

  if (ae-ab < W) return (0);

  Compute_Trace_PTS(align,work,100,GREEDIEST);

#ifdef DEBUG_CONFIRM
  printf("  ALIGN (%d,%d) - (%d,%d)",ab,bb,ae,be);
  printf("\n");
  Print_Alignment(stdout,align,work,4,100,10,0,5);
  printf("\n");
#endif

  if (path->diffs > MAX_ERATE * (ae-ab)) return (0);

  hit->qvmat = (100 * path->diffs) / (ae-ab);
  if (ae+be <= rlen)
    { hit->gap   = rlen - (ae+be);
#ifdef DEBUG_CONFIRM
      printf("  SHORT -> (%d,%d)",rlen-be,rlen-ae);
#endif
      ae = rlen - bb;
      be = rlen - ab;
#ifdef DEBUG_CONFIRM
      printf(" - (%d,%d)\n",ae,be);
#endif
    }
  else
    hit->gap = 0;

  if (ae-ab >= 2*W)
    { hit->midpt = (ab+ae)/2;
      hit->width = ae-ab;
      return (1);
    }

  return (0);
}

static int pal_filter(Alignment *align, Work_Data *work, Align_Spec *spec, Pal_Hit *hit)
{ WORD  *e, *tran[SIGMA], _tran[SIGMA*PLENG];
  Scell *s, S[PLENG], *SE = S+(PLENG-1);
  int    i, bestv, besti;
  char  *b;

  int    rlen = align->alen;
  char  *aseq = align->aseq;
  char  *bseq = align->bseq;

  bestv = PalThr;

  for (s = S; s <= SE; s++)
    { s->P = ALL;
      s->M = 0;
      s->V = ((s-S)+1)*W;
    }

  b = (bseq + rlen) - (PalMax + PalOff);
  setup_search(b,PalMax,tran,_tran);

  for (i = 0; i < rlen-PalOff; i++)
    { int  a, x, c;

      b -= 1;
 
      if (b < bseq)
        x = SIGMA;
      else
        x = *b;
      for (a = 0; a < SIGMA; a++)
        { WORD *t;
          int   g, h;
 
          h = (x == a);
          t = tran[a];
          for (s = S; s < SE; s++, t++)
            { g = ((*t & WBIT) != 0);
              *t <<= 1;
              if (h)
                *t |= 1;
              h = g;
            }
          *t <<= 1;
          if (h)
            *t |= 1;
        }

      c = -1;
      for (s = S; s <= SE; s++)
        { int d;

          if (s->P & WBIT)
            d = 1;
          else if (s->M & WBIT)
            d = -1;
          else
            d = 0;
          s->V -= d;
          s->P <<= 1;
          s->M <<= 1;
          if (c > 0)
            s->P |= 1;
          else if (c < 0)
            s->M |= 1;
          c = d;
        }

#ifdef DEBUG_DP
      { int k;
        printf("  ");
        for (s = S; s <= SE; s++)
          { for (k = 0; k < W; k++)
              if (s->P & (((long) 1) << k))
                printf("+");
              else if (s->M & (((long) 1) << k))
                printf("-");
              else
                printf("0");
            printf(" ");
          }
        printf("[%d]\n",SE->V);
      }
#endif

      e = tran[(int) (aseq[i])];
      c = -1;
      for (s = S; s <= SE; s++)
        { ADVANCE_BLOCK(s->P,s->M,*e++,c)
          s->V += c;
        }

#ifdef DEBUG_DP
      { int k;
      
        printf("%c ",Letter[(int) aseq[i]]);
        e = tran[(int) (aseq[i])];
        for (s = S; s <= SE; s++, e++)
          { for (k = 0; k < W; k++)
              if (*e & (0x1ll<<k))
                printf("1");
              else
                printf("0");
            printf(" ");
          }
        printf("\n ");
        for (s = S; s <= SE; s++)
          { for (k = 0; k < W; k++)
              if (s->P & (((long) 1) << k))
                printf("+");
              else if (s->M & (((long) 1) << k))
                printf("-");
              else
                printf("0");
            printf(" ");
          }
	printf("[%d]\n",SE->V);
      }
#endif

      // for each interval of columns scoring less than MAX_ERATE error
      //   on anti-diagonal rlen-PalOff, call pal on the best one, and if get
      //   a hit return with said.

      if (SE->V < PalThr)
        { if (SE->V < bestv)
            { besti = i+1;
              bestv = SE->V;
            }
        }
      else
        { if (bestv < PalThr)
            { if (pal_confirm(align,work,spec,besti,hit))
                return (1);
              bestv = PalThr;
            }
        }
    }

  return (0);
}


Pal_Hit *Pal_Finder(HITS_DB *db, int i)
{ static Pal_Hit    _hit, *hit = &_hit;
  static Alignment  _align, *align = &_align;
  static Path       _path,  *path  = &_path;
  static Work_Data  *work;
  static Align_Spec *spec;
  static char       *read1, *read2;
  static int         firstime = 1;
  static int         lasti = -1;

  int rlen, off, sym;

  if (firstime)
    { firstime = 0;

      spec = New_Align_Spec( .70, 100, db->freq);
      work = New_Work_Data();
      align->path = path;

      read1 = New_Read_Buffer(db);
      read2 = New_Read_Buffer(db);
    }

  rlen = db->reads[i].rlen;

  if (lasti == i)
    off = hit->midpt + (hit->width - PalMax)/2;
  else
    off = 0;

#ifdef SHOW
  if (lasti != i)
    { printf("Analyzing %d(%x) %d %d\n",i,(db->reads[i].flags & DB_BEST) != 0,off,rlen);
      fflush(stdout);
    }
#endif

  rlen -= off;
  if (rlen < W)
    { lasti = -1;
      return (NULL);
    }

  if (lasti != i)
    { if (Load_Read(db,i,read1,0))
        exit (1);
      complement_read(read2,read1,rlen);
    }

  sym = read1[off-1];
  read1[off-1] = 4;

  align->aseq = read1+off;
  align->bseq = read2;
  align->alen = rlen;
  align->blen = rlen;
  if (pal_filter(align,work,spec,hit))
    { hit->midpt += off;
      lasti = i;
#ifdef SHOW
      printf("  HIT %5d %5d %3d%% : %3d\n",hit->midpt,hit->width,hit->qvmat,hit->gap);
#endif
      read1[off-1] = sym;
      return (hit);
    }

  read1[off-1] = sym;
  lasti = -1;
  return (NULL);
}


/****************************************************************************************\
*                                                                                        *
*  SATELLITE REPEAT FINDER                                                               *
*                                                                                        *
\****************************************************************************************/

static int sat_period(Alignment *align)
{ int  *trace, tlen;
  char *a, *b;
  int   i, j, k;
  int   c, p, n;
  int   num, den;

  trace = (int *) align->path->trace;
  tlen  = align->path->tlen;

  a = align->aseq-1;
  b = align->bseq-1;
  i = n = align->path->abpos;
  j = align->path->bbpos;
  k = i-j;
  num = den = 0;

#ifdef SHOW_UNITS
  printf("(%d,%d) = %d\n",i,j,k);
#endif
  for (c = 0; c < tlen; c++)
    if ((p = trace[c]) < 0)
      { p = -p;
        while (i < p)
          { if (j == n)
              { num += k;
                den += 1;
#ifdef SHOW_UNITS
                printf("(%d,%d) = %d\n",i,j,k);
#endif
                n = i;
              }
            i += 1;
            j += 1;
          }
        if (j == n)
          { num += k;
            den += 1;
#ifdef SHOW_UNITS
            printf("(%d,%d) = %d\n",i,j,k);
#endif
            n = i;
          }
        j += 1;
        k -= 1;
      }
    else
      { while (j < p)
          { if (j == n)
              { num += k;
                den += 1;
#ifdef SHOW_UNITS
                printf("(%d,%d) = %d\n",i,j,k);
#endif
                n = i;
              }
            i += 1;
            j += 1;
          }
        i += 1;
        k += 1;
      }
  p = align->path->aepos;
  while (i < p)
    { if (j == n)
        { num += k;
          den += 1;
#ifdef SHOW_UNITS
          printf("(%d,%d) = %d\n",i,j,k);
#endif
          n = i;
        }
      i += 1;
      j += 1;
    }
  num += k;
  den += 1;
#ifdef SHOW_UNITS
  printf("(%d,%d) = %d\n",i,j,k);
  printf("Final: %d/%d = %d\n",num,den,num/den);
#endif
  return (num/den);
}

#ifdef DEBUG_SEEDS

static void print_vector(int *val, int vlen)
{ int i, j, max, min;

  min =  INT32_MAX;
  max = -INT32_MAX;
  for (i = 1; i <= vlen; i++)
    { if (val[i] < min)
        min = val[i];
      if (val[i] > max)
        max = val[i];
    }
  for (j = max; j >= min; j--)
    { printf("    %2d: ",j);
      for (i = 1; i <= vlen; i++)
        if (val[i] == j)
          printf("*");
        else
          printf(" ");
      printf("\n");
    }
}

#endif

static void sort_vector(int *val, int *sort, int vlen)
{ int i, p, q;
  int cnt[W+1];

  for (i = 0; i <= W; i++)
    cnt[i] = 0;
  for (i = 1; i <= vlen; i++)
    cnt[val[i]] += 1; 
  p = 0;
  for (i = 0; i <= W; i++)
    { q = p+cnt[i];
      cnt[i] = p;
      p = q;
    }
  for (i = vlen; i >= 1; i--)
    sort[cnt[val[i]]++] = i;
}

typedef struct
  { int    left;
    int    rght;
    int    chil;
    int    sibs;
    int    min;
    int    val;
    double score;
  } Level;

static int level_tree(int *val, int *sort, Level *tree, int vlen)
{ int   i, j;
  int   l, r, v, m;
  int   nl, nr, c, nset;
  int   in[SatMax+2];

  nset = 0;
  for (i = 0; i <= vlen+1; i++)
    in[i] = -1;

  for (i = 0; i < vlen; i++)
    { l = r = sort[i];
      v = val[r];
      for (j = i+1; j < vlen && sort[j] == l-1 && val[l-1] == v; j++)
        l = l-1;
      i += (r-l);

      nl = in[l-1];
      nr = in[r+1];
      if (nl >= 0)
        if (nr >= 0)
          { l = tree[nl].left;
            r = tree[nr].rght;
            if (tree[nr].val < v)
              { tree[c = nset++] = tree[nr];
                in[r] = c;
                tree[c].chil = nr;
                tree[nr].sibs = -1;
              }
            else
              c = nr;
            in[tree[c].left = l] = c;
            if (tree[c].min > tree[nl].min)
              tree[c].min = tree[nl].min;
            tree[nl].sibs = tree[c].chil;
            tree[c].chil = nl;
          }
        else
          { l = tree[nl].left;
            if (tree[nl].val < v)
              { tree[c = nset++] = tree[nl];
                in[l] = c;
                tree[c].chil = nl;
                tree[nl].sibs = -1;
              }
            else
              c = nl;
            in[tree[c].rght = r] = c;
          }
      else
        if (nr >= 0)
          { r = tree[nr].rght;
            if (tree[nr].val < v)
              { tree[c = nset++] = tree[nr];
                in[r] = c;
                tree[c].chil = nr;
                tree[nr].sibs = -1;
              }
            else
              c = nr;
            in[tree[c].left = l] = c;
          }
        else
          { c = nset++;
            tree[c].left = l;
            tree[c].rght = r;
            tree[c].chil = -1;
            tree[c].min  = v;
            in[l] = in[r] = c;
          }
      tree[c].val  = v;
      if (l == 1)
        l -= (v-val[l]);
      if (r == vlen)
        r += (v-val[r]);
      m = tree[c].min;
      tree[c].score = 2. * ((v-m)+1.) / ((r-l)+2.);
    }

  return (nset-1);
}

#ifdef DEBUG_SEEDS

static void print_tree(int s, Level *set, int deep)
{ printf("%*s",deep,"");
  printf(" [%3d,%3d] %.2f %d %d\n",set[s].left,set[s].rght,set[s].score,set[s].val,set[s].min);
  for (s = set[s].chil; s >= 0; s = set[s].sibs)
    print_tree(s,set,deep+2);
}

#endif

typedef struct
  { int row;
    int lcol, hcol;
    int val;
  } Seed;

static Seed *append_micro(int s, int *val, Level *tree, Seed *cand, int row, int vlen)
{ int j;

  if (tree[s].left <= MICROMAX)
    { if (tree[s].chil < 0)
        { if (tree[s].val <= SatThr)
            for (j = tree[s].left; j <= tree[s].rght && j <= MICROMAX; j++)
              { cand->hcol = j;
                cand->lcol = j;
                cand->row  = row;
                cand->val  = tree[s].min;
                cand += 1;
              }
        }
      else
        for (s = tree[s].chil; s >= 0; s = tree[s].sibs)
          cand = append_micro(s,val,tree,cand,row,vlen);
    }
  return (cand);
}

static Seed *append_minis(int s, int *val, Level *tree, Seed *cand, int row, int vlen)
{ int k, h, v;

  if (tree[s].score >= .5)
    { if (tree[s].val - tree[s].min >= 8 && tree[s].min <= SatThr)
        { cand->lcol = SatMax;
          cand->val  = tree[s].min;
          cand->row  = row;
          for (k = tree[s].left; k <= tree[s].rght; k = h)
            { v = val[k];
              for (h = k+1; h <= vlen; h++)
                if (val[h] != v)
                  break;
              if (v < val[k-1] && v < val[h] && v < tree[s].min+3)
                { if (cand->lcol > k)
                    cand->lcol = k;
                  cand->hcol = h-1;
                }
            }
          if (cand->hcol > MICROMAX)
            { if (cand->lcol <= MICROMAX)
                cand->lcol = MICROMAX+1;
              cand += 1;
            }
        }
    }
  else
    for (s = tree[s].chil; s >= 0; s = tree[s].sibs)
      cand = append_minis(s,val,tree,cand,row,vlen);
  return (cand);
}

#ifdef DEBUG_SEEDS

static void print_seeds(Seed *cand, int cbot, int ctop)
{ int   c, j, k;

  if (cbot == ctop) return;

  j = 1;
  printf("        ");
  for (c = cbot; c < ctop; c++)
    { k = cand[c].lcol;
      for ( ; j < k; j++)
        printf(" ");
      k = cand[c].hcol;
      printf("+");
      j += 1;
      for ( ; j <= k; j++)
        printf("-");
    }
  printf("\n");

  for (c = cbot; c < ctop; c++)
    printf("  (%4d,%3d-%3d) %d\n",cand[c].row,cand[c].lcol,cand[c].hcol,cand[c].val);
}

#endif

static void sort_seeds(Seed *cand, int ctop, int *sort)
{ int c, i, p, q;
  int cnt[W+1];

  for (i = 0; i <= W; i++)
    cnt[i] = 0;
  for (c = 0; c < ctop; c++)
    cnt[cand[c].val] += 1; 
  p = 0;
  for (i = 0; i <= W; i++)
    { q = p+cnt[i];
      cnt[i] = p;
      p = q;
    }
  for (c = 0; c < ctop; c++)
    sort[cnt[cand[c].val]++] = c;
}
 
static int SATSORT(const void *x, const void *y)
{ Sat_Hit *l = (Sat_Hit *) x;
  Sat_Hit *r = (Sat_Hit *) y;

  if (l->beg != r->beg)
    return (l->beg - r->beg);
  else
    return (l->end - r->end);
}

static Seed *merge_micro(int *val, Seed *cand, Seed *ctop, int row, int vlen)
{ int   hist[SatMax+2];
  int   i, j, v, m, tot;
  Seed *c, *n;

  if (vlen < MICRORNG)
    return (cand);

  for (i = 0; i <= MICROMAX+2; i++)
    hist[i] = 0;

  tot = m = 0;
  for (i = 1; i <= MICRORNG; i = j)
    { v = val[i];
      for (j = i+1; j <= vlen; j++)
        if (val[j] != v)
          break;
      if (v < val[i-1] && v < val[j])
        { v = (i+j)/2;
          if (v-m < MICROMAX+2)
            hist[v-m] += 1;
          tot += 1;
          m = v;
        }
    }
  
  m = .4 * tot;
  c = n = cand;
  for (i = 1; i <= MICROMAX; i++)
    if (c < ctop && c->lcol == i)
      { n->hcol = i;
        c += 1;
        n += 1;
      }
    else if (hist[i] > 0)
      if ((hist[i] + hist[i-1] >= m || hist[i] + hist[i+1] >= m) && val[i] <= SatThr)
        { n->hcol = i;
          n += 1;
        }

  for (c = cand; c < n; c++)
    { c->lcol = c->hcol;
      c->val  = val[c->lcol];
      c->row  = row;
    }

#ifdef DEBUG_SEEDS

  printf("\n");
  for (i = MICROMAX+1; i >= 1; i--)
    if (hist[i] > 0)
      { printf("  %3d: %3d",i,hist[i]);
        if ((hist[i] + hist[i-1] >= m || hist[i] + hist[i+1] >= m) && val[i] <= SatThr)
          printf(" XXX");
        printf("\n");
      }
  printf("     %3d (%d)\n",tot,m);

#endif

  return (n);
}

static Sat_Hit *sat_filter(Alignment *align, Align_Spec *spec, Work_Data *work)
{ static Seed *Cand = NULL;
  static int   Cmax = 0;
  int          ctop;

  static Sat_Hit *Match = NULL;
  static int      Mmax  = 0;
  int             mtop;

  int   Val[SatMax+2];
  WORD *_tran[SIGMA], tran[SIGMA];

  char *aseq = align->aseq;
  int   rlen = align->alen;

  int   lrow, crow, hit;

  Val[0] = Val[SatMax+1] = W+1;

  ctop = 0;
  for (lrow = 0; lrow < rlen-W; lrow += W)
    { crow = lrow+W;

      setup_search(aseq+lrow,W,_tran,tran);

      { int    i, t;
        WORD   e, d, p, m;
        int    v, c, *val;

        p = ALL;
        m = 0;
        v = W;
        d = 1;
        for (i = lrow+1; i <= crow; i++)
          { e = tran[(int) (aseq[i])] & d;
            d = (d << 1) | 1;
            c = 0;
            ADVANCE_BLOCK(p,m,e,c)
            v += c;

#ifdef DEBUG_DP
            { int k, h;

              printf("  %c:",Letter[(int) aseq[i]]);
              for (k = 0; k < W; k++)
                { if (e & (0x1ll << k))
                    printf(" +");
                  else
                    printf("  ");
                  printf("%c",Letter[(int) aseq[r+k]]);
                }
              printf("\n    ");
              h = 0;
              for (k = 0; k < W; k++)
                { if (p & (((long) 1) << k))
                    h += 1;
                  else if (m & (((long) 1) << k))
                    h -= 1;
                  printf(" %2d",h);
                }
              printf(" [%d]\n",v);
            }
#endif
          }
      
        hit = 0;
        val = Val-crow;
        t   = crow+SatMax;
        if (t > rlen)
          t = rlen;
        for (; i <= t; i++)
          { e = tran[(int) (aseq[i])];
            c = 0;
            ADVANCE_BLOCK(p,m,e,c)
            v += c;
            val[i] = v;
            if (v <= SatThr)
              hit = 1;

#ifdef DEBUG_DP
            { int k, h;

              printf("  %c:",Letter[(int) aseq[i]]);
              for (k = 0; k < W; k++)
                printf("  %c",Letter[(int) aseq[r+k]]);
              printf("\n    ");
              h = 0;
              for (k = 0; k < W; k++)
                { if (p & (((long) 1) << k))
                    h += 1;
                  else if (m & (((long) 1) << k))
                    h -= 1;
                  printf(" %2d",h);
                }
              printf(" [%d]\n",v);
            }
#endif
          }
      }

      if (hit)
        { int   Sort[SatMax+2];
          Level Tree[SatMax+2];
          int   root, ntop, vlen;

          if (ctop + SatMax > Cmax)
            { Cmax = 1.2*(ctop+SatMax) + 5000;
              Cand = (Seed *) Realloc(Cand,sizeof(Seed)*Cmax,"Allocating Candidate Vector");
            }

          vlen = rlen - crow;
          if (vlen >= SatMax)
            vlen = SatMax;
          else
            Val[vlen+1] = W;

          sort_vector(Val,Sort,vlen);
          root = level_tree(Val,Sort,Tree,vlen);

	  ntop = append_micro(root,Val,Tree,Cand+ctop,crow,vlen) - Cand;
          ntop = merge_micro(Val,Cand+ctop,Cand+ntop,crow,vlen) - Cand;
          ntop = append_minis(root,Val,Tree,Cand+ntop,crow,vlen) - Cand;
#ifdef DEBUG_SEEDS
          printf("  %5d:\n",crow);
          print_vector(Val,vlen);
          if (ntop > ctop)
            print_seeds(Cand,ctop,ntop);
          // print_tree(root,Tree,0);
#endif
          ctop = ntop;
        }
    }

  { int   Sort[ctop];
    Path *apath, *bpath;
    int   c;

#ifdef DEBUG_EXPLORE

    for (c = 0; c < ctop; c++)
      printf("  (%4d,%3d-%3d) %d\n",Cand[c].row,Cand[c].lcol,Cand[c].hcol,Cand[c].val);

#endif

    sort_seeds(Cand,ctop,Sort);

    mtop  = 0;
    apath = align->path;
    for (c = 0; c < ctop; c++)
      { int low, hgh, row, len;
        int lbord, hbord;
        int per;

        if (Cand[Sort[c]].val >= W)
          continue;

        low = Cand[Sort[c]].lcol;
        hgh = Cand[Sort[c]].hcol;
        row = Cand[Sort[c]].row;

        if (low <= MICROMAX)
          lbord = Sat_Width[low];
        else
          lbord = SAT_BORD*low;
        if (hgh <= MICROMAX)
          hbord = Sat_Width[hgh];
        else
          hbord = SAT_BORD*hgh;

        bpath = Local_Alignment(align,work,spec,low,hgh,2*row+(low+hgh)/2,lbord,hbord);

        len = apath->aepos - apath->abpos;
        if (len >= low && row-apath->bbpos >= W && apath->diffs <= MAX_ERATE * len)
          { Compute_Trace_PTS(align,work,W,GREEDIEST);
            per = sat_period(align);

#ifdef DEBUG_EXPLORE

            printf("    %5d x %3d-%3d: ",row,low,hgh);
            printf(" (%5d,%5d)->(%5d,%5d) %5d/%5d = %.2f%%\n",
                    apath->abpos,apath->bbpos,apath->aepos,apath->bepos,
                    apath->diffs,len,(apath->diffs*100.)/len);

            // printf(" %5d - %5d @ %4d : %5d/%5d = %.2f%%\n",
                    // apath->bbpos,apath->aepos,per,
                    // apath->diffs,len,(apath->diffs*100.)/len);

            // Print_Reference(stdout,align,work,4,100,10,0,5);

#endif

            if (hgh > MICROMAX)
              { int ab, bb;
                int d, i, p;
                int w = .1*per;

                ab = (bpath->abpos / W) * W;
                bb = bpath->bbpos;
                d  = 0;
                while (Cand[d].row <= ab)
                  d += 1;
                for (i = 0; i < bpath->tlen; i += 2)
                  { bb += ((uint16 *) (bpath->trace))[i+1];
                    ab += W;
                    if (ab > bpath->aepos)
                      { ab = bpath->aepos;
                        break;
                      }
                    p = bb-ab;
                    while (d < ctop && Cand[d].row == ab)
                      { if (Cand[d].hcol > MICROMAX)
                          { if (p-w <= Cand[d].hcol && Cand[d].lcol <= p+w)
                              Cand[d].val = W;
                          }
                        d += 1;
                      }
                  }
              }
            else
              { int ab, bb;
                int d, i, p;

                ab = (bpath->abpos / W) * W;
                bb = bpath->bbpos;
                d  = 0;
                while (Cand[d].row <= ab)
                  d += 1;
                for (i = 0; i < bpath->tlen; i += 2)
                  { bb += ((uint16 *) (bpath->trace))[i+1];
                    ab += W;
                    if (ab > bpath->aepos)
                      { ab = bpath->aepos;
                        break;
                      }
                    p = bb-ab;
                    while (d < ctop && Cand[d].row == ab)
                      { if (hgh == p &&  Cand[d].lcol == p)
                          Cand[d].val = W;
                        d += 1;
                      }
                  }
              }

            if (mtop >= Mmax)
              { Mmax  = 1.2*mtop + 1000;
                Match = (Sat_Hit *) Realloc(Match,sizeof(Sat_Hit)*Mmax,"Allocating Match List");
              }

            Match[mtop].beg = apath->bbpos;
            Match[mtop].end = apath->aepos;
            mtop += 1;

#ifdef DEBUG_EXPLORE
            { int m;

              for (m = 0; m < ctop; m++)
                if (Cand[m].val < W)
                  printf("  (%4d,%3d-%3d) %d\n",Cand[m].row,Cand[m].lcol,Cand[m].hcol,Cand[m].val);
                else
                  printf("  (%4d,%3d-%3d) XXX\n",Cand[m].row,Cand[m].lcol,Cand[m].hcol);
            }
#endif
          }
      }
  }

  if (mtop == 0)
    return (NULL);

  { int i, j;
    int max;
    int ntop;

    qsort(Match,mtop,sizeof(Sat_Hit),SATSORT);

    Match[mtop].beg = rlen+1;
    max  = Match[0].end;
    j    = 0;
    ntop = 0;
    for (i = 1; i <= mtop; i++)
      { if (Match[i].beg >= max)
          { // printf("   Range %5d - %5d\n",Match[j].beg,max);
            Match[ntop].beg = Match[j].beg;
            Match[ntop].end = max;
            ntop += 1;
            j = i;
          }
        if (Match[i].end > max)
          max = Match[i].end;
      }
    Match[ntop].beg = -1;

    // printf("\n");
    // for (i = 0; Match[i].beg >= 0; i++)
      // printf("  %5d - %5d\n",Match[i].beg,Match[i].end);
  }

  return (Match);
}


/****************************************************************************************\
*                                                                                        *
*  SATELLITE REPEAT MODELER                                                              *
*                                                                                        *
\****************************************************************************************/

/*

static void sat_aligner(Alignment *align)
{ int   Val[SatMax+W+1], Hal[SatMax+W+1];
  WORD *_tran[SIGMA], tran[SIGMA];

  char *aseq = align->aseq;
  int   abeg = align->path->abpos;
  int   aend = align->path->aepos;

  int   lrow, crow;
  int   hrow, hlen, hlast;
  int  *val, *nal;

  hrow = 1;
  if (aend - abeg < SatMax)
    hlen = aend-abeg;
  else
    hlen = SatMax;
  hlast = 0;

  { int i;

    for (i = 1; i <= SatMax; i++)
      Val[i] = 0;
    for (i = SatMax+1; i <= SatMax+W; i++)
      Val[i] = 1;
  }

  val  = Val - abeg;
  for (lrow = abeg; lrow < aend; lrow += W)
    { int    i, t, u, c;
      WORD   e, d, p, m;

      crow = lrow+W;
      nal  = val;

      setup_search(aseq+lrow,W,_tran,tran);

      p = ALL;
      m = 0;
      d = 1;
      for (i = lrow+1; i <= crow; i++)
        { e = tran[(int) (aseq[i])] & d;
          d = (d << 1) | 1;
          c = nal[i];
          ADVANCE_BLOCK(p,m,e,c)

#ifdef DEBUG_DP
          { int k, h;

            printf("  %c:",Letter[(int) aseq[i]]);
            for (k = 0; k < W; k++)
              { if (e & (0x1ll << k))
                  printf(" +");
                else
                  printf("  ");
                printf("%c",Letter[(int) aseq[r+k]]);
              }
            printf("\n    ");
            h = 0;
            for (k = 0; k < W; k++)
              { if (p & (((long) 1) << k))
                  h += 1;
                else if (m & (((long) 1) << k))
                  h -= 1;
                printf(" %2d",h);
              }
            printf(" [%d]\n",v);
          }
#endif
        }
      
      val = Val-crow;
      t = u = crow+SatMax;
      if (t > aend)
        u = aend;
      for (; i <= u; i++)
        { e = tran[(int) (aseq[i])];
          c = nal[i];
          ADVANCE_BLOCK(p,m,e,c)
          val[i] = c;

#ifdef DEBUG_DP
          { int k, h;

            printf("  %c:",Letter[(int) aseq[i]]);
            for (k = 0; k < W; k++)
              printf("  %c",Letter[(int) aseq[r+k]]);
            printf("\n    ");
            h = 0;
            for (k = 0; k < W; k++)
              { if (p & (((long) 1) << k))
                  h += 1;
                else if (m & (((long) 1) << k))
                  h -= 1;
                printf(" %2d",h);
              }
            printf(" [%d]\n",v);
          }
#endif
        }

      if (t > u)
        { int del = t-u;
          if (del < W)
            del = W-del;
          else
            del = 0;
          for (i = del; i < W; i++)
            { if (p & (((long) 1) << i))
                hlast += 1;
              else if (m & (((long) 1) << i))
                hlast -= 1;
              Hal[hrow++] = hlast;
            }
        }
    }

  { int   Sort[SatMax+2];
    Level Tree[SatMax+2];
    int   root;

    Hal[0] = Hal[1]+1;
    Hal[hlen+1] = Hal[hlen]+1;
    // sort_vector(Hal,Sort,hlen);
    // root = level_tree(Hal,Sort,Tree,hlen);

    print_vector(Hal,hlen);
    // print_tree(root,Tree,0);
#ifdef DEBUG_SEEDS
#endif
  }
}

*/

Sat_Hit *Sat_Finder(HITS_DB *db, int i)
{ static Alignment  _align, *align = &_align;
  static Path       _path,  *path  = &_path;
  static Work_Data  *work;
  static Align_Spec *spec;
  static char       *read;
  static int         firstime = 1;

  int      rlen;
  Sat_Hit *hits;

  if (firstime)
    { firstime = 0;

      spec = New_Align_Spec( .70, W, db->freq);
      work = New_Work_Data();
      align->path = path;

      db->maxlen += W;
      read = New_Read_Buffer(db);
      db->maxlen -= W;
    }

  rlen = db->reads[i].rlen;           //  Fetch read

  // printf("\nAnalyzing %d(%x) %d\n",i+1,(db->reads[i].flags & DB_BEST) != 0,rlen);
  // fflush(stdout);

  if (rlen < 2*W)
    return (NULL);

  if (Load_Read(db,i,read,0))
    exit (1);

  align->aseq = read;
  align->bseq = read;
  align->alen = rlen;
  align->blen = rlen;

  hits = sat_filter(align,spec,work);

/*
  { Sat_Hit *h;

    for (h = hits; h->beg >= 0; h++)
      { align->path->abpos = h->beg;
        align->path->aepos = h->end;
        sat_aligner(align);
      }
  }
*/

  return (hits);
}
