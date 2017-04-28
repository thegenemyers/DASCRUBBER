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
 *  Map and extend every overlap to the patched read framework.
 *
 *  Author:  Gene Myers
 *  Date  :  March 2015
 *
 *******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>
#include <stdint.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "DB.h"
#include "align.h"

#undef    OUTLINE
#undef    SHOW_MAP
#undef    TRACE
#undef    SHOW_FINAL
#undef    SHOW_ALIGNMENTS

static char *Usage = " [-v] [-l<int(800)>] <block1:db> <block2:db> <source:las> <target:las>";

static int     TRACE_SPACING;  //  Trace spacing (from .las file)
static int     TBYTES;         //  Bytes per trace segment (from .las file)
static int     SMALL;          //  Trace points can fit in a byte
static int     MIN_LEN;        //  Minimum piece length

static HITS_DB _ADB, *ADB = &_ADB;    //  A-read database
static HITS_DB _BDB, *BDB = &_BDB;    //  B-read database

static int ADB_ofirst, ADB_olast;
static int BDB_ofirst, BDB_olast;

static int AFIRST, BFIRST;

static int64 *AMAP_IDX;         //  Map to originals for A-reads
static int   *AMAP;
static int64 *BMAP_IDX;         //  Map to orignals for B-reads
static int   *BMAP;

static int   *IAMAP, *IBMAP;    //  Inverse map, old x -> new IMAP[x]..IMAP[x+1]-1

static FILE   *OUTPUT;          //  The new set of overlaps
static int64   WOVLS;

static int     VERBOSE;

int lowTK(int a)
{ return (a/TRACE_SPACING); }

int hghTK(int a)
{ return ((a+(TRACE_SPACING-1))/TRACE_SPACING); }

int lowTP(int a)
{ return ((a/TRACE_SPACING)*TRACE_SPACING); }


/*******************************************************************************************
 *
 *  Finger iterator: allows one to map the next trace point of a read to its
 *     patched read as the trace points are examined in order.
 *
 *******************************************************************************************/

typedef struct
  { int cidx;
    int lidx;
    int dist;
    int last;
    int blen;
    int *map;
  } Finger;

  //  A finger is initialized with init_finger where cur is suppled by the user, and the
  //     patch sequence is in GRIM[gb..ge].

static inline void init_finger(Finger *f, int *map, int mb, int me, int blen)
{ if (blen == 0)
    { f->cidx = mb+3;
      f->lidx = me;
      f->last = map[mb];
    }
  else
    { f->cidx = me-4;
      f->lidx = mb-1;
      f->last = blen - map[me-1];
    }
  f->dist = 0;
  f->blen = blen;
  f->map  = map;
}

  //  Advance finger to position pos and return position in patched read, if known, -1 otherwise

static inline int good(Finger *cur, int pos)
{ int blen, *map;

  map  = cur->map;
  blen = cur->blen;
  if (blen == 0)
    { while (cur->cidx < cur->lidx && pos >= map[cur->cidx])
        { cur->dist += (map[cur->cidx-2] - cur->last) + map[cur->cidx-1];
          cur->last  = map[cur->cidx];
          cur->cidx += 3; 
        }
      if (pos <= map[cur->cidx-2])
        { if (pos < cur->last)
            return (-1);
          else
            return (cur->dist + (pos-cur->last));
        }
    }
  else
    { while (cur->cidx > cur->lidx && pos >= blen - map[cur->cidx])
        { cur->dist += ((blen - map[cur->cidx+2]) - cur->last) + map[cur->cidx+1];
          cur->last  = blen - map[cur->cidx];
          cur->cidx -= 3; 
        }
      if (pos <= blen - map[cur->cidx+2])
        { if (pos < cur->last)
            return (-1);
          else
            return (cur->dist + (pos-cur->last));
        }
    }
  return (-1);
}

  //  Advance finger to position pos, and return best estimate of position in patched read,
  //    or -1 if outside the bounds of the patched read.  acc points at the distance the
  //    estimate is from a non-patched segment (0 if mapped).

static inline int where(Finger *cur, int pos, int *acc)
{ int blen, *map;

  map  = cur->map;
  blen = cur->blen;
  if (blen == 0)
    { while (cur->cidx < cur->lidx && pos >= map[cur->cidx])
        { cur->dist += (map[cur->cidx-2] - cur->last) + map[cur->cidx-1];
          cur->last  = map[cur->cidx];
          cur->cidx += 3; 
        }
      if (pos <= map[cur->cidx-2])
        { if (pos < cur->last)
            return (-1);
          else
            { *acc = 0;
              return (cur->dist + (pos-cur->last));
            }
        }
      if (cur->cidx >= cur->lidx)
        return (-1);
      else
        { int ab, ae;
    
          ab = map[cur->cidx-2];
          ae = map[cur->cidx];
          if (pos-ab < ae-pos)
            *acc = pos-ab;
          else
            *acc = ae-pos;
          return (cur->dist + (ab-cur->last) + ((1.*(pos-ab))/(ae-ab)) * map[cur->cidx-1]);
        }
    }
  else
    { while (cur->cidx > cur->lidx && pos >= blen - map[cur->cidx])
        { cur->dist += ((blen - map[cur->cidx+2]) - cur->last) + map[cur->cidx+1];
          cur->last  = blen - map[cur->cidx];
          cur->cidx -= 3; 
        }
      if (pos <= blen - map[cur->cidx+2])
        { if (pos < cur->last)
            return (-1);
          else
            { *acc = 0;
              return (cur->dist + (pos-cur->last));
            }
        }
      if (cur->cidx <= cur->lidx)
        return (-1);
      else
        { int ab, ae;
    
          ab = blen - map[cur->cidx+2];
          ae = blen - map[cur->cidx];
          if (pos-ab < ae-pos)
            *acc = pos-ab;
          else
            *acc = ae-pos;
          return (cur->dist + (ab-cur->last) + ((1.*(pos-ab))/(ae-ab)) * map[cur->cidx+1]);
        }
    }
}


/*******************************************************************************************
 *
 *  Trace point mapping:
 *     recon makes a trace with all the mapable pairs
 *     when there are no mapable pairs, estimate finds the estimated point position closest
 *       to a mapable region.
 *
 *******************************************************************************************/

static int recon(Path *image, Path *path, Finger *afinger, Finger *bfinger)
{ static int     tmax = -1;
  static uint16 *itrace = NULL;

  int ae, be;
  int al, bl;
  int an, bn;
  int t, tl, df;
  uint16 *strace = ((uint16 *) (path->trace));

  if (path->tlen > tmax)
    { tmax = 1.2*path->tlen + 100;
      itrace = (uint16 *) Realloc(itrace,sizeof(uint16)*tmax,"Reallocating image trace");
    }
  image->trace = itrace;

#ifdef TRACE
  printf("      Backbone:\n");
  fflush(stdout);
#endif
  df = 0;
  tl = -1;
  al = bl = 0;
  ae = lowTP(path->abpos);
  be = path->bbpos;
  if ((an = good(afinger,path->abpos)) >= 0 && (bn = good(bfinger,be)) >= 0)
    { image->abpos = al = an;
      image->bbpos = bl = bn;
      tl = 0;
#ifdef TRACE
      printf("        %5d,%5d -> %5d,%5d\n",path->abpos,be,an,bn);
      fflush(stdout);
#endif
    }
  for (t = 1; t < path->tlen; t += 2)
    { ae += TRACE_SPACING;
      be += strace[t];
      if (ae > path->aepos)
        ae = path->aepos;
      if (tl >= 0)
        df += strace[t-1];
      if ((an = good(afinger,ae)) >= 0 && (bn = good(bfinger,be)) >= 0)
        { if (tl < 0)
            { image->abpos = an;
              image->bbpos = bn;
              tl = 0;
            }
          else
            { itrace[tl]   = an-al;
              itrace[tl+1] = bn-bl;
              tl += 2;
            }
          image->aepos = al = an;
          image->bepos = bl = bn;
          image->diffs = df;
#ifdef TRACE
          printf("        %5d,%5d -> %5d,%5d\n",ae,be,an,bn);
          fflush(stdout);
#endif
        }
    }

  image->tlen = tl;
  if (tl <= 0)
    return (0);
  else
    return (1);
}

static int estimate(Path *path, Finger *afinger, Finger *bfinger, int *bsta, int *bstb, int *acc)
{ int ae, be;
  int an, bn;
  int best, adst, bdst;
  int t;
  uint16 *strace = ((uint16 *) (path->trace));

  *bsta = *bstb = -1;
  best  = INT32_MAX;

#ifdef TRACE
  printf("      Point Estimate:\n");
  fflush(stdout);
#endif
  ae = lowTP(path->abpos);
  be = path->bbpos;
  if ((an = where(afinger,path->abpos,&adst)) >= 0 && (bn = where(bfinger,be,&bdst)) >= 0)
    { best = adst + bdst;
      *bsta = an;
      *bstb = bn;
#ifdef TRACE
      printf("        %5d,%5d -> %5d(%d),%5d(%d)\n",path->abpos,be,an,adst,bn,bdst);
      fflush(stdout);
#endif
    }
  for (t = 1; t < path->tlen; t += 2)
    { ae += TRACE_SPACING;
      be += strace[t];
      if (ae > path->aepos)
        ae = path->aepos;
      if ((an = where(afinger,ae,&adst)) >= 0 && (bn = where(bfinger,be,&bdst)) >= 0)
        { if (adst + bdst < best)
            { best = adst + bdst;
              *bsta = an;
              *bstb = bn;
            }
#ifdef TRACE
          printf("        %5d,%5d -> %5d(%d),%5d(%d)\n",ae,be,an,adst,bn,bdst);
          fflush(stdout);
#endif
        }
    }
  *acc = best;

  return (*bsta >= 0);
}

#ifdef SHOW_MAP

static void print_map(int *map, int mb, int me, int clen)
{ int b, dist;

  if (clen == 0)
    { printf(" n");
      for (b = mb; b < me; b += 3)
        { printf(" [%5d,%5d]",map[b],map[b+1]);
          if (b+2 < me)
            printf(" %5d",map[b+2]);
        }
      printf("\n");

      printf("  ");
      dist = 0;
      for (b = mb; b < me; b += 3)
        { printf(" [%5d,%5d]",dist,dist+(map[b+1]-map[b]));
          dist += map[b+1]-map[b];
          if (b+2 < me)
            { printf(" %5d",map[b+2]);
              dist += map[b+2];
            }
        }
      printf("\n");
    }
  else
    { printf(" c");
      for (b = me; b >= mb; b -= 3)
        { printf(" [%5d,%5d]",clen-map[b-1],clen-map[b-2]);
          if (b-3 > mb)
            printf(" %5d",map[b-3]);
        }
      printf("\n");

      printf("  ");
      dist = 0;
      for (b = me; b >= mb; b -= 3)
        { printf(" [%5d,%5d]",dist,dist+(map[b-1]-map[b-2]));
          dist += map[b-1]-map[b-2];
          if (b-3 > mb)
            { printf(" %5d",map[b-3]);
              dist += map[b-3];
            }
        }
      printf("\n");
    }
}

#endif

#ifdef SHOW_FINAL

static void show_overlap(Overlap *ovl)
{ int     i, a, b;
  uint16 *t;

  t = (uint16 *) (ovl->path.trace);
  a = ovl->path.abpos;
  b = ovl->path.bbpos;
  for (i = 0; i < ovl->path.tlen; i += 2)
    { a += t[i];
      b += t[i+1];
      printf("        %5d %5d :: %5d %5d\n",t[i],t[i+1],a,b);
      fflush(stdout);
    }
}

#endif

static void convert_trace(Path *path)
{ int ab, ae;
  int t;
  uint16 *trace = ((uint16 *) (path->trace));
  
  ae = lowTP(path->abpos);
  ab = path->abpos;
  for (t = 0; t < path->tlen; t += 2)
    { ae += TRACE_SPACING;
      if (ae > path->aepos)
        ae = path->aepos;
#ifdef TRACE
      printf("        %5d,%5d -> %5d,%5d\n",trace[t],trace[t+1],ae-ab,trace[t+1]);
      fflush(stdout);
#endif
      trace[t] = ae-ab;
      ab = ae;
    }
}


//  Produce the concatentation of path1 and path2 where they are known to meet at
//    the trace point with coordinate ap. Place this result in a big growing buffer,
//    that gets reset when fusion is called with path1 = NULL

static void fusion(Path *path1, Path *path2, int wch)
{ static uint16  *paths = NULL;
  static int      pmax  = 0;
  static int      ptop  = 0;

  int     k;
  int     len;
  uint16 *trace;

  if (path1 == NULL)
    { ptop = 0;
      return;
    }

  len = path1->tlen + path2->tlen;

  if (ptop + len >= pmax)
    { pmax = 1.2*(ptop+len) + 1000;
      paths = (uint16 *) Realloc(paths,sizeof(uint16)*pmax,"Allocating paths");
      if (paths == NULL)
        exit (1);
    }
  trace = paths+ptop;
  ptop += len;

  len  = 0;
  if (path1->tlen > 0)
    { uint16 *t = (uint16 *) (path1->trace);
      for (k = 0; k < path1->tlen; k += 2)
        { trace[len++] = t[k];
          trace[len++] = t[k+1];
        }
    }
  if (path2->tlen > 0)
    { uint16 *t = (uint16 *) (path2->trace);
      for (k = 0; k < path2->tlen; k += 2)
        { trace[len++] = t[k];
          trace[len++] = t[k+1];
        }
    }

  if (wch == 1)
    { path1->aepos  = path2->aepos;
      path1->bepos  = path2->bepos;
      path1->diffs += path2->diffs;
      path1->trace  = trace;
      path1->tlen   = len;
    }
  else
    { path2->abpos  = path1->abpos;
      path2->bbpos  = path1->bbpos;
      path2->diffs += path1->diffs;
      path2->trace  = trace;
      path2->tlen   = len;
    }
}

static void EXTENDER(int aread, Overlap *ovls, int novl)
{ Finger      _afinger, *afinger = &_afinger;
  Finger      _bfinger, *bfinger = &_bfinger;

  static Overlap     _ovla, *ovla = &_ovla;
  static Path        *ipath = &_ovla.path;

  static Path         rpath, fpath;
  static Alignment   _ralign, *ralign = &_ralign;
  static Alignment   _falign, *falign = &_falign;

  static Work_Data  *work = NULL;
  static Align_Spec *spec;

  int ap, alast;

  if (aread < ADB_ofirst || aread >= ADB_olast)
    return;

  if (work == NULL)
    { spec = New_Align_Spec(.70,100,ADB->freq);
      work = New_Work_Data();
      ralign->path = &rpath;
      falign->path = &fpath;
    }

  alast = IAMAP[aread+1];
  for (ap = IAMAP[aread]; ap < alast; ap++)
    { int mb, me;
      int aend, abeg, alen;

      mb = AMAP_IDX[ap]+2;
      me = AMAP_IDX[ap+1];

      abeg = AMAP[mb];
      aend = AMAP[me-1];
      if (aend - abeg < MIN_LEN)
        continue;

      ralign->aseq = falign->aseq = ((char *) ADB->bases) + ADB->reads[ap].boff;
      ralign->alen = falign->alen = alen = ADB->reads[ap].rlen;

#ifdef OUTLINE
      printf("AREAD %d -> %d [%d,%d]\n",aread,ap,abeg,aend);
      fflush(stdout);
#endif
#ifdef SHOW_MAP
      print_map(AMAP,mb,me,0);
#endif

      { int   o, ob, oe;
        Path *path;
        int   bread;
        int   bp, blast;

        for (ob = 0; ob < novl; ob = oe)
          { bread = ovls[ob].bread;
            for (oe = ob+1; oe < novl && ovls[oe].bread == bread; oe += 1)
              ;
            if (bread < BDB_ofirst || bread >= BDB_olast)
              continue;

            blast = IBMAP[bread+1];
            for (bp = IBMAP[bread]; bp < blast; bp++)
              { int hb, he;
                int bend, bbeg, blen;
                int alpos, clen;

                hb = BMAP_IDX[bp]+2;
                he = BMAP_IDX[bp+1];

                bbeg = BMAP[hb];
                bend = BMAP[he-1];
                if (bend - bbeg < MIN_LEN)
                  continue;

#ifdef OUTLINE
		printf("  BREAD %d->%d [%d,%d]\n",bread,bp,bbeg,bend);
                fflush(stdout);
#endif
#ifdef SHOW_MAP
		print_map(BMAP,hb,he,0);
#endif

                alpos = -1;
                for (o = ob; o < oe; o++)
                  { path = &(ovls[o].path);
#ifdef OUTLINE
                    printf("    OVL %d: [%d,%d] [%d,%d]\n",o,
                           path->abpos,path->aepos,path->bbpos,path->bepos);
                    fflush(stdout);
#endif

                    if (path->abpos <= aend-MIN_LEN && path->aepos >= abeg+MIN_LEN &&
                        path->bbpos <= bend-MIN_LEN && path->bepos >= bbeg+MIN_LEN)

                      { ralign->bseq  = falign->bseq  = ((char *) BDB->bases) + BDB->reads[bp].boff;
                        ralign->blen  = falign->blen  = blen = BDB->reads[bp].rlen;
                        ralign->flags = falign->flags = ovls[o].flags;

                        if (COMP(ralign->flags))
                          { clen = BMAP[hb-1];
                            Complement_Seq(ralign->bseq,blen);
                          }
                        else
                          clen = 0;

#ifdef SHOW_MAP
                        print_map(BMAP,hb,he,clen);
#endif

                        init_finger(afinger,AMAP,mb,me,0);
                        init_finger(bfinger,BMAP,hb,he,clen);
                        if ( ! recon(ipath,path,afinger,bfinger))
                          { int apos, bpos, acc, len, diag;

                            init_finger(afinger,AMAP,mb,me,0);
                            init_finger(bfinger,BMAP,hb,he,clen);
                            if (estimate(path,afinger,bfinger,&apos,&bpos,&acc))
                              if (apos > alpos)
                                { diag = apos-bpos;

                                  acc /= 2;
                                  if (apos + acc > alen)
                                    acc = alen-apos;
                                  if (bpos + acc > blen)
                                    acc = blen-bpos;
                                  if (apos < acc)
                                    acc = apos;
                                  if (bpos < acc)
                                    acc = bpos;
                                  if (acc > 500)
                                    acc = 500;
                                  acc *= 2;

#ifdef OUTLINE
                                  printf("      Trying: %d,%d + %d\n",apos,bpos,acc);
                                  fflush(stdout);
#endif
                                  Local_Alignment(ralign,work,spec,
                                                  diag-acc,diag+acc,apos+bpos,-1,-1);
#ifdef OUTLINE
                                  printf("      Local: <%d,%d> -> <%d,%d>\n",
                                         ralign->path->abpos,ralign->path->bbpos,
                                         ralign->path->aepos,ralign->path->bepos);
                                  fflush(stdout);
#endif
                                  len = ralign->path->aepos - ralign->path->abpos;
                                  if (len >= MIN_LEN && ralign->path->diffs <= .35*len)
                                    { ovla->aread = ap;
                                      ovla->bread = bp;
                                      ovla->flags = ralign->flags;
                                      _ovla.path  = *(ralign->path);
                                      convert_trace(ralign->path);
#ifdef OUTLINE
                                      printf("      Final: %d[%d..%d] vs %d[%d..%d] %c d=%d\n",
                                             ovla->aread,ovla->path.abpos,ovla->path.aepos,
                                             ovla->bread,ovla->path.bbpos,ovla->path.bepos,
                                             (COMP(ovla->flags) ? 'c' : 'n'), ovla->path.diffs);
                                      fflush(stdout);
#endif
#ifdef SHOW_FINAL
                                      show_overlap(ovla);
#endif
#ifdef SHOW_ALIGNMENTS
                                      Compute_Trace_IRR(ralign,work,GREEDIEST);
                                      Print_Alignment(stdout,ralign,work,4,80,10,0,6);
#else
                                      Write_Overlap(OUTPUT,ovla,sizeof(uint16));
                                      WOVLS += 1;
#endif

                                      alpos  = ralign->path->aepos;
#ifdef OUTLINE
                                      printf("        ACCEPT\n");
#endif
                                    }
#ifdef OUTLINE
                                  else
                                    printf("        REJECT\n");
#endif
                                }
#ifdef OUTLINE
                              else
                                printf("        SKIP\n");
                            else
                              printf("        NO OVERLAP\n");
                            fflush(stdout);
#endif
                          }

                        else if (ipath->aepos > alpos)
                          { int ab, bb;
                            int ae, be;
                            int ar, br;
                            int af, bf;

                            ab = ipath->abpos;
                            bb = ipath->bbpos;
                            ae = ipath->aepos;
                            be = ipath->bepos;

                            ar = ab;
                            br = bb;
                            if (ab > 0 && bb > 0)
                              { Find_Extension(ralign,work,spec,ab-bb,ab+bb,-1,-1,1);
                                ar = ralign->path->abpos;
                                br = ralign->path->bbpos;
#ifdef OUTLINE
                                printf("      Rev: (%d,%d)",ab,bb);
                                printf(" -> (%d,%d)",ar,br);
                                printf(" %d",ralign->path->diffs);
                                fflush(stdout);
#endif
                                if (ar == 0 || br == 0 || ralign->path->diffs <= .35*(ab-ar))
                                  {
#ifdef OUTLINE
                                    printf(" OK\n");
                                    fflush(stdout);
#endif
                                    if (ab - 10 < ar)
                                      { uint16 *trace = (uint16 *) ipath->trace;
                                        int     tlen  = ipath->tlen;

                                        if (tlen > 0)
                                          { trace[0] += ab - ar;
                                            trace[1] += bb - br;
                                          }
                                        ipath->abpos = ar;
                                        ipath->bbpos = br;
                                      }
                                    else
                                      { convert_trace(ralign->path);
                                        fusion(ralign->path,ipath,2);
                                      }
                                  }
                                else
                                  { ar = ab;
                                    br = bb;
#ifdef OUTLINE
                                    printf(" NOTOK\n");
                                    fflush(stdout);
#endif
                                  }
                              }

                            af = ae;
                            bf = be;
                            if (ae < alen && be < blen)
                              { Find_Extension(falign,work,spec,ae-be,ae+be,-1,-1,0);
                                af = falign->path->aepos;
                                bf = falign->path->bepos;
#ifdef OUTLINE
                                printf("      Fow: (%d,%d)",ae,be);
                                printf(" -> (%d,%d)",af,bf);
                                printf(" %d",falign->path->diffs);
                                fflush(stdout);
#endif
                                if (af == alen || bf == blen || falign->path->diffs <= .35*(af-ae))
                                  {
#ifdef OUTLINE
                                    printf(" OK\n");
                                    fflush(stdout);
#endif
                                    if (ae + 10 > af)
                                      { uint16 *trace = (uint16 *) ipath->trace;
                                        int     tlen  = ipath->tlen;

                                        if (tlen > 0)
                                          { trace[tlen-2] += af-ae;
                                            trace[tlen-1] += bf-be;
                                          }
                                        ipath->aepos = af;
                                        ipath->bepos = bf;
                                      }
                                    else
                                      { convert_trace(falign->path);
                                        fusion(ipath,falign->path,1);
                                      }
                                  }
                                else
                                  { af = ae;
                                    bf = be;
#ifdef OUTLINE
                                    printf(" NOTOK\n");
                                    fflush(stdout);
#endif
                                  }
                              }

                            alpos = af;
                            if (af-ar >= MIN_LEN) 
                              { ovla->aread = AFIRST + ap;
                                ovla->bread = BFIRST + bp;
                                ovla->flags = ralign->flags;
#ifdef OUTLINE
                                printf("      Final: %d[%d..%d] vs %d[%d..%d] %c d=%d\n",
                                       ovla->aread,ovla->path.abpos,ovla->path.aepos,
                                       ovla->bread,ovla->path.bbpos,ovla->path.bepos,
                                       (COMP(ovla->flags) ? 'c' : 'n'), ovla->path.diffs);
                                fflush(stdout);
#endif
#ifdef SHOW_FINAL
                                show_overlap(ovla);
#endif
#ifdef SHOW_ALIGNMENTS
				fpath = *ipath;
                                Compute_Trace_IRR(falign,work,GREEDIEST);
                                Print_Alignment(stdout,falign,work,4,80,10,0,6);
#else
                                Write_Overlap(OUTPUT,ovla,sizeof(uint16));
                                WOVLS += 1;
#endif
                              }
#ifdef OUTLINE
                            else
                              printf("        NO OVERLAP\n");
#endif
                            fusion(NULL,NULL,0);
                          }

#ifdef OUTLINE
                        else
                          printf("     SKIP\n");
                        fflush(stdout);
#endif

                        if (COMP(ralign->flags))
                          Complement_Seq(ralign->bseq,blen);
                      }
                  }
              }
          }
      }
    }
}


static int make_a_pass(FILE *input, void (*ACTION)(int, Overlap *, int), int trace)
{ static Overlap *ovls  = NULL;
  static int      omax  = 500;
  static uint16  *paths = NULL;
  static int      pmax  = 100000;

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
    { TBYTES = sizeof(uint8);
      SMALL  = 1;
    }
  else
    { TBYTES = sizeof(uint16);
      SMALL  = 0;
    }

  Read_Overlap(input,ovls);
  if (trace)
    { if (ovls[0].path.tlen > pmax)
        { pmax = 1.2*(ovls[0].path.tlen)+10000;
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

  pcur = 0;
  n = max = 0;
  for (j = ovls[0].aread; j < ADB_olast; j++)
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

      if (j >= ADB_ofirst)
        ACTION(j,ovls,n);
    }

  return (max);
}


int main(int argc, char *argv[])
{ HITS_TRACK *map1, *map2;

  //  Process arguments

  { int  i, j, k;
    int  flags[128];
    char *eptr;

    ARG_INIT("DASrealign")

    MIN_LEN = 800;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("v")
            break;
          case 'l':
            ARG_POSITIVE(MIN_LEN,"Minimum piece length to recompute")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];

    if (argc != 5)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }
  }

  { int status;

    status = Open_DB(argv[1],ADB);
    if (status < 0)
      exit (1);
    if (status == 1)
      { fprintf(stderr,"%s: Cannot be called on a .dam index: %s\n",Prog_Name,argv[1]);
        exit (1);
      }
    if ( ! ADB->part)
      { fprintf(stderr,"%s: Must be called on a block: %s\n",Prog_Name,argv[1]);
        exit (1);
      }
    Trim_DB(ADB);

    Read_All_Sequences(ADB,0);
  }

  if (strcmp(argv[1],argv[2]) == 0)
    BDB = ADB;
  else
    { int status;

      status = Open_DB(argv[2],BDB);
      if (status < 0)
        exit (1);
      if (status == 1)
        { fprintf(stderr,"%s: Cannot be called on a .dam index: %s\n",Prog_Name,argv[1]);
          exit (1);
        }
      if ( ! BDB->part)
        { fprintf(stderr,"%s: Must be called on a block: %s\n",Prog_Name,argv[1]);
          exit (1);
        }
      Trim_DB(BDB);


      Read_All_Sequences(BDB,0);
    }
  AFIRST = ADB->tfirst;
  BFIRST = BDB->tfirst;

  map1 = Load_Track(ADB,"map");
  if (map1 != NULL)
    { int i, o, q, n; 
 
      AMAP_IDX = (int64 *) map1->anno;
      AMAP     = (int *) map1->data;
      for (i = 0; i <= ADB->nreads; i++)
        AMAP_IDX[i] /= sizeof(int);

      ADB_ofirst = AMAP[AMAP_IDX[0]];
      ADB_olast  = AMAP[AMAP_IDX[ADB->nreads-1]]+1;
      IAMAP = (int *) Malloc(sizeof(int)*((ADB_olast-ADB_ofirst)+1),"Inverse map") - ADB_ofirst;

      IAMAP[q = ADB_olast] = n = ADB->nreads;
      for (i = ADB->nreads-1; i >= 0; i--)
        { o = AMAP[AMAP_IDX[i]];
          if (q > o)
            while (--q > o)
              IAMAP[q] = n;
          IAMAP[o] = n = i;
        }
    }
  else
    { fprintf(stderr,"%s: Must have a 'map' track, run DASedit\n",Prog_Name);
      exit (1);
    }

  if (BDB == ADB)
    { map2 = map1;
      BMAP_IDX = AMAP_IDX;
      BMAP = AMAP;
      BDB_ofirst = ADB_ofirst;
      BDB_olast  = ADB_olast ;
      IBMAP = IAMAP;
    }
  else
    { map2 = Load_Track(BDB,"map");
      if (map2 != NULL)
        { int i, o, q, n; 
     
          BMAP_IDX = (int64 *) map2->anno;
          BMAP     = (int *) map2->data;
          for (i = 0; i <= BDB->nreads; i++)
            BMAP_IDX[i] /= sizeof(int);

          BDB_ofirst = BMAP[BMAP_IDX[0]];
          BDB_olast  = BMAP[BMAP_IDX[BDB->nreads-1]]+1;
          IBMAP = (int *) Malloc(sizeof(int)*((BDB_olast-BDB_ofirst)+1),"Inverse map") - BDB_ofirst;

          IBMAP[q = BDB_olast] = n = BDB->nreads;
          for (i = BDB->nreads-1; i >= 0; i--)
            { o = BMAP[BMAP_IDX[i]];
              if (q > o)
                while (--q > o)
                  IBMAP[q] = n;
              IBMAP[o] = n = i;
            }
        }
      else
        { fprintf(stderr,"%s: Must have a 'map' track, run DASedit\n",Prog_Name);
          exit (1);
        }
    }

  ADB_ofirst = AMAP[AMAP_IDX[0]];
  BDB_ofirst = BMAP[BMAP_IDX[0]];
  ADB_olast  = AMAP[AMAP_IDX[ADB->nreads-1]]+1;
  BDB_olast  = BMAP[BMAP_IDX[BDB->nreads-1]]+1;

  //  Open .las and process piles therein output new piles to F.las

  { FILE *input;
    char *las, *pwd;
    char *lasT, *pwdT;
    int64 novl;

    las = Root(argv[3],".las");
    pwd = PathTo(argv[3]);

    lasT = Root(argv[4],".las");
    pwdT = PathTo(argv[4]);

    if (strcmp(las,lasT) == 0 && strcmp(pwd,pwdT) == 0)
      { fprintf(stderr,"%s: source and target are the same !\n",Prog_Name);
        exit (1);
      }

    input  = Fopen(Catenate(pwd,"/",las,".las"),"r");
    OUTPUT = Fopen(Catenate(pwdT,"/",lasT,".las"),"w");
    if (input == NULL || OUTPUT == NULL)
      exit (1);
    free(pwd);
    free(las);

    WOVLS = 0;
    TRACE_SPACING = 0;
    fwrite(&WOVLS,sizeof(int64),1,OUTPUT);
    fwrite(&TRACE_SPACING,sizeof(int),1,OUTPUT);

    fread(&novl,sizeof(int64),1,input);
    fread(&TRACE_SPACING,sizeof(int),1,input);

    make_a_pass(input,EXTENDER,1);

    rewind(OUTPUT);
    fwrite(&WOVLS,sizeof(int64),1,OUTPUT);
    fclose(OUTPUT);
  }

  exit (0);
}
