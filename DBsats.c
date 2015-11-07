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
 *  Scan every read looking for micro- and min-satellites (up to 256bp module) and record
 *    any found intervals in a "sats" track.
 *
 *  Author:  Gene Myers
 *  Date  :  March 2015
 *
 *******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "DB.h"
#include "align.h"
#include "seqpats.h"

#ifdef HIDE_FILES
#define PATHSEP "/."
#else
#define PATHSEP "/"
#endif

#undef SAT_DEBUG

static char *Usage = " [-v] <source:db|dam>";

static int     VERBOSE;

static HITS_DB _DB, *DB  = &_DB;   //  Data base

static FILE   *SAT_AFILE;   //  .sats.anno
static FILE   *SAT_DFILE;   //  .sats.data
static int64   SAT_INDEX;   //  Current index into .sats.data file

//  Statistics

static int64 nsats, nbps;

int main(int argc, char *argv[])
{ char *root, *pwd;
  int   status;

  //  Process arguments

  { int  i, j, k;
    int  flags[128];

    ARG_INIT("DBsats")

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

    if (argc != 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }

    status = Open_DB(argv[1],DB);
    if (status < 0)
      exit (1);
  }

  //   Set up QV and trimming tracks

  pwd   = PathTo(argv[1]);
  if (status)
    root = Root(argv[1],".dam");
  else
    root = Root(argv[1],".db");

  SAT_AFILE = Fopen(Catenate(pwd,PATHSEP,root,".sats.anno"),"w");
  SAT_DFILE = Fopen(Catenate(pwd,PATHSEP,root,".sats.data"),"w");
  if (SAT_AFILE == NULL || SAT_DFILE == NULL)
    exit (1);

  free(pwd);
  free(root);

  { int size;

    size = 0;
    fwrite(&(DB->nreads),sizeof(int),1,SAT_AFILE);
    fwrite(&size,sizeof(int),1,SAT_AFILE);
    SAT_INDEX = 0;
    fwrite(&SAT_INDEX,sizeof(int64),1,SAT_AFILE);
  }

  //  Initialize statistics gathering

  if (VERBOSE)
    { nsats = 0;
      nbps  = 0;
    }

  //  Process each read

  { int      i, j;
    Sat_Hit *hit;

    for (i = 0; i < DB->nreads; i++)
      {
#ifdef SAT_DEBUG
        int rlen = DB->reads[i].rlen;

        printf("Analyzing %d(%x) %d\n",i+1,(DB->reads[i].flags & DB_BEST) != 0,rlen);
#endif

        hit = Sat_Finder(DB,i);
        if (hit != NULL)
          for (j = 0; hit[j].beg >= 0; j++)
            { int v;
              v = hit[j].beg;
              fwrite(&v,sizeof(int),1,SAT_DFILE);
              v = hit[j].end;
              fwrite(&v,sizeof(int),1,SAT_DFILE);
              SAT_INDEX += 2*sizeof(int);
              nsats += 1;
              nbps  += hit[j].end - hit[j].beg;
#ifdef SAT_DEBUG
              printf(" %5d - %5d\n",hit[j].beg,hit[j].end);
#endif
            }
        fwrite(&SAT_INDEX,sizeof(int64),1,SAT_AFILE);
      }
  }

  //  If verbose output statistics summary to stdout

  if (VERBOSE)
    { printf("  Found %lld satellites, totaling ",nsats);
      Print_Number(nbps,0,stdout);
      printf(" base pairs\n");
    }

  //  Clean up

  fclose(SAT_AFILE);
  fclose(SAT_DFILE);

  free(Prog_Name);

  exit (0);
}
