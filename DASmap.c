/*******************************************************************************************
 *
 *  Display a specified set of reads of a database in fasta format.
 *
 *  Author:  Gene Myers
 *  Date  :  September 2013
 *  Mod   :  With DB overhaul, made this a routine strictly for printing a selected subset
 *             and created DB2fasta for recreating all the fasta files of a DB
 *  Date  :  April 2014
 *  Mod   :  Added options to display QV streams
 *  Date  :  July 2014
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "DB.h"

#ifdef HIDE_FILES
#define PATHSEP "/."
#else
#define PATHSEP "/"
#endif

static char Usage[] = "[-p] <path:db|dam> [ <reads:FILE> | <reads:range> ... ]";

#define LAST_READ_SYMBOL   '$'
#define MAX_BUFFER       10001

typedef struct
  { FILE  *input;
    int    lineno;
    int    read;
    int    beg;
    int    end;
  } File_Iterator;

File_Iterator *init_file_iterator(FILE *input)
{ File_Iterator *it;

  it = Malloc(sizeof(File_Iterator),"Allocating file iterator");
  it->input = input;
  it->lineno = 1;
  rewind(input);
  return (it);
}

int next_read(File_Iterator *it)
{ static char nbuffer[MAX_BUFFER];

  char *eol;
  int   x;

  if (fgets(nbuffer,MAX_BUFFER,it->input) == NULL)
    { if (feof(it->input))
        return (1);
      SYSTEM_READ_ERROR;
    }
  if ((eol = index(nbuffer,'\n')) == NULL)
    { fprintf(stderr,"%s: Line %d in read list is longer than %d chars!\n",
                     Prog_Name,it->lineno,MAX_BUFFER-1);
      return (1);
    }
  *eol = '\0';
  x = sscanf(nbuffer," %d %d %d",&(it->read),&(it->beg),&(it->end));
  if (x == 1)
    it->beg = -1;
  else if (x != 3)
    { fprintf(stderr,"%s: Line %d of read list is improperly formatted\n",Prog_Name,it->lineno);
      return (1);
    }
  it->lineno += 1;
  return (0);
}

int main(int argc, char *argv[])
{ DAZZ_DB     _db, *db = &_db;
  DAZZ_TRACK *map;

  int            reps, *pts;
  int            input_pts;
  File_Iterator *iter = NULL;
  FILE          *input;

  int PRETTY;

  //  Process arguments

  { int   i, j, k;
    int   flags[128];

    ARG_INIT("DASmap")

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        ARG_FLAGS("p")
      else
        argv[j++] = argv[i];
    argc = j;

    PRETTY = flags['p'];

    if (argc <= 1)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -p: Pretty print (vs easy to parse).\n");
        exit (1);
      }
  }

  //  Open DB or DAM, and if a DAM open also .hdr file

  { int   status, kind;

    status = Open_DB(argv[1],db);
    if (status < 0)
      exit (1);
    if (status == 1)
      { fprintf(stderr,"%s: Cannot be called on a .dam index: %s\n",Prog_Name,argv[1]);
        exit (1);
      }
    if (db->part)
      { fprintf(stderr,"%s: Cannot be called on a block: %s\n",Prog_Name,argv[1]);
        exit (1);
      }

    status = Check_Track(db,"map",&kind);
    if (status == -2)
      { fprintf(stderr,"%s: 'map' track not found.\n",Prog_Name);
        exit (1);
      }
    else if (status == -1)
      { fprintf(stderr,"%s: Warning: 'map' track not sync'd with db.\n",Prog_Name);
        exit (1);
      }
    map = Load_Track(db,"map");
  }

  //  Process read index arguments into a list of read ranges

  input_pts = 0;
  if (argc == 3)
    { if (argv[2][0] != LAST_READ_SYMBOL || argv[2][1] != '\0')
        { char *eptr, *fptr;
          int   b, e;

          b = strtol(argv[2],&eptr,10);
          if (eptr > argv[2] && b > 0)
            { if (*eptr == '-')
                { if (eptr[1] != LAST_READ_SYMBOL || eptr[2] != '\0')
                    { e = strtol(eptr+1,&fptr,10);
                      input_pts = (fptr <= eptr+1 || *fptr != '\0' || e <= 0);
                    }
                }
              else
                input_pts = (*eptr != '\0');
            }
          else
            input_pts = 1;
        }
    }

  if (input_pts)
    { input = Fopen(argv[2],"r");
      if (input == NULL)
        exit (1);

      iter = init_file_iterator(input);
    }
  else
    { pts  = (int *) Malloc(sizeof(int)*2*(argc-1),"Allocating read parameters");
      if (pts == NULL)
        exit (1);

      reps = 0;
      if (argc > 2)
        { int   c, b, e;
          char *eptr, *fptr;

          for (c = 2; c < argc; c++)
            { if (argv[c][0] == LAST_READ_SYMBOL)
                { b = db->nreads;
                  eptr = argv[c]+1;
                }
              else
                b = strtol(argv[c],&eptr,10);
              if (eptr > argv[c])
                { if (b <= 0)
                    { fprintf(stderr,"%s: %d is not a valid index\n",Prog_Name,b);
                      exit (1);
                    }
                  if (*eptr == 0)
                    { pts[reps++] = b;
                      pts[reps++] = b;
                      continue;
                    }
                  else if (*eptr == '-')
                    { if (eptr[1] == LAST_READ_SYMBOL)
                        { e = db->nreads;
                          fptr = eptr+2;
                        }
                      else
                        e = strtol(eptr+1,&fptr,10);
                      if (fptr > eptr+1 && *fptr == 0 && e > 0)
                        { pts[reps++] = b;
                          pts[reps++] = e;
                          if (b > e)
                            { fprintf(stderr,"%s: Empty range '%s'\n",Prog_Name,argv[c]);
                              exit (1);
                            }
                          continue;
                        }
                    }
                }
              fprintf(stderr,"%s: argument '%s' is not an integer range\n",Prog_Name,argv[c]);
              exit (1);
            }
        }
      else
        { pts[reps++] = 1;
          pts[reps++] = db->nreads;
        }
    }

  //  Display each read (and/or QV streams) in the active DB according to the
  //    range pairs in pts[0..reps) and according to the display options.

  { int         c, b, e, i;
    int64      *anno;
    int        *data;
    int64       s, f, j;

    anno = (int64 *) map->anno;
    data =   (int *) map->data;

    c = 0;
    while (1)
      { if (input_pts)
          { if (next_read(iter))
              break;
            e = iter->read;
            b = e-1;
          }
        else
          { if (c >= reps)
              break;
            b = pts[c]-1;
            e = pts[c+1];
            if (e > db->nreads)
              e = db->nreads;
            c += 2;
          }

        if (PRETTY)
          for (i = b; i < e; i++)
            { s = (anno[i] >> 2);
              f = (anno[i+1] >> 2);
              printf(" %d -> %d(%d)",i+1,data[s]+1,data[s+1]);
              for (j = s+2; j < f; j += 3)
                { printf(" [%d,%d]",data[j],data[j+1]);
                  if (j+2 < f)
                    printf(" %d",data[j+2]);
                }
              printf("\n");
            }
        else
          for (i = b; i < e; i++)
            { s = (anno[i] >> 2);
              f = (anno[i+1] >> 2);
              printf(" %d %d %d %lld",i+1,data[s]+1,data[s+1],(f-s)-2);
              for (j = s+2; j < f; j += 3)
                { printf(" %d %d",data[j],data[j+1]);
                  if (j+2 < f)
                    printf(" %d",data[j+2]);
                }
              printf("\n");
            }
      }
  }

  if (input_pts)
    { fclose(input);
      free(iter);
    }
  else
    free(pts);

  Close_DB(db);

  exit (0);
}
