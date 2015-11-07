/*******************************************************************************************
 *
 *  Sequence pattern finder module:
 *     Routines that check individual reads for interesting patterns such as micro-satellities,
 *       palindromes, flip patterns, etc.
 *
 *  Author:  Gene Myers
 *  Date  :  October 2014
 *
 ********************************************************************************************/

#ifndef _PS_MODULE

#define _PS_MODULE

//  Each of the pattern finders below operates on the i'th read of the specified db.


//  Flip_Finder searches for a suffix or prefix of db[i] that is a palindrome with an interior
//    gap that is not too big.  This pattern is characteristic of either
//           (a) a missed adapter call,
//           (b) the polymerase jumps from one strand to another (real?)
//        or (c) truly palindromic sequence at the end or beginning of the read
//
//   Flip_Finder always returns a pointer to an array of *nhit Flip_Hit records, in order
//   of division points.

typedef struct
  { int divpt;      // The imputed mid-point of the palindromic transition
    int qvmat;      // Alignment % of matching palindromes
    int gap;        // A region of gap bases symmetrically around the mid-point does not align
    int qvgap;      // Alignment % of gapped central part (only def'd if gap > 0)
  } Flip_Hit;

Flip_Hit *Flip_Finder(HITS_DB *db, int i, int *nhit, int beg, int end);


//  Pal_Finder scans left to right looking for palindromes that may have a small interior
//    gap (e.g. 30bp or so).  This pattern may also recognize a polymerase flip / missing adapter
//    in cases where the interior gap of the flip is small.
//
//  Pal_Finder returns 0 if it does not find a palindrome and 1 otherwise.  If it finds
//    a palindrome then calling it again with the same read (i.e. same value of i) will cause
//    it to search for the next non-overlapping palindrome in the remainder of the sequence beyond
//    the last reported hit.  Thus one can easily deliver all the palindromes in a read by
//    repeatedly calling Palindrome Finder until it returns a 0.

typedef struct
  { int midpt;      // The imputed mid-point of the palindromic transition
    int qvmat;      // Alignment % of matching palindromes
    int gap;        // The center gap bases did not align well
    int width;      // Width of the palindrome
  } Pal_Hit;

Pal_Hit  *Pal_Finder(HITS_DB *db, int i);


//  Sat_Finder scans left to right looking for micro/mini-satellites up to length SatMax
//    (defined in seqpats.c).
//
//  Sat_Finder retuns NULL if no satellites were found and a pointer to an array of Sat_Hit
//    intervals ordered left to right.  The array is terminated by a record with 'unit = 0'


typedef struct
  { int beg;
    int end;
  } Sat_Hit;

Sat_Hit *Sat_Finder(HITS_DB *db, int i);

#endif // _PS_MODULE
