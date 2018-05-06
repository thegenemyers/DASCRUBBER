
# Dascrubber: The Dazzler Read Scrubbing Suite

## _Author:  Gene Myers_
## _First:   March 27, 2016_

For typeset documentation, examples of use, and design philosophy please go to
my [blog](https://dazzlerblog.wordpress.com/command-guides/dascrubber-command-guide).

This is still a preliminary release.
The current set of commands provide a pipeline that one can use to scrub reads and if desired
to scrub the alignment piles (with DASrealign).  Ultimately DASpatch/DASedit and DASrealign will
be replaced with more powerful programs that correct reads and not only scrub alignment piles,
but also remove haplotype and repeat induced overlaps, prior to assembly via a string graph
method.

The goal of scrubbing is to produce a set of edited reads that are guaranteed to
(a) be continuous stretches of the underlying genome (i\.e\. no unremoved adapters
and not chimers), and (b) have no very low quality stretches (i\.e\. the error rate
never exceeds some reasonable maximum, 20% or so in the case of Pacbio data).  The
secondary goal of scrubbing is to do so with the minimum removal of data and splitting
of reads.

Note carefully that the current scrubbing pipeline requires that one has employed
repeat-masking in the daligner run as per the DAMASKER module described in this
[post](https://dazzlerblog.wordpress.com/2016/04/01/detecting-and-soft-masking-repeats).

The current \"DAS\" suite consists of a pipeline of several programs that in sequence accomplish
the task of scrubbing: DAScover &rarr; DASqv &rarr; DAStrim &rarr; DASpatch &rarr; DASedit.
For the commands, the \<source\> argument must always refer to the entire DB, and only
the \<overlaps\> arguments can involve a block number.  If \<overlaps\> involves a block
number, e\.g\. Ecoli.2.las, then the .las file is expected to contain all the overlaps
where the A-read is in block 2 of the underlying database.  The HPC.daligner scripts in the
DALIGNER module produce such .las files as their final result.  Parameters are propoagated
down the pipeline to subsequent phases via the annotation  tracks so one need not specify
the same parameter over and over again, i.e, the parameters/flags -H, -c, -g, and -b.

All programs add suffixes (e.g. .db, .las) as needed.
For the commands that take multiple .las block files as arguments, e.g. DAScover, DASqv, ...,
one can place a @-sign in the name, which is then interpreted as the sequence of files
obtained by replacing the @-sign by 1, 2, 3, ... in sequence until a number is reached for
which no file matches.  One can also place a @-sign followed by an integer, say, i, in which
case the sequence starts at i.  Lastly, one can also place @i-j where i and j are integers, in
which case the sequence is from i to j, inclusive.

```
1. DAScover [-v] [-H<int>] [-m<track>]+ <subject:db> <overlaps:las> ...
```

This command takes as input a database \<source\> and a sequence of sorted local alignments
blocks, \<overlaps\>, produced by an overlap/daligner run for said database.
Note carefully that \<source\> must always refer to the entire DB, only the \<overlaps\> can
involve block numbers.

Using the local
alignment-pile for each A-read, DAScover produces a histogram of the depth of coverage of
each trace point tile that is not within one of the intervals of the optionally specified tracks.
It places this histogram in a .covr track for the bock and these block tracks are merged
later with Catrack.  If the -v option is set, the histogram for each block is displayed and an
estimate of the coverage of the underlying target genome is output.

If the overlap file contains a block number
then the track files also contain a block number, e\.g\. \"DAScovr DB OVL.2\" will result in
the track files DB.2.covr.[anno,data].  Furthermore, if DAScovr is run on .las blocks, then
once it has been run on all the blocks of the DB, the block tracks must be concatenated
into a single track for the entire database with Catrack in preparation for the next phase
of scrubbing by DAStrim.

```
2. DASqv [-v] [-c<int>] <subject:db> <overlaps:las> ...
```

This command takes as input a database \<source\> and a sequence of sorted local alignments
blocks, \<overlaps\>, produced by an overlap/daligner run for said database.
A .covr track obtained by running DASqv and Catrack must be present for the entire data base.
Note carefully that \<source\> must always refer to the entire DB, only the \<overlaps\> can
involve a block number.

Using the local alignment-pile for each A-read, DASqv produces a QV value for each complete
segment of TRACE_SPACING bases (e\.g\. 100bp, the -s parameter to daligner). The quality value
of ecah trace tile is the average of the best 25-50% of the estimated coverage alignment matches,
where the estimated coverage is computed from the histogram of the .covr track.
If one supplies the -c parameter, than this value explicitly overrides the estimated coverage
produced by default.  All quality values over 50 are clipped to 50.

The -v option prints out a histogram of the segment align matches, and the quality values
produced.  This histogram is useful in assessing, for a given data set, what constitutes the
threshold -g and -b, to be used by down stream commands, for what is definitely a good segment
and what is definitely a bad segment.  The -H option is for HGAP-based assembly (see the -H
option of daligner) wherein only reads longer than the -H parameter are considered for overlap,
scrubbing, and assembly.  With this option set, DASqv and all subsequent commands in the
scrubbing pipeline, only perform their functions on reads of length -H or more.  All other
reads are used in the overlap piles for H-reads to help assess and scrub the H-reads, but
are themselves not scrubbed.

The quality values are written to a .qual track, that can be viewed by calling DBdump with
the -i option set (\"i\" for \"intrinsic QV\").
Like DAScovr and all other
scrubber modules, block tracks are produced in response to block .las files and these must be
concatenated with Catrack into a single .qual track for the entire DB in preparation for the
next phase of scrubbing by DAStrim.

```
3. DAStrim [-v] [-g<int>] [-b<int>] <subject:db> <overlaps:las> ...
```

A DB-wide .qual track produced by DASqv and Catrack are required as input to this command.
This command further takes as input a database \<source\> and a sequence of sorted local alignments
blocks, \<overlaps\>, produced by an overlap/daligner run for said database.
 A .qual track
obtained by running DASqv must be present for the entire data base.  Note carefully that
\<source\> must always refer to the entire DB, only \<overlaps\> can involve block numbers.

Using the local alignment-pile for each A-read and the QV\'s for all the reads in the pile,
DAStrim (1) finds and breaks all chimeric reads, (2) finds all missed adaptamers and retains
only the longest subread between missed adaptaers, and (3) identifies all low-quality regions
that should be improved/replaced by better sequence.  It makes these inherently heuristic
decisions conservatively so that what remains is very highly likely to be free of chimers,
adaptamers, and undetected low-quality sequence segments.  Some of these artifacts may still
get through, but at very low odds, less than 1 in 10,000 in our experience.  The decision
process is guided by the parameters -g and -b which indicate the thresholds for considering
intrinsic QV values good, bad, or unknown.  By default these parameters are automatically set
to be the 80'th and 93'rd percentiles of the qv-histograms hidden in the .qual track.  They
may however be explicitly set at the command line to over-rule this default choice.

The -v option prints out a report of how many chimer and adaptamer breaks were detected, how
much sequence was trimmed, how many low-quality segments were spanned by alignments, and how
many were rescued by many pairs of local alignments spanning the gap indicued by the
low-quality region, and so on.

The retained high-quality intervals for each read are written to a .trim track, in
left-to-right order with an indicator of whether the gap between two such intervals is spanned
by local alignments or by span-consistent pairs of local alignments.
Like DAScovr and all other
scrubber modules, block tracks are produced in response to block .las files and these must be
concatenated with Catrack into a single .trim track for the entire DB in preparation for the
next phase of scrubbing by DASpatch.

```
4. DASpatch [-v] <subject:db> <overlaps:las> ...
```

This command takes as input a database \<source\> and a sequence of sorted local alignments
blocks, \<overlaps\>, produced by an overlap/daligner run for said database.  A .qual track
and a .trim track obtained by running DASqv and DAStrim must be present for the entire data
base.  Note carefully that \<source\> must always refer to the entire DB, only \<overlaps\> can
involve block numbers.

Using the local alignment-pile for each A-read, the QV\'s for all the reads in the pile, and
the hiqh-quality segments annotated by DAStrim, DASpatch selects a high-quality B-read segment
with which to patch every intervening low-quality segment of an A-read.  Given that these gaps
are annotated before each read is trimmed by DAStrim, it may be the case in this second
examination, that the gap is no longer spanned by the now trimmed B-reads in which case the
span/patch can fail.  This is very rare but does occur and is the number of such events is
reported by DASpatch when the -v option is set.

The B-read segments for each patch (or a special \"failure patch\") are written to a .patch
track, in left-to-right order.  Like DAScovr and all other scrubber modules, block tracks are
produced in response to block .las files and these must be concatenated with Catrack into a
single .patch track for the entire DB in preparation for the next phase of scrubbing by DASedit.

```
5. DASedit [-v] [-x<int>] <source:db> <target:db>
```

This command takes as input a database \<source\> for which a .trim, and .patch tracks have
produced by DAStrim and DASpatch in sequence (and perforce DAScover and DASqv before them).
Using the
information in the two tracks, DASedit produces a new database \<target\> whose reads are the
patched, high-quality sub-reads.  That is, every low quality segment is patched with the relevant
sequence of another read, and some reads give rise to two or more reads if deemed chimers, or
no reads if the entire read was deemed junk.  This command can take considerable time as the
access pattern to read sequences (for the patching) is not sequential or localized, implying
poor cache performance.

The new database does not have a .qvs or .arr component, that is, it is a a sequence or
S-database (see the original Dazzler DB post).  Very importantly, the new database has exactly
the same block divisions as the original.  That is, all patched subreads in a block of the new
database have been derived from reads of the same block in the original database, and only
from those reads.  The new database does have a .map track that for each read encodes the
original read in the source DB that was patched and the segments of that read that were
high-quality (i\.e\. not patched).  The program DASmap below can be used to output this
information in either an easy-to-read or an easy-to-parse format.

```
6. DASmap [-p] <path:db> [ <reads:FILE> |<reads:range> ... ]
```

This command takes as input a database of patched reads \<source\> produced by DASedit, and for
the specified reads outputs a line for each showing the source read index and length in the
originating DB, as well as annotating which segments were original and which were patched.
The convention on interpreting the read arguments is as for DBshow and DBdump.  As an example,
with the -p option (pretty print) set one might see:

```
     55 -> 57(2946) [400,2946]
     56 -> 58(11256) [700,1900]
     57 -> 58(11256) [6600,9900] 83 [10000,11256]
     58 -> 59(12282) [400,4100] 88 [4200,9400] 97 [9500,12282]
```

The first line indicates that read 55 in the patched database was derived from read 57 in the
original database and is the segment from [400,2946] of that read.  Reads 56 and 57 were both
derived from read 58 in the original DB, and read 57 consists of segments [6600,9900] and
[10000,11256] of read 58 with a patch between them of 83bp (but the source of the patch data
is not given).  The read length of each original read is given for convenience.  With the -p
option off, the output consists of space separated integers encoding the same information where
the 4th field is always the number of integers in the segment description (always 3n+2 for
some n):

```
     55 57 2946 2 400 2946
     56 58 11256 2 700 1900
     57 58 11256 5 6600 9900 83 10000 11256
     58 59 12282 8 400 4100 88 4200 9400 97 9500 12282
```

```
7. DASrealign [-v] [-l<int:(800)>] <block1> <block2> <source:las> <target:las>
```

This command takes as input two blocks a patched database \<source\> created by DASedit, and the
original .las file \<overlap\> for the block pair.  That is the .las file that was produced for
the blocks when daligner was run on the original database blocks.  DASrealign then produces
the set of alignments (inferrable from those in the original) between the new patched reads,
placing them in the file \<target\>.  These new .las files can then be merged
with LAmerge to form block .las files for the new database.

The idea of this program is to avoid having to run daligner again on the patched reads, but
rather to simply refine the alignments already computed with respect to the new read set.  It
has the draw back that there is some small chance that there are previously undetected overlaps
between reads now that they are patched, and the patched trace point encoding, while stillable
to deliver alignments are no longer usable for quality estimation as the tracepoint spacing in
the A-read becomes irregular.  This is contrasted with the speedup of the new process which
is roughly 40X faster than the original daligner run and the collection of overlaps in the
output can be feed directly into a string graph construction process.

```
8. REPcover <subject:db> ...
```

This command takes as input a sequence of databases or blocks \<source\> and for each outputs
a histogram of the coverage of the unmasked portions of the reads in the source along with
a recommendation of the -c value with which to run DASqv.  The .covr track produced by
DAScover must be present for all sources referred to.  The command is a quick way get the
-v output of DAScover at any time after producing the coverage histograms and to get the
histogram and coverage estimate for the entire data base (as opposed to a block of the database.

```
9. REPqv <subject:db> ...
```

This command takes as input a sequence of databases or blocks \<source\> and for each outputs
a histogram of the intrinsic quality values of the reads in the source along with a
recommendation of the -g and -b values with which to run DAStrim.  The .qual track produced
by DASqv must be present for all sources referred to.  The command is a quick way to get
the -v output of DASqv at any time after producing the intrinsic quality values and to get
the histograms for the entire data base (as opposed to a block of the database).

```
10. REPtrim <subject:db> ...
```

This command takes as input a sequence of databases or blocks \<source\> and for each outputs
the scrubbing statistics for the source, i.e. the same report produced by DAStrim with
the -v option set.  The .trim track produced by DAStrim must be present for all sources
referred to.  The command is a quick way to get the -v output of DAStrim at any time after
the fact and to get the statistics for the entire data base (as opposed to a block of the
database).
