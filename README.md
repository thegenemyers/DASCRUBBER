
# Dascrubber: The Dazzler Read Scrubbing Suite

## _Author:  Gene Myers_
## _First:   March 27, 2016_

For typeset documentation, examples of use, and design philosophy please go to
my [blog](https://dazzlerblog.wordpress.com/command-guides/dascrubber-command-guide).

This is still an incomplete release.
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
of reads.  The \"DAS\" suite consists of a pipeline of several programs that in sequence accomplish
the task of scrubbing.

```
1. DASqv [-v] [-H<int>] -c<int> <subject:db> <overlaps:las> ...
```

This command takes as input a database \<source\> and a sequence of sorted local alignments
blocks, \<overlaps\>, produced by an overlap/daligner run for said database.  It is recommended
that one employ repeat-masking in the overlap run as per the new DAMASKER module described in
this post.  Note carefully that \<source\> must always refer to the entire DB, only \<overlaps\>
can involve a block number.  If \<overlaps\> involves a block number, e\.g\. Ecoli.2.las, then
the .las file is expected to contain all the overlaps where the A-read is in block 2 of the
underlying database.  The HPC.daligner scripts in the DALIGNER module produce such .las files
as their final result.

Using the local alignment-pile for each A-read, DASqv produces a QV value for each complete
segment of TRACE_SPACING bases (e\.g\. 100bp, the -s parameter to daligner). The quality value
of the average percentile of the best 25-50% alignment matches covering it depending on the
coverage estimate -c.  One must supply the -c parameter to the expected coverage of the genome
in question.  All quality values over 50 are clipped to 50.

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
the -i option set (\"i\" for \"intrinsic QV\").  If the overlap file contains a block number
then the track files also contain a block number, e\.g\. \"DASqv -c50 DB OVL.2\" will result in
the track files DB.2.qual.[anno,data].  Furthermore, if DASqv is run on .las blocks, then
once it has been run on all the blocks of the DB, the block tracks must be concatenated
into a single track for the entire database with Catrack in preparation for the next phase
of scrubbing by DAStrim.

```
2. DAStrim [-v] -g<int> -b<int> <subject:db> <overlaps:las> ...
```

This command takes as input a database \<source\> and a sequence of sorted local alignments
blocks, \<overlaps\>, produced by an overlap/daligner run for said database.  A .qual track
obtained by running DASqv must be present for the entire data base.  Note carefully that
\<source\> must always refer to the entire DB, only \<overlaps\> can involve a block number.

Using the local alignment-pile for each A-read and the QV\'s for all the reads in the pile,
DAStrim (1) finds and breaks all chimeric reads, (2) finds all missed adaptamers and retains
only the longest subread between missed adaptaers, and (3) identifies all low-quality regions
that should be improved/replaced by better sequence.  It makes these inherently heuristic
decisions conservatively so that what remains is very highly likely to be free of chimers,
adaptamers, and undetected low-quality sequence segments.  Some of these artifacts may still
get through, but at very low odds, less than 1 in 10,000 in our experience.  The decision
process is guided by the parameters -g and -b which indicate the thresholds for considering
intrinsic QV values good, bad, or unknown.  In our experience, -g should be set to the value
for which 80% or more of the QV\'s produced by DASqv are better than this value, and -b should
be set to the value for which 5-7% or more of the QV\'s are worse than this value.  Consult
the histogram produced by DASqv to set these parameter.  We further recommend doing this
for each project on a case-by-case basis.

The -v option prints out a report of how many chimer and adaptamer breaks were detected, how
much sequence was trimmed, how many low-quality segments were spanned by alignments, and how
many were rescued by many pairs of local alignments spanning the gap indicued by the
low-quality region, and so on.

The retained high-quality intervals for each read are written to a .trim track, in
left-to-right order with an indicator of whether the gap between two such intervals is spanned
by local alignments or by span-consistent pairs of local alignments.  Like DASqv and all other
scrubber modules, block tracks are produced in response to block .las files and these must be
concatenated with Catrack into a single .trim track for the entire DB in preparation for the
next phase of scrubbing by DASpatch.

```
3. DASpatch [-v] <subject:db> <overlaps:las> ...
```

This command takes as input a database \<source\> and a sequence of sorted local alignments
blocks, \<overlaps\>, produced by an overlap/daligner run for said database.  A .qual track
and a .trim track obtained by running DASqv and DAStrim must be present for the entire data
base.  Note carefully that \<source\> must always refer to the entire DB, only \<overlaps\> can
involve a block number.

Using the local alignment-pile for each A-read, the QV\'s for all the reads in the pile, and
the hiqh-quality segments annotated by DAStrim, DASpatch selects a high-quality B-read segment
with which to patch every intervening low-quality segment of an A-read.  Given that these gaps
are annotated before each read is trimmed by DAStrim, it may be the case in this second
examination, that the gap is no longer spanned by the now trimmed B-reads in which case the
span/patch can fail.  This is very rare but does occur and is the number of such events is
reported by DASpatch when the -v option is set.

The B-read segments for each patch (or a special \"failure patch\") are written to a .patch
track, in left-to-right order.  Like DASqv and all other scrubber modules, block tracks are
produced in response to block .las files and these must be concatenated with Catrack into a
single .patch track for the entire DB in preparation for the next phase of scrubbing by DASedit.

```
4. DASedit [-v] [-x<int>] <source:db> <target:db>
```

This command takes as input a database \<source\> for which a .trim, and .patch tracks have
produced by DAStrim and DASpatch in sequence (and perforce DASqv initially).  Using the
information in the two tracks, DASedit produces a new database \<target\> whose reads are the
patched, high-quality sub-reads.  That is, every low quality segment is patched with the relevant
sequence of another read, and some reads give rise to two or more reads if deemed chimers, or
no reads if the entire read was deemed junk.  This command can take considerable time as the
access pattern to read sequences (for the patching) is not sequential or localized, implying
poor cache performance.

The new database does not have a .qvs or .arr component, that is it is a a sequence or
S-database (see the original Dazzler DB post).  Very importantly, the new database has exactly
the same block divisions as the original.  That is, all patched subreads in a block of the new
database have been derived from reads of the same block in the original database, and only
from those reads.  The new database does have a .map track that for each read encodes the
original read in the source DB that was patched and the segments of that read that were
high-quality (i\.e\. not patched).  The program DASmap below can be used to output this
information in either an easy-to-read or an easy-to-parse format.

```
5. DASmap [-p] <path:db> [ <reads:FILE> |<reads:range> ... ]
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
6. DASrealign [-v] [-l<int:(800)>] <block1> <block2> <source:las> <target:las>
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
7. REPqv <subject:db> ...
```

This command takes as input a sequence of databases or blocks \<source\> and for each outputs
a histogram of the intrinsic quality values of the reads in the source along with a
recommendation of the -g and -b values with which to run DAStrim.  The .qual track produced
by DASqv must be present for all sources referred to.  The command is a quick way to get
the -v output of DASqv at any time after producing the intrinsic quality values and to get
the histogram for the entire data base (as opposed to a block of the database).

```
8. REPtrim <subject:db> ...
```

This command takes as input a sequence of databases or blocks \<source\> and for each outputs
the scrubbing statistics for the source, i.e. the same report produced by DAStrim with
the -v option set.  The .trim track produced by DAStrim must be present for all sources
referred to.  The command is a quick way to get the -v output of DAStrim at any time after
the fact and to get the statistics for the entire data base (as opposed to a block of the
database).
