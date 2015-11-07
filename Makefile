CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

ALL = DASqv DBpaln DBsats

all: $(ALL)

DASqv: DASqv.c align.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o DASqv DASqv.c align.c DB.c QV.c -lm

DBpaln: DBpaln.c seqpats.c seqpats.h align.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o DBpaln DBpaln.c seqpats.c align.c DB.c QV.c

DBsats: DBsats.c seqpats.c seqpats.h align.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o DBsats DBsats.c seqpats.c align.c DB.c QV.c

clean:
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f scrubber.tar.gz

install:
	cp $(ALL) ~/bin

package:
	make clean
	tar -zcf scrubber.tar.gz README *.h *.c Makefile
