DEST_DIR = ~/bin

CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

ALL = DASqv DAStrim

all: $(ALL)

DASqv: DASqv.c align.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o DASqv DASqv.c align.c DB.c QV.c -lm

DAStrim: DAStrim.c align.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o DAStrim DAStrim.c align.c DB.c QV.c -lm

clean:
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f scrubber.tar.gz

install:
	cp $(ALL) $(DEST_DIR)

package:
	make clean
	tar -zcf scrubber.tar.gz README Makefile *.h *.c
