CC = gcc

#PROFILING=-g -pg

WARN=-Wall
CFLAGS=-O3 $(WARN) $(PROFILING)
LIBS=-lm
DEPS=fid.h

all: tester-avx tester-sse tester

#%.o : %.c $(DEPS)
#	$(CC) $(CFLAGS) $(WARN) -o $@ $<

utils.o : utils.c $(DEPS)
	$(CC) -c $(WARN) $(CFLAGS) -o $@ $<

tester: fid1.c utils.o
	$(CC) $(CFLAGS) $(WARN) $+ $(LIBS) -o $@

tester-sse: fid1_sse.c utils.o
	$(CC) $(CFLAGS) -mssse3 $(WARN) $+ $(LIBS) -o $@

tester-avx: fid1_avx.c utils.o
	$(CC) $(CFLAGS) -mavx $(WARN) $+ $(LIBS) -o $@
clean:
	rm tester-sse tester-avx tester utils.o

