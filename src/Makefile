CC = gcc

#PROFILING=-g -pg

WARN=-Wall
CFLAGS=-O3 $(WARN) $(PROFILING)
LIBS=-lm
DEPS=fid.h

all: tester-avx tester-avx-fma tester-sse tester

#%.o : %.c $(DEPS)
#	$(CC) $(CFLAGS) $(WARN) -o $@ $<

utils.o : utils.c $(DEPS)
	$(CC) -c $(CFLAGS) -o $@ $<

tester: fid1.c utils.o
	$(CC) $(CFLAGS) $+ $(LIBS) -o $@

tester-sse: fid1_sse.c utils.o
	$(CC) $(CFLAGS) -mssse3 $+ $(LIBS) -o $@

tester-avx: fid1_avx.c utils.o
	$(CC) $(CFLAGS) -mavx $+ $(LIBS) -o $@

tester-avx-fma: fid1_avx.c utils.o
	$(CC) $(CFLAGS) -DHAVE_FMA -mavx -mfma $+ $(LIBS) -o $@
clean:
	rm tester-sse tester-avx tester-avx-fma tester utils.o

