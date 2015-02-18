CC = gcc

#PROFILING=-g -pg

WARN=-Wall
CFLAGS=-O3 -g $(WARN) $(PROFILING)
LIBS=-lm
DEPS=fid.h

all: tester-avx2 tester-avx tester-avx-fma tester-sse tester-sse-fma tester tester-avx-fmath tester-fmath tester-fid1-emit tester-fid1-diag tester-fid1-diag-scale tester-fid1-diag-sse3

#%.o : %.c $(DEPS)
#	$(CC) $(CFLAGS) $(WARN) -o $@ $<

utils.o : utils.c $(DEPS)
	$(CC) -c $(CFLAGS) -o $@ $<

tester: fid1.c utils.o
	$(CC) $(CFLAGS) $+ $(LIBS) -o $@

tester-fmath: fid1.c utils.o
	$(CC) $(CFLAGS) -ffast-math $+ $(LIBS) -o $@

tester-sse: fid1_sse.c utils.o
	$(CC) $(CFLAGS) -mssse3 $+ $(LIBS) -o $@

tester-sse-fma: fid1_sse.c utils.o
	$(CC) $(CFLAGS) -DHAVE_FMA -mssse3 -mfma $+ $(LIBS) -o $@

tester-avx: fid1_avx.c utils.o
	$(CC) $(CFLAGS) -mavx $+ $(LIBS) -o $@

tester-avx-fmath: fid1_avx.c utils.o
	$(CC) $(CFLAGS) -ffast-math -mavx $+ $(LIBS) -o $@

tester-avx2: fid1_avx.c utils.o
	$(CC) $(CFLAGS) -DHAVE_AVX2 -mavx2 $+ $(LIBS) -o $@

tester-fid1-emit: fid1_emit.c utils.o
	$(CC) $(CFLAGS) $+ $(LIBS) -o $@

tester-fid1-diag: fid1_diagonal.c utils.o
	$(CC) $(CFLAGS) $+ $(LIBS) -o $@

tester-fid1-diag-scale: fid1_diagonal.c utils.o
	$(CC) -DSCALING $(CFLAGS) $+ $(LIBS) -o $@

tester-fid1-diag-sse3: fid1_diagonal_sse.c utils.o
	$(CC) -DSCALING -mssse3 $(CFLAGS) $+ $(LIBS) -o $@

tester-avx-fma: fid1_avx.c utils.o
	$(CC) $(CFLAGS) -DHAVE_FMA -mavx -mfma $+ $(LIBS) -o $@
clean:
	rm tester-sse tester-sse-fma tester-avx tester-avx-fma tester-avx2 tester utils.o tester-fmath tester-fid1-emit tester-fid1-diag-scale tester-fid1-diag-sse3

