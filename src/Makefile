CC = gcc

#PROFILING=-g -pg

WARN=-Wall
CFLAGS=-O3 $(WARN) $(PROFILING)
DEPS=fid.h

all: fid1

fid1: fid1.c utils.c
	$(CC) -Wall -O3 $+ -o $@
