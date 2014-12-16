#FID (Fragment Insertions and Deletions)

## Introduction

The aim is to develop fast, vectorized functions for computing the FID model. The implementation should be:

* open source with an appropriate open source license
* 64-bit design supporting the most recent and popular SIMD architectures

## List of available implementations

As a start, we should implement two FID models, namely FID-1 and FID-2.

#### FID-1

* A column-based implementation that requires less memory, but backtracking is not possible as only the last column is stored.
* A full matrix implementation, that requires more memory (the size of the matrix) but allows backtracking.
* An SSE-3 and AVX based vectorized implementation of the above two methods.

#### FID-2

* A column-based implementation that requires less memory, but backtracking is not possible as only the last column is stored.
* A full matrix implementation, that requires more memory (the size of the matrix) but allows backtracking.
* An SSE-3 and AVX based vectorized implementation of the above two methods.

## Code

The code is written in C. 

    File     | Description
-------------|------
**fid1.c** | FID-1 implementations.
**utils.c** | Various common utility functions.

## Bugs

The source code has not been tested comprehensively yet. All bug reports are highly appreciated.
