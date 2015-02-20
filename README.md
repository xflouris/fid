#FID (Fragment Insertions and Deletions)

## Introduction

The aim is to develop fast, vectorized functions for computing the FID model. The implementation should be:

* open source with an appropriate open source license
* 64-bit design supporting the most recent and popular SIMD architectures

## List of available implementations

As a start, we should implement two FID models, namely FID-1 and FID-2.

#### FID-1

* A diagonal-based double-precision implementation with scaling, allowing both emission probabilities and homologies.
* A diagonal-based double-precision implementation with scaling, disallowing both emission probabilities and homologies.
* A diagonal-based double-precision implementation with scaling, considering emission probabilities and disallowing homologies.
* A SSE-3 vectorized diagonal-based double-precision implementation with scaling, allowing emission probabilities and homologies.
* A SSE-3 vectorized diagonal-based double-precision implementation with scaling, disallowing both emission probabilities and homologies.
* A SSE-3 vectorized diagonal-based double-precision implementation with scaling, considering emission probabilities and disallowing homologies.

## Code

The code is written in C. 

    File     | Description
-------------|------
**fid.h** | Header definitions
**fid1_diagonal.c** | Diagonal double-precision implementation with scaling, allowing emissions and homologies.
**fid1_diagonal_sse.c** | SSE-3 Vectorized diagonal double-precision implementation with scaling, allowing emissions and homologies.

## Testing

Test scripts are located in the /test directory. Run as

`./test.scale.sh > scale.txt`
`./test.sse3.sh > sse3.txt`

This will produce all possible instances of matrices with the available data and output all outputs to a file. The two files can then be compared with:

`diff scale.txt sse3.txt`

## Bugs

The source code has not been tested comprehensively yet. All bug reports are highly appreciated.
