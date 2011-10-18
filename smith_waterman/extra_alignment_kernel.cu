/*
 * Copyright 1993-2007 NVIDIA Corporation.  All rights reserved.
 *
 * NOTICE TO USER:
 *
 * This source code is subject to NVIDIA ownership rights under U.S. and
 * international Copyright laws.  Users and possessors of this source code
 * are hereby granted a nonexclusive, royalty-free license to use this code
 * in individual and commercial software.
 *
 * NVIDIA MAKES NO REPRESENTATION ABOUT THE SUITABILITY OF THIS SOURCE
 * CODE FOR ANY PURPOSE.  IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR
 * IMPLIED WARRANTY OF ANY KIND.  NVIDIA DISCLAIMS ALL WARRANTIES WITH
 * REGARD TO THIS SOURCE CODE, INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY, NONINFRINGEMENT, AND FITNESS FOR A PARTICULAR PURPOSE.
 * IN NO EVENT SHALL NVIDIA BE LIABLE FOR ANY SPECIAL, INDIRECT, INCIDENTAL,
 * OR CONSEQUENTIAL DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
 * OF USE, DATA OR PROFITS,  WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
 * OR OTHER TORTIOUS ACTION,  ARISING OUT OF OR IN CONNECTION WITH THE USE
 * OR PERFORMANCE OF THIS SOURCE CODE.
 *
 * U.S. Government End Users.   This source code is a "commercial item" as
 * that term is defined at  48 C.F.R. 2.101 (OCT 1995), consisting  of
 * "commercial computer  software"  and "commercial computer software
 * documentation" as such terms are  used in 48 C.F.R. 12.212 (SEPT 1995)
 * and is provided to the U.S. Government only as a commercial end item.
 * Consistent with 48 C.F.R.12.212 and 48 C.F.R. 227.7202-1 through
 * 227.7202-4 (JUNE 1995), all U.S. Government End Users acquire the
 * source code with only those rights set forth herein.
 *
 * Any use of this source code in individual and commercial software must
 * include, in the user documentation and internal comments to the code,
 * the above Disclaimer and U.S. Government End Users Notice.
 */

/* alignment project which demonstrates the basics on how to setup a project 
 * example application.
 * Device code.
 */

#ifndef _ALIGNMENT_KERNEL_H_
#define _ALIGNMENT_KERNEL_H_

#include <stdio.h>

#define MAX_SHARED_MEM 1000;
//#define SDATA( index)      cutilBankChecker(sdata, index)

/**
 *\brief init_sequence
 *
 * Initialize the sequences in device space
 */
__global__ void init_sequence()
{

}


__device__ inline int score(char ch1, char ch2)
{
  return (ch1 == '-' || ch2 == '-')? 0 : 2;
}



__global__ void alignmentKernel( char* query, int query_size, char* seqbase, int* seqbase_idx, int* score_list )
{
  // each thread compute a single alignment from the database
  int idx = threadIdx.x + blockIdx.x * blockDim.x ;

  // the host will determine the number of threads in each block
  // by the buffer size and query size
  char* sequence = seqbase + seqbase_idx[idx];
  int seq_size = seqbase_idx[idx+1] - seqbase[idx];
  
  // store the alignment buffer in the shared memory
  // each alignment will use a single buffer
  // whose size is the same as the length of the query
  // TODO: resolve the bank conflict issue
  __shared__ int g_buffer[MAX_SHARED_MEM];
  for (int i=0; i<seq_size; g_buffer[threadIdx.x * seq_size + i] = 0, ++i);

  int* buffer = g_buffer[threadIdx.x * seq_size];
  int max_val = -1;
  for ( int i = 0; i < seq_size; ++i)
    {
      int prev = 0;
      for ( int j = 1; j < query_size; ++j )
  	{
  	  char ch1 = sequence[i];
  	  char ch2 = query[j];
  	  int vl = prev + score(ch1, '-');
  	  int vu = buffer[j] + score(ch1, '-');
  	  int vd = buffer[j-1] + score(ch1, ch2);
  	  buffer[j-1] = prev;
  	  prev = max(vl, max(vd, vu));
  	  max_val = max(max_val, prev);
  	}
    }
  score_list[idx] = max_val;
}



#endif // #ifndef _alignment_KERNEL_H_
