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

/* Template project which demonstrates the basics on how to setup a project 
 * example application.
 * Host code.
 */

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
//#include <cutil_inline.h>

// includes, kernels
//#include <alignment_kernel.cu>

//********************** global variables ***************************
// BLOSUM62 matrix and its size
// these should be copied to constant memory in the device
/* const int size_blosum  = 24; */
/* int BLOSUM62[size_blosum*size_blosum] = { */
/*   4,  -1,  -2,  -2,   0,  -1,  -1,   0,  -2,  -1,  -1,  -1,  -1,  -2,  -1,   1,   0,  -3,  -2,   0,  -2,  -1,   0,  -4,    */
/*   -1,   5,   0,  -2,  -3,   1,   0,  -2,   0,  -3,  -2,   2,  -1,  -3,  -2,  -1,  -1,  -3,  -2,  -3,  -1,   0,  -1,  -4,    */
/*   -2,   0,   6,   1,  -3,   0,   0,   0,   1,  -3,  -3,   0,  -2,  -3,  -2,   1,   0,  -4,  -2,  -3,   3,   0,  -1,  -4,    */
/*   -2,  -2,   1,   6,  -3,   0,   2,  -1,  -1,  -3,  -4,  -1,  -3,  -3,  -1,   0,  -1,  -4,  -3,  -3,   4,   1,  -1,  -4,    */
/*   0,  -3,  -3,  -3,   9,  -3,  -4,  -3,  -3,  -1,  -1,  -3,  -1,  -2,  -3,  -1,  -1,  -2,  -2,  -1,  -3,  -3,  -2,  -4,    */
/*   -1,   1,   0,   0,  -3,   5,   2,  -2,   0,  -3,  -2,   1,   0,  -3,  -1,   0,  -1,  -2,  -1,  -2,   0,   3,  -1,  -4,    */
/*   -1,   0,   0,   2,  -4,   2,   5,  -2,   0,  -3,  -3,   1,  -2,  -3,  -1,   0,  -1,  -3,  -2,  -2,   1,   4,  -1,  -4,    */
/*   0,  -2,   0,  -1,  -3,  -2,  -2,   6,  -2,  -4,  -4,  -2,  -3,  -3,  -2,   0,  -2,  -2,  -3,  -3,  -1,  -2,  -1,  -4,    */
/*   -2,   0,   1,  -1,  -3,   0,   0,  -2,   8,  -3,  -3,  -1,  -2,  -1,  -2,  -1,  -2,  -2,   2,  -3,   0,   0,  -1,  -4,    */
/*   -1,  -3,  -3,  -3,  -1,  -3,  -3,  -4,  -3,   4,   2,  -3,   1,   0,  -3,  -2,  -1,  -3,  -1,   3,  -3,  -3,  -1,  -4,    */
/*   -1,  -2,  -3,  -4,  -1,  -2,  -3,  -4,  -3,   2,   4,  -2,   2,   0,  -3,  -2,  -1,  -2,  -1,   1,  -4,  -3,  -1,  -4,    */
/*   -1,   2,   0,  -1,  -3,   1,   1,  -2,  -1,  -3,  -2,   5,  -1,  -3,  -1,   0,  -1,  -3,  -2,  -2,   0,   1,  -1,  -4,    */
/*   -1,  -1,  -2,  -3,  -1,   0,  -2,  -3,  -2,   1,   2,  -1,   5,   0,  -2,  -1,  -1,  -1,  -1,   1,  -3,  -1,  -1,  -4,    */
/*   -2,  -3,  -3,  -3,  -2,  -3,  -3,  -3,  -1,   0,   0,  -3,   0,   6,  -4,  -2,  -2,   1,   3,  -1,  -3,  -3,  -1,  -4,    */
/*   -1,  -2,  -2,  -1,  -3,  -1,  -1,  -2,  -2,  -3,  -3,  -1,  -2,  -4,   7,  -1,  -1,  -4,  -3,  -2,  -2,  -1,  -2,  -4,    */
/*   1,  -1,   1,   0,  -1,   0,   0,   0,  -1,  -2,  -2,   0,  -1,  -2,  -1,   4,   1,  -3,  -2,  -2,   0,   0,   0,  -4,    */
/*   0,  -1,   0,  -1,  -1,  -1,  -1,  -2,  -2,  -1,  -1,  -1,  -1,  -2,  -1,   1,   5,  -2,  -2,   0,  -1,  -1,   0,  -4,    */
/*   -3,  -3,  -4,  -4,  -2,  -2,  -3,  -2,  -2,  -3,  -2,  -3,  -1,   1,  -4,  -3,  -2, 11,   2,  -3,  -4,  -3,  -2,  -4,    */
/*   -2,  -2,  -2,  -3,  -2,  -1,  -2,  -3,   2,  -1,  -1,  -2,  -1,   3,  -3,  -2,  -2,   2,   7,  -1,  -3,  -2,  -1,  -4,    */
/*   0,  -3,  -3,  -3,  -1,  -2,  -2,  -3,  -3,   3,   1,  -2,   1,  -1,  -2,  -2,   0,  -3,  -1,   4,  -3,  -2,  -1,  -4,    */
/*   -2,  -1,   3,   4,  -3,   0,   1,  -1,   0,  -3,  -4,   0,  -3,  -3,  -2,   0,  -1,  -4,  -3,  -3,   4,   1,  -1,  -4,    */
/*   -1,   0,   0,   1,  -3,   3,   4,  -2,   0,  -3,  -3,   1,  -1,  -3,  -1,   0,  -1,  -3,  -2,  -2,   1,   4,  -1,  -4,    */
/*   0,  -1,  -1,  -1,  -2,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -2,   0,   0,  -2,  -1,  -1,  -1,  -1,  -1,  -4,    */
/*   -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,   1 */
/* }; */

////////////////////////////////////////////////////////////////////////////////
// declaration, forward
void alignment();

// external definition of the alignment kernel
extern "C"
void alignmentKernel( int* g_idata, int* g_index, int* g_odata );



////////////////////////////////////////////////////////////////////////////////
//! Run a simple test for CUDA
////////////////////////////////////////////////////////////////////////////////
//void runTest( int argc, char** argv) 
void alignment()
{
  /* // use command-line specified CUDA device, otherwise use device with highest Gflops/s */
  /* if( cutCheckCmdLineFlag(argc, (const char**)argv, "device") ) */
  /*   cutilDeviceInit(argc, argv); */
  /* else */
  /*   cudaSetDevice( cutGetMaxGflopsDeviceId() ); */

  /* // setup the timer */
  /* unsigned int timer = 0;  */
  /* cutilCheckError( cutCreateTimer( &timer)); */
  /* cutilCheckError( cutStartTimer( timer)); */

  /* // input */
  /* char *a = "atcggtcaattggcc"; */
  /* char *b = "atgccctgatgggccga"; */

  /* // dimension of multi-threading */
  /* unsigned int num_threads = 1; */
  /* unsigned int mem_size = sizeof( char ) * (num_threads + 1) * 1000; */
  /* unsigned int index_size = (num_threads + 2) * sizeof(int); */
  /* unsigned int out_size = num_threads * sizeof(int); */
  /* unsigned int shared_mem_size =  num_threads * strlen(a) * sizeof(char) * 2 ; */

  /* // allocate host memory */
  /* char* h_idata = (char*) malloc( mem_size ); */
  /* int* h_index = (int*) malloc( index_size ); */
  /* // initalize the memory */
  /* /\* for( unsigned int i = 0; i < num_threads; ++i)  */
  /*    { */
  /*    h_idata[i] = (float) i; */
  /*    }*\/ */
	
  /* *(h_index) = strlen(a); */
  /* *(h_index+1) = strlen(a) + strlen(b); */
  /* memcpy( h_idata, a, strlen(a) * sizeof(char) ); */
  /* memcpy( h_idata + strlen(a) * sizeof(char), b, strlen(b) * sizeof(char) ); */

  /* // allocate device memory */
  /* int* d_idata; */
  /* int* d_index; */
  /* cutilSafeCall( cudaMalloc( (void**) &d_idata, mem_size)); */
  /* cutilSafeCall( cudaMalloc( (void**) &d_index, index_size ) ); */
  /* // copy host memory to device */
  /* cutilSafeCall( cudaMemcpy( d_idata, h_idata, mem_size, cudaMemcpyHostToDevice) ); */
  /* cutilSafeCall( cudaMemcpy( d_index, h_index, index_size, cudaMemcpyHostToDevice ) ); */

  /* // allocate device memory for result */
  /* int* d_odata; */
  /* cutilSafeCall( cudaMalloc( (void**) &d_odata, out_size ) ); */

  /* // setup execution parameters */
  /* dim3  grid( 1, 1, 1); */
  /* dim3  threads( num_threads, 1, 1); */

  /* // execute the kernel */
  /* alignmentKernel<<< grid, threads, shared_mem_size >>>( d_idata, d_index,d_odata); */

  /* // check if kernel execution generated and error */
  /* cutilCheckMsg("Kernel execution failed"); */

  /* // allocate mem for the result on host side */
  /* int* h_odata = (int*) malloc( out_size ); */
  /* // copy result from device to host */
  /* cutilSafeCall( cudaMemcpy( h_odata, d_odata, sizeof( float) * num_threads, */
  /* 			     cudaMemcpyDeviceToHost) ); */

  /* // compute the time spent for device computation */
  /* cutilCheckError( cutStopTimer( timer)); */
  /* printf( "Processing time: %f (ms)\n", cutGetTimerValue( timer)); */
  /* cutilCheckError( cutDeleteTimer( timer)); */

  /* // calculate the final result */
  /* int result = 0; */
  /* for ( int i =0; i< num_threads; ++i ) */
  /*   { */
  /*     result += *( h_odata + i ); */
  /*   } */

  /* printf("result is %d\n", result); */

  /* // cleanup memory */
  /* free( h_idata); */
  /* free( h_index ); */
  /* free( h_odata); */
  /* cutilSafeCall(cudaFree(d_idata)); */
  /* cutilSafeCall(cudaFree(d_odata)); */
  /* cutilSafeCall(cudaFree(d_index)); */

  /* cudaThreadExit(); */
}
