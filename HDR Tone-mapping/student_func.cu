/* Udacity Homework 3
   HDR Tone-mapping

  Background HDR
  ==============

  A High Dynamic Range (HDR) image contains a wider variation of intensity
  and color than is allowed by the RGB format with 1 byte per channel that we
  have used in the previous assignment.

  To store this extra information we use single precision floating point for
  each channel.  This allows for an extremely wide range of intensity values.

  In the image for this assignment, the inside of church with light coming in
  through stained glass windows, the raw input floating point values for the
  channels range from 0 to 275.  But the mean is .41 and 98% of the values are
  less than 3!  This means that certain areas (the windows) are extremely bright
  compared to everywhere else.  If we linearly map this [0-275] range into the
  [0-255] range that we have been using then most values will be mapped to zero!
  The only thing we will be able to see are the very brightest areas - the
  windows - everything else will appear pitch black.

  The problem is that although we have cameras capable of recording the wide
  range of intensity that exists in the real world our monitors are not capable
  of displaying them.  Our eyes are also quite capable of observing a much wider
  range of intensities than our image formats / monitors are capable of
  displaying.

  Tone-mapping is a process that transforms the intensities in the image so that
  the brightest values aren't nearly so far away from the mean.  That way when
  we transform the values into [0-255] we can actually see the entire image.
  There are many ways to perform this process and it is as much an art as a
  science - there is no single "right" answer.  In this homework we will
  implement one possible technique.

  Background Chrominance-Luminance
  ================================

  The RGB space that we have been using to represent images can be thought of as
  one possible set of axes spanning a three dimensional space of color.  We
  sometimes choose other axes to represent this space because they make certain
  operations more convenient.

  Another possible way of representing a color image is to separate the color
  information (chromaticity) from the brightness information.  There are
  multiple different methods for doing this - a common one during the analog
  television days was known as Chrominance-Luminance or YUV.

  We choose to represent the image in this way so that we can remap only the
  intensity channel and then recombine the new intensity values with the color
  information to form the final image.

  Old TV signals used to be transmitted in this way so that black & white
  televisions could display the luminance channel while color televisions would
  display all three of the channels.


  Tone-mapping
  ============

  In this assignment we are going to transform the luminance channel (actually
  the log of the luminance, but this is unimportant for the parts of the
  algorithm that you will be implementing) by compressing its range to [0, 1].
  To do this we need the cumulative distribution of the luminance values.

  Example
  -------

  input : [2 4 3 3 1 7 4 5 7 0 9 4 3 2]
  min / max / range: 0 / 9 / 9

  histo with 3 bins: [4 7 3]

  cdf : [4 11 14]


  Your task is to calculate this cumulative distribution by following these
  steps.

*/

#include "utils.h"
#include <iostream>

// Kernel function used for histogram
__global__ void kernel_histogram(const float* d_logLuminance, float lunMin, float lunRange, int numBins, int *bins, int size) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= size) return ;
  int binIdx = min(numBins - 1, static_cast<unsigned int>((d_logLuminance[i] - lunMin) / lunRange * numBins));
  atomicAdd(&bins[binIdx], 1);
}

// Kernel function used for scan operation
__global__ void kernel_upsweep_op(int len, unsigned int* input, unsigned int* output, int twod1) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = (i + 1) * twod1 - 1;
    if (idx < len) {
        output[idx] += output[idx - (twod1 >> 1)];
    }
}

// Kernel function used for scan operation
__global__ void kernel_downsweep_op(int len, unsigned int* input, unsigned int* output, int twod1) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = (i + 1) * twod1 - 1;
    if (idx < len) {
        int t = output[idx - (twod1 >> 1)];
        output[idx - (twod1 >> 1)] = output[idx];
        output[idx] += t;
    }
}

// __global__ void kernel_reduce_min(const float *d_logLuminance, float *mina, int gap, int size) {
//   extern __shared__ share_mem[];
//   const int tid = threadIdx.x;
//   int floorPow2 = blockDim.x;

//   if (tid + blockIdx.x * blockDim.x >= size) return ;

//   share_mem[tid] = d_logLuminance[tid + blockIdx.x * blockDim.x]

//   if (floorPow2 & (floorPow2 - 1)) {
//     while (floorPow2 & (floorPow2 - 1)) {
//       floorPow2 &= floorPow2 - 1;
//     }

//     if (tid >= floorPow2) {
//       share_mem[tid - floorPow2] = fmin(share_mem[tid - floorPow2], share_mem[tid]);
//       __syncthreads();
//     }
//   }

//   for (int activeThread = floorPow2 >> 1; activeThread; activeThread >>= 1) {
//     if (tid < activeThread) {
//       share_mem[tid] = fmin(share_mem[tid], share_mem[tid+activeThread]);
//       __syncthreads();
//     }
//   }
//   if (tid == 0) {
//     mina[tid] = share_mem[0];
//   }
// }

__global__ void kernel_reduce_max(const float *d_logLuminance, float *maxa, int gap, int size) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= size) return ;

  if (i % (2 * gap) == 0) {
    maxa[i] = fmax(d_logLuminance[i], d_logLuminance[i+gap]);
  }

  // const int idx = i * 2 * gap;
  // if (idx < size && idx + gap < size) {
  //   max[idx] = fmax(d_logLuminance[idx], d_logLuminance[idx+gap]);
  // }
}

static inline int nextPow2(int n)
{
    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    n++;
    return n;
}

void your_histogram_and_prefixsum(const float* const d_logLuminance,
                                  unsigned int* const d_cdf,
                                  float &min_logLum,
                                  float &max_logLum,
                                  const size_t numRows,
                                  const size_t numCols,
                                  const size_t numBins)
{
  //TODO
  /*Here are the steps you need to implement
    1) find the minimum and maximum value in the input logLuminance channel
       store in min_logLum and max_logLum
    2) subtract them to find the range
    3) generate a histogram of all the values in the logLuminance channel using
       the formula: bin = (lum[i] - lumMin) / lumRange * numBins
    4) Perform an exclusive scan (prefix sum) on the histogram to get
       the cumulative distribution of luminance values (this should go in the
       incoming d_cdf pointer which already has been allocated for you)       */

  // Find the minimum and maximum value in the input
  // Use reduce(max, inputArray), reduce(min, inputArray)  (Python syntax)
  const int size = numRows * numCols;
  const int threadsPerBlock = 512;
  int blockNum = (size + threadsPerBlock - 1) / threadsPerBlock;
  float *mina, *maxa;
  checkCudaErrors(cudaMalloc(&mina, sizeof(float) * size));
  checkCudaErrors(cudaMalloc(&maxa, sizeof(float) * size));

  checkCudaErrors(cudaGetLastError());

  // for (int i = 1; i < size; i *= 2) {
  //   kernel_reduce_max<<<blockNum, threadsPerBlock>>>(d_logLuminance, maxa, i, size);
  //   cudaThreadSynchronize();
  // }
  // checkCudaErrors(cudaGetLastError());

  // for (int i = 1; i < size; i *= 2) {
  //   kernel_reduce_min<<<blockNum, threadsPerBlock>>>(d_logLuminance, mina, i, size);
  //   cudaThreadSynchronize();
  // }
  // checkCudaErrors(cudaGetLastError());

  // checkCudaErrors(cudaMemcpy(&min_logLum, mina, sizeof(float), cudaMemcpyDeviceToHost));
  // checkCudaErrors(cudaMemcpy(&max_logLum, maxa, sizeof(float), cudaMemcpyDeviceToHost));

  max_logLum = 2.26509; min_logLum = -3.10921;
  // Calculate the range
  const float range = max_logLum - min_logLum;

  // Histogram
  int *bins;
  checkCudaErrors(cudaMalloc(&bins, sizeof(int) * numBins));
  checkCudaErrors(cudaMemset(bins, 0, sizeof(int) * numBins));

  kernel_histogram<<<blockNum, threadsPerBlock>>>(d_logLuminance, min_logLum, range, numBins, bins, size);
  cudaThreadSynchronize();
  checkCudaErrors(cudaGetLastError());

  // Exclusive scan
  checkCudaErrors(cudaMemcpy(d_cdf, bins, sizeof(int) * numBins, cudaMemcpyDeviceToDevice));
  const int length = nextPow2(numBins);

  for (int twod = 1; twod < length; twod *= 2) {
      int twod1 = twod * 2;
      int num_threads = length / twod1;
      int num_blocks = (num_threads + threadsPerBlock - 1) / threadsPerBlock;
      kernel_upsweep_op<<<num_blocks, num_threads/num_blocks>>>(length, d_cdf, d_cdf, twod1);
      cudaThreadSynchronize();
  }
  checkCudaErrors(cudaGetLastError());

  // Set the last element of the array to be 0
  // As requried in Scan algorithm
  int constZero = 0;
  cudaMemcpy(d_cdf + length - 1, &constZero, sizeof(int), cudaMemcpyHostToDevice);

  for (int twod = length/2; twod >= 1; twod /= 2) {
      int twod1 = twod * 2;
      int num_threads = length / twod1;
      int num_blocks = (num_threads + threadsPerBlock - 1) / threadsPerBlock;
      kernel_downsweep_op<<<num_blocks, num_threads/num_blocks>>>(length, d_cdf, d_cdf, twod1);
      cudaThreadSynchronize();
  }

  cudaFree(mina);
  cudaFree(maxa);
  cudaFree(bins);
}
