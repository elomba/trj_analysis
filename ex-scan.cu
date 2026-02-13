//===============================================================================
// File: ex-scan.cu
//===============================================================================
// Purpose:
//   C++ wrapper functions for NVIDIA Thrust library parallel prefix sum
//   (exclusive scan) operations. Provides C-compatible interface for calling
//   Thrust functions from Fortran CUDA code via mod_thrust module.
//
// Key Functionality:
//   - Exclusive prefix sum (scan) for GPU device arrays
//   - Supports integer, single precision, and double precision data
//   - Wraps raw device pointers with thrust::device_ptr
//   - Provides extern "C" linkage for Fortran interoperability
//
// Functions:
//   scan_int_wrapper()    - Integer exclusive scan starting from 'start'
//   scan_float_wrapper()  - Single precision exclusive scan
//   scan_double_wrapper() - Double precision exclusive scan
//
// Algorithm:
//   Exclusive scan: output[i] = sum(input[0:i]) + start
//   Example: [1,0,2,2,1,3] with start=0 → [0,1,1,3,5,6]
//
// Usage:
//   Called from Fortran via mod_thrust interfaces
//   Used in cluster labeling algorithms (connected components)
//   Critical for efficient parallel prefix computation on GPU
//
// Compilation:
//   nvcc -c -arch=sm_XX ex-scan.cu
//   where XX is the GPU compute capability (e.g., sm_75, sm_80)
//
// Notes:
//   - Requires NVIDIA Thrust library (part of CUDA toolkit)
//   - Uses thrust::execution_policy for device execution
//   - Arrays must already be allocated on GPU device
//===============================================================================
// Filename: csort.cu
// nvcc -c -arch sm_13 csort.cu
// #include <thrust/device_vector.h>
// #include <thrust/sort.h>
#include <thrust/device_ptr.h>
#include <thrust/execution_policy.h>
#include <thrust/scan.h>
extern "C" {
// Sort for integer arrays
void scan_int_wrapper(int *data, int N, int start) {
  // Wrap raw pointer with a device_ptr
  thrust::device_ptr<int> dev_ptr(data);
  // Use device_ptr in Thrust sort algorithm
  //  int data[6] = {1, 0, 2, 2, 1, 3};
  thrust::exclusive_scan(thrust::device, dev_ptr, dev_ptr + N, dev_ptr, start);
}

// Sort for float arrays
void scan_float_wrapper(float *data, int N, float start) {
  thrust::device_ptr<float> dev_ptr(data);
  thrust::exclusive_scan(thrust::device, dev_ptr, dev_ptr + N, dev_ptr, start);
}

// Sort for double arrays
void scan_double_wrapper(double *data, int N, double start) {
  thrust::device_ptr<double> dev_ptr(data);
  thrust::exclusive_scan(thrust::device, dev_ptr, dev_ptr + N, dev_ptr, start);
}
}
