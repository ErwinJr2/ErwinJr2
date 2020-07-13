An auto-correlation calculator based on FFT
==========================================

This is a plain C implementation for calculating auto-correlation of real 
double inputs using FFT.

## Motivation

Often times people calculate FFT for convolution and auto-correlation rather
than if they really need the Fourier transform. There have already been a lot of
wheels for FFT but all the famous ones are not designed in a way that is simple
(yet still efficient) enough to integrate into codes just for auto-correlation
purposes: they usually require a large planner and spend a lot of time dealing
with complicated dividing-and-conquer strategies, which is not necessary for
auto-correlation purposes because we have the freedom for choosing the padding
length. 

## Usage

For a single shot of auto-correlation calculation, `autocorr` function is a one
stop solution. For multiple runs, this is also a good guideline for how to break
the routine down and to provide a plan and a memory pool to optimize the
calculation.

- `auto_plan` contains information about the dividing strategy for FFT. 
   It's thread safe, meaning a plan can be shared among different threads, 
   because everything in it is determined by the length of the array and is
   constant during FFT. 
- Besides an `auto_plan`, `autocorr_mem` also requires a `mempool`, whose length
  is determined by `mem_len(plan)`. This is the memory needed for FFT. The
  memory is **not** thread-safe and need to be created separately for different
  thread, but can be shared for different run of the same thread with out the
  need of initialization again.

## Compiler 

This code require gcc >= 5.4 or clang >= 7.0. Smaller version is not tested.

## Acknowledgment

The FFT part of the program was modified from
[pocketfft](https://gitlab.mpcdf.mpg.de/mtr/pocketfft).

This is copied from https://github.com/CareF/fftautocorr
