# signal-analysis
C++/octave code for analyzing signals on VLF,LF,MF,HF

----

## Prerequisites
* make
* C++ compiler
* octave (+signal package)

## Installation
* `make` builds the octave .oct files
* signal octave package: `pkg install -forge signal`, see also [http://wiki.octave.org/Octave-Forge](http://wiki.octave.org/Octave-Forge)

## Contents
### FSK bit stream extraction
`m/analyze_fsk.m`

### LFSR detection
`m/lfsr_test.m`

### PLL
`m/pll_init.m` initialization of PLL state to be used together
with `src/pll_cc.cc`

### Viterbi decoder (rate 1/2 constraint length k=N)
`include/viterbi2.hpp` using add-compare-select butterfly operations
without preprocessor macros but is within 15% of the KA9Q
implementation.

`include/viterbi2_simple.hpp` a more simple implementation which is
slower than the one above but may be more suited for SIMD

### STANAG 4285 decoder
Implementation along the [STANAG 4285
specification](http://www.n2ckh.com/MARS_ALE_FORUM/s4285.PDF)

`include/s4285_frame_detector.hpp` synch sequence detection and
doppler and carrier offset estimation

`include/s4285_channel_estimator.hpp` adaptive LMS filter with
decision feedback for channel estimation and equalization

`include/s4285_bitstream_decoder.hpp` deinterleaver and viterbi
decoder with soft symbol input

`src/s4285_cc.cc` octave interface

Currently the 600,1200,2400 baud modes with long interleaver are
supported. The 2400 baud mode is untested since I did not find any
signal on air using this mode
