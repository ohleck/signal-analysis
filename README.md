# signal-analysis
C++/octave code for analyzing signals on VLF,LF,MF,HF

----

## Prerequisites
* make
* C++ compiler
* octave (+signal package)

## Contents
### FSK bit stream extraction
`m/analyze_fsk.m`

### LFSR detection
`m/lfsr_test.m`

### PLL
`m/pll_init.m` initialization of PLL state to be used together with `src/pll_cc.cc`
