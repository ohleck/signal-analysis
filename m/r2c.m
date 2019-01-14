## Copyright (C) 2019  -*- octave -*-

## real -> complex conversion
##
## Usage: z=r2c(fn, fc)
##  a new wav files with complex IQ samples is created
#
## input:  fn   ... filename containing real samples
##         fc   ... center frequency
## output: z    ... complex samples shifted by -fc
##         y    ... real samples
##         fs   ... sampling frequency


function [z,y,fs]=r2c(fn, fc)
  [y,fs] = audioread(fn);
  f   = @(x,df) x.*exp(-2i*pi*[1:numel(x)]'/fs*df);
  b   = fir1(128, fc/fs*2);
  z   = filter(b, 1, f(y, -fc));
  fns = strrep(fn, '.wav', '_iq.wav');
  audiowrite(fns, [real(z) imag(z)], fs);

  subplot(2,1,1);
  plot(fs*(mod([0:numel(y)-1]/numel(y)+0.5,1)-0.5), abs(fft(y)))
  title(strrep(fn, '_', '\_'))
  xlabel 'frequency (Hz)'
  ylabel 'FFT(y)'

  subplot(2,1,2);
  plot(fs*(mod([0:numel(z)-1]/numel(z)+0.5,1)-0.5), abs(fft(z)))
  title(strrep(fns, '_', '\_'))
  xlabel 'frequency (Hz)'
  ylabel 'FFT(z)'
endfunction
