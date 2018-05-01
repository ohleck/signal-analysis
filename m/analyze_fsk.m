## Copyright (C) 2018 Christoph Mayer -*- octave -*-

## usage: [synch,shift,s,sf]=analyze_fsk(z, fs, baud)
##
## input: z - IQ samples, fs- sampling rate (Hz), baud
## output: synch - bit synchronization, sync.bits are the extracted bits, shift (Hz)

function [synch,shift,s,sf]=analyze_fsk(z,fs,baud)
  if nargin != 3
    print_usage;
  endif

  ## fm demodulation
  fm = @(z) arg(conj(z(1:end-1)).*z(2:end));
  s  = fm(z);

  ## estimation of mark and shift frequencies and recenter
  t      = estimate_tones(s);
  shift  = diff(t)/2/pi*fs; ## in Hz
  center = mean(t);
  z     .*= exp(-1i*center*reshape([1:length(z)], size(z)));

  ## normalize to +-1
  s = 2*(s-t(1))/diff(t)-1;

  ## filter and resample
  b  = fir1(100, 0.0225*4);
  sf = filter(b,1,s);

  thr = 20.0;
  s(s>+thr) = +thr;
  s(s<-thr) = -thr;

  fsr = fs/16;
  sr  = sf(1:16:end);

  ## bit synchronization and extraction
  synch.period  = fsr/baud;
  synch.alpha   = 0.5;
  synch.history = zeros(1,2*synch.period+1);
  synch.counter = 0;
  synch.t       = synch.period;
  synch.err     = 0;

  [synch, bits] = early_late_cc(synch, sr);
  synch.bits    = bits;
endfunction

function t=estimate_tones(s)
  t0 = mean(s);
  t  = [mean(s(s<t0)) mean(s(s>t0))];
  t0 = mean(t);
  for i=1:10
    t  = [mean(s(s<t0)) mean(s(s>t0))];
    t0 = mean(t);
  end
endfunction

