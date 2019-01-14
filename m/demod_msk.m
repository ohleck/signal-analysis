## Copyright (C) 2019  -*- octave -*-

## MSK demodulation
##   https://www.wescottdesign.com/articles/MSK/msk.pdf
##   MIL-STD-188-140.010193.pdf
##
## usage: [bX,bY,b]=demod_msk(fn, baud)
##
## input:  fn   ... filename containing IQ-mode samples; the signal is assumed to be centered around 0Hz
##         baud ... baud
## output: bX   ... demodulated bits cos(modulation carrier)
##         bY   ... demodulated bits sin(modulation carrier)
##         b    ... xor(bX,bY) used, e.g., for DGPS

function [bX,bY,b,z,zr,t]=demod_msk(fn, baud)
  [y,fs] = audioread(fn);
  z   = y(:,1)+1i*y(:,2);
  zr  = resample(z,1000,fs);
  fprintf(stdout, 'resampling done\n')
  fflush(stdout)
  zr  = zr(200:end);
  fsr = 1000;

  f  = @(z,df) z.*exp(-2i*pi*[1:numel(z)]'/fsr*df);

  s1 = pll_init(1/sqrt(2), 1,  baud/2, fsr);
  s2 = pll_init(1/sqrt(2), 1, -baud/2, fsr);

  theta1 = pll_cc(s1, zr.**2);
  theta2 = pll_cc(s2, zr.**2);

  zr .*= exp(-1i*(theta1+theta2)/4);

  t = (theta1-theta2)/4;
  fprintf(stdout, 'PLL done\n')
  fflush(stdout)

  figure 1;
  subplot(2,2,2);
  plot(sin(t), imag(zr), '.');
  title 'X'
  xlabel('sin(t)')
  ylabel('imag(z)');
  xlim([-1 1]); ylim(.8*[-1 1]);
  subplot(2,2,1);
  plot(cos(t), real(zr), '.');
  title 'Y'
  xlabel('cos(t)')
  ylabel('real(z)');
  xlim([-1 1]); ylim(.8*[-1 1]);

  bc  = [];
  zrc = real(zr) .* cos(t);
  b   = cos(t)>0;
  mu = 5*4/pi;
  while ~isempty(b)
    idx = find(b!=b(1), 1);
    if isempty(idx);
      bc(end+1) = mu*mean(zrc(1:end));
      b(1:end) = [];
    else
      bc(end+1) = mu*mean(zrc(1:idx(1)-1));
      b  (1:idx-1) = [];
      zrc(1:idx-1) = [];
    end
  end
  fprintf(stdout, 'bc done\n')
  fflush(stdout)

  bs  = [];
  zrs = imag(zr) .* sin(t);
  b   = sin(t)>0;
  while ~isempty(b)
    idx = find(b!=b(1), 1);
    if isempty(idx)
      bs(end+1) = mu*mean(zrs(1:end));
      b(1:end) = [];
    else
      bs(end+1) = mu*mean(zrs(1:idx(1)-1));
      b  (1:idx-1) = [];
      zrs(1:idx-1) = [];
    end
  end
  fprintf(stdout, 'bs done\n')
  fflush(stdout)

  subplot(2,2,3); hist(bc, [-2.2:.1:2.2]); grid on; title 'X'
  subplot(2,2,4); hist(bs, [-2.2:.1:2.2]); grid on; title 'Y'

  bX=bc>0;
  bY=bs>0;
  figure 2;
  b1 = bX(1:2:end);
  b2 = bX(2:2:end);
  b3 = bY(1:2:end);
  b4 = bY(2:2:end);
  f = @(x,k) reshape(x(1:floor(numel(x)/k)*k), k, floor(numel(x)/k))';
  k = 14;
  m = 32;
  a1 = f(b1, k);
  a2 = f(b2, k);
  a3 = f(b3, k);
  a4 = f(b4, k);

  subplot(2,2,1);
  frac=[]; _a={}; for i=1:k; for j=4:m; [frac(i,j),_a{i,j}] = lfsr_test (a1(:,i)',j); end; end; imagesc(frac', [0 1])
  title 'X(1:2:end)'
  xlabel 'start bit offset'
  ylabel 'LFSR length'
  title(colorbar, 'probability')

  subplot(2,2,2);
  frac=[]; _a={}; for i=1:k; for j=4:m; [frac(i,j),_a{i,j}] = lfsr_test (a2(:,i)',j); end; end; imagesc(frac', [0 1])
  title 'X(2:2:end)'
  xlabel 'start bit offset'
  ylabel 'LFSR length'
  title(colorbar, 'probability')

  subplot(2,2,3);
  frac=[]; _a={}; for i=1:k; for j=4:m; [frac(i,j),_a{i,j}] = lfsr_test (a3(:,i)',j); end; end; imagesc(frac', [0 1])
  title 'Y(1:2:end)'
  xlabel 'start bit offset'
  ylabel 'LFSR length'
  title(colorbar, 'probability')

  subplot(2,2,4);
  frac=[]; _a={}; for i=1:k; for j=4:m; [frac(i,j),_a{i,j}] = lfsr_test (a4(:,i)',j); end; end; imagesc(frac', [0 1])
  title 'Y(2:2:end)'
  xlabel 'start bit offset'
  ylabel 'LFSR length'
  title(colorbar, 'probability')

  bbc = reshape([bX;bX], 1, 2*numel(bX));
  bbs = reshape([bY;bY], 1, 2*numel(bY));
  bbc(1) = [];
  n = min([numel(bbc) numel(bbs)]);
  b = bitxor(bbc(1:n), bbs(1:n));
endfunction
