## Copyright (C) 2019  -*- octave -*-

## MSK demodulation
##   https://www.wescottdesign.com/articles/MSK/msk.pdf
##   http://everyspec.com/MIL-STD/MIL-STD-0100-0299/MIL-STD-188-140A_10192
##
## usage: [bX,bY,b]=demod_msk(fn, baud)
##
## input:  fn   ... filename containing IQ-mode samples; the signal is assumed to be centered around 0Hz
##         baud ... baud
## output: bX   ... demodulated bits cos(modulation carrier)
##         bY   ... demodulated bits sin(modulation carrier)
##         b    ... xor(bX,bY) used, e.g., for DGPS

function [bX,bY,b]=demod_msk(fn, baud)
  [y,fs] = audioread(fn);
  y   = single(y);
  z   = y(:,1)+1i*y(:,2);
#  z = z(2e6+[1:1000e3]);
  fsr = 5*baud
  zr  = single(resample(z,fsr,fs));
  fprintf(stdout, 'resampling done\n')
  fflush(stdout)

  s1 = pll_init(1/sqrt(2), 15,  baud/2, fsr);
  s2 = pll_init(1/sqrt(2), 15, -baud/2, fsr);

  theta1 = pll_cc(s1, zr.**2);
  theta2 = pll_cc(s2, zr.**2);

  zr .*= exp(-1i*(theta1+theta2)/4);
  t    = (theta1-theta2)/4;

  fprintf(stdout, 'PLL done\n')
  fflush(stdout)

  [bX,bY] = demod_msk_helper_cc(zr, t);

  idx=ceil(fsr/2) + [1:ceil(fsr/2)];

  figure 1;
  subplot(2,2,1);
  plot(cos(t(idx)), real(zr(idx)), '.');
  title 'X'
  xlabel('cos(t)')
  ylabel('real(z)');
  xlim([-1 1]); ylim(1.2*max(abs(real(z)))*[-1 1]);

  subplot(2,2,2);
  plot(sin(t(idx)), imag(zr(idx)), '.');
  title 'Y'
  xlabel('sin(t)')
  ylabel('imag(z)');
  xlim([-1 1]); ylim(1.2*max(abs(imag(z)))*[-1 1]);

  dX = 1.2*max(abs(bX));
  dY = 1.2*max(abs(bY));
  subplot(2,2,3); hist(bX, [-dX:.01:dX]); grid on; title 'X'; xlabel 'X soft decisions'
  subplot(2,2,4); hist(bY, [-dY:.01:dY]); grid on; title 'Y'; xlabel 'Y soft decisions'

  bX = bX > 0;
  bY = bY > 0;

  if false
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
  end
  bbc = reshape([bX;bX], 1, 2*numel(bX));
  bbs = reshape([bY;bY], 1, 2*numel(bY));
  bbc(1) = [];
  n = min([numel(bbc) numel(bbs)]);
  b = bitxor(bbc(1:n), bbs(1:n));
endfunction
