## Copyright (C) 2018 Christoph Mayer -*- octave -*-

## usage: [frac,a]=lfsr_test(bits, n)
##
## input:  sequence of bits and length of LFSR
## output: vector of taps 'a' and the fraction of bits compatile with the LFSR

function [frac,a]=lfsr_test(bits,n)
  if nargin != 2
    print_usage;
  endif

  m = zeros(n,n);
  for j=0:n-1
    m(1+j, 1:n) = bits(1+j:n+j);
  end
  u = bits(1+n:2*n)';

  a = gauss_mod2_cc(m, u);
  frac = 0;

  if !isempty(a)
    bits_pred = bits; ## predicted bits
    a = a';
    for i=n+1:length(bits)
      bits_pred(i) = mod(sum(bits_pred(i-n:i-1).*a), 2);
    end
    frac = sum(bits(2*n+1:end) == bits_pred(2*n+1:end))/(length(bits) - 2*n);
  end
endfunction

