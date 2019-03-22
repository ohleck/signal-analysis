## Copyright (C) 2019 Christoph Mayer -*- octave -*-

## usage: conv_test
##
##  test a given bit sequence for convolutional encoding
##
## intput:
##  u ... 1st output
##  v ... 2nd output
##  a ... polyomial coefficients for 1st output
##  b ... polyomial coefficients for 2nd output
##
## output
##  success ... true if conv. encoding was detected, otheriwse false

function [success,z]=conv_test(u,v, a,b, pu, pv)
  narginchk(4,6);
  if nargin == 6
    assert(all(size(pu) == size(pv)))
    k = numel(pu);
  else
    pu = pv = 1;
    k = 1
  end
  pu
  pv
  k
  assert(all(size(a) == size(b)))
  m = numel(a);
  z = zeros(2,numel(u)-m, 'uint8');
  for i=1:numel(u)-m
    idx    = i+[0:m-1];
    z(1,i) = bitand(dot(b, u(idx).*pv(1+mod(idx-1, k))), 1);
    z(2,i) = bitand(dot(a, v(idx).*pu(1+mod(idx-1, k))), 1);
  end
  success = all(z(1,:) == z(2,:));
  z
endfunction
