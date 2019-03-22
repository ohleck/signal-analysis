## Copyright (C) 2018 Christoph Mayer -*- octave -*-

## usage: conv_encode
##
## input::
##  x ... random binary sequence
##  a ... polyomial coefficients for 1st output
##  b ... polyomial coefficients for 2nd output
##
## output:
##  y ... concolutionally encoded x

function y=conv_encode(x, a,b, pu, pv)
  narginchk(3,5);
  if nargin == 5
    assert(all(size(pu) == size(pv)))
    k = numel(pu)
  end
  assert(all(size(a) == size(b)))

  m = numel(a);
  y = zeros(1, 2*(numel(x)-m));
  for i=1:numel(x)-m
    idx = i+[0:m-1];
    y(2*i-1) = bitand(dot(a,x(idx)), 1);
    y(2*i)   = bitand(dot(b,x(idx)), 1);
  end
  ## apply puncturing
  if nargin == 5
    k = numel(pu);
    n = numel(y(1:2:end));
    y(1:2:end) .*= pu(1+mod([0:n-1], k));
    n = numel(y(2:2:end));
    y(2:2:end) .*= pv(1+mod([0:n-1], k));
  end
endfunction
