## Copyright (C) 2018 Christoph Mayer -*- octave -*-

## usage: test_conv
##
## blind estimation of convolutional decoder properties
##
## this script finds the coefficients of a k=7 r=1/2 convolutional decoder
##
## output:
##  x ... random binary sequence
##  y ... concolutionally encoded x
##  z ... contains the result of the test of the property used to find a,b
##  s ... state space after all iterations

function [x,y,s,z]=test_conv()

  ## make up a random binary sequences
  M = 1000;
  x = rand(1,M) > 0.5;

  ## k=7 r=1/2 convolutional decoder (e.g. used in STANAG 4285)
  a = [1 0 1 1 0 1 1];
  b = [1 1 1 1 0 0 1];

  ## encode x
  y = encode(x, a,b);

  ## test the property used to find a,b
  z = test_property(y,a,b);
  __test_passed__ = all(z(1,:) == z(2,:))

  ## find a,b using a sieve
  s = find_ab(y);
endfunction

function s=find_ab(y)
  ## make up the search space
  s.a = make_state;
  s.b = make_state;

  ## extract even and odd bits
  u = y(1:2:end-1);
  v = y(2:2:end);

  counter = 0;
  ## for each iteration
  for j=1:30
    ## - append to the state the result of sum(a(i) u(i+j)), sum(b(i) v(i+j))
    for i=1:size(s.a,1)
      s.a(i,7+j) = mod(dot(s.a(i,1:7) , u(j+[0:6])), 2);
    end
    for i=1:size(s.b,1)
      s.b(i,7+j) = mod(dot(s.b(i,1:7) , v(j+[0:6])), 2);
    end

    ## - prune
    s.a = prune(s.a, s.b);
    s.b = prune(s.b, s.a);

    printf('iteration %2d:   number of states = [%3d, %3d]\n', ...
           j, size(s.a,1), size(s.b,1));

    if size(s.a,1)==0 || size(s.b,1)==0
      printf('\nno convolutional decoder detected\n');
      break;
    end

    counter += size(s.a,1) == 1 && size(s.b,1)==1;
    ## - stop if the state space was of size 1 for at least 5 iterations
    if counter > 5
      printf('\nconvolutional decoder detected:\n\ta = [%s]\n\tb = [%s]\n', ...
             num2str(s.b(1:7)), num2str(s.a(1:7)));
      break;
    end
  end
endfunction

function s1=prune(s1,s2)
  idx = zeros(1,size(s1,1));
  for i=1:size(s1,1)
    ## avoided an explicit 2nd loop by using any(all(...))
    idx(i) = any(all(s1(i,8:end) == s2(:,8:end), 2));
  end
  ## remove all rows in s1 for which no row in s2 was found
  s1(idx==0, :) = [];
endfunction

function z=test_property(y, a,b)
  u = y(1:2:end-1);
  v = y(2:2:end);
  z = [];
  for i=1:numel(u)-7
    z(1,i) = mod(sum(b .* u(i+[0:6])), 2);
    z(2,i) = mod(sum(a .* v(i+[0:6])), 2);
  end
endfunction

function s=make_state()
  for i=1:127
    s(i,:) = uint8(dec2bin(i,7)=='1');
  end
endfunction

function y=encode(x, a,b)
  y=[];
  for i=1:numel(x)-7
    y(end+1) = mod(sum(a .* x(i+[0:6])), 2);
    y(end+1) = mod(sum(b .* x(i+[0:6])), 2);
  end
endfunction
