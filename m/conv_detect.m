## -*- octave -*-

function [success,a,b]=conv_detect(u,v, m, pu,pv)
  narginchk(3,5);
  if nargin == 5
    assert(all(size(pu) == size(pv)))
    k = numel(pu);
  else
    pu = pv = 1;
    k = 1;
  end

  ## make up the search space
  s.a = make_state(m);
  s.b = make_state(m);

  success = false;
  a = []; b = [];
  counter = 0;
  ## for each iteration
  for j=1:30
    ## - append to the state the result of sum(a(i) u(i+j)), sum(b(i) v(i+j))
    idx = j+[0:m-1];
    for i=1:size(s.a,1)
      _b   = s.a(i,1:m);
      _b .*= pu(1+mod(idx-1, k));
      s.a(i,m+j) = bitand(dot(_b, u(idx)), 1);
    end
    for i=1:size(s.b,1)
      _a = s.b(i,1:m);
      _a .*= pv(1+mod(idx-1, k));
      s.b(i,m+j) = bitand(dot(_a, v(idx)), 1);
    end

    ## - prune
    s.a = prune(s.a, s.b, m);
    s.b = prune(s.b, s.a, m);

    printf('iteration %2d:   number of states = [%3d, %3d]\n', ...
           j, size(s.a,1), size(s.b,1));

    if size(s.a,1)==0 || size(s.b,1)==0
      printf('\nno convolutional decoder detected\n');
      break;
    end

    counter += size(s.a,1) == 1 && size(s.b,1)==1;
    ## - stop if the state space was of size 1 for at least 5 iterations
    if counter > 10
      success = conv_test(u,v,s.b(1:m), s.a(1:m));
      if success
        a = s.b(1:m);
        b = s.a(1:m);
        printf('\nconvolutional decoder detected:\n\ta = [%s]\n\tb = [%s]\n', ...
               num2str(a), num2str(b));
        else
          printf('\nno convolutional decoder detected\n');
      end
      break;
    end
  end
endfunction

function s1=prune(s1,s2,m)
  idx = zeros(1,size(s1,1));
  for i=1:size(s1,1)
    ## avoided an explicit 2nd loop by using any(all(...))
    idx(i) = any(all(s1(i,m+1:end) == s2(:,m+1:end), 2));
  end
  ## remove all rows in s1 for which no row in s2 was found
  s1(idx==0, :) = [];
endfunction


function s=make_state(m)
  for i=1:(2**m)-1
    s(i,:) = uint8(dec2bin(i,m)=='1');
  end
endfunction

