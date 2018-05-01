## -*- octave -*-

function b=viterbi_decode_by_hand(y)

  pol1 = [1,3,4,6,7];
  pol2 = [1,2,3,4,7];

  s = init;

  b = [];
  for i=1:2:length(y)-1
    s = prune(s, real(y(i))<0, real(y(i+1))<0, pol1, pol2);
    n = size(s,1);
    printf("STEP: n=%d\n", n);
    if n==1
      if isempty(b)
        b = s;
      else
        b(end+1) = s(end);
      endif
    end

    if n == 0
      s = init;
    end
    s = advance(s);
  end
endfunction

function s=prune(s,t1,t2, pol1, pol2)
  n = size(s,1);
  b = zeros(n,1);
  for i=1:n
    _t1 = mod(sum(s(i,pol1)), 2);
    _t2 = mod(sum(s(i,pol2)), 2);
    b(i) = ((_t1 == t1) && (_t2 == t2));
  end
  s = s(b==1,:);
endfunction

function s=advance(s)
  s_old = s(:, 1:6);
  n = size(s,1);
  s = [ zeros(n,1)==1 s_old; ones(n,1)==1 s_old];

endfunction

function s=init
  s=zeros(128,7);
  for i=1:128
    s(i,:) = dec2bin(i-1,7)=='1';
  end
endfunction
