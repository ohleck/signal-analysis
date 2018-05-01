## -*- octave -*-

function y=s4285_test_deinterleave(M, data)
  N=32;

  de = {};
  for i=1:N
    de{i} = -ones(1,M*(N-i+0));
    length(de{i})
  end

  rows = 1+mod(9*[0:31], 32)

  y = [];
  for i=1:32:length(data)

    ## input
    for k=1:32
      de{k} = [data(i+k-1) de{k}];
    end

    ## output and shift
    for k=1:32
      l = rows(k);
      y(end+1) = de{l}(end);
      de{l}     = de{l}(1:end-1);
    end

  end
endfunction
