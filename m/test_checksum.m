## Copyright (C) 2018  -*- octave -*-

## usage: [counter,sum_data,checksum]=test_checksum(bits, k, m)
##
## input:  bits ... vector of bits \in {0,1}
##         k    ... frame length
##         m    ... length of checkum: a frame consists of k-m data bits and of m checksum bits
##
## output: counter  ... number of collisions for each frame start offset
##         sum_data ... data word for the offset with minimum number of collisions
##         checksum ... checksum for the offset with minimum number of collisions

function [counter,sum_data,checksum]=test_checksum(bits, k, m)

  for start=1:k
    test_bits       = bits(start:end);
    n               = floor(length(test_bits)/k);
    frames{start}   = reshape(test_bits(1:n*k),k,n)';

    sum_data{start}        = frames{start}(:,1:k-m) * 2.**([k-m-1:-1:0])';
    checksum{start}        = frames{start}(:,k-m+1:end)* 2.**([m-1:-1:0])';

    xx={};
    for i=1:2**m
      xx{i} = sum_data{start}(checksum{start}==i-1);
    end

    counter(start) = 0;
    for i=1:2**m
      for j=i+1:2**m
        counter(start) += nnz(xx{i}==xx{j}');
      end
    end
  end

  [m,i]    = min(counter);
  sum_data = sum_data{i};
  checksum = checksum{i};

  clf;
  stairs(counter);
  xlim([1 k]);
  xlabel 'frame offset'
  ylabel 'number of checksum colisions'
  title(sprintf('test for [data(1:%d) checksum(1:%d)] frames', k-m, m));
endfunction
