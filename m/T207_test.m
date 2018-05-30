## Copyright (C) 2018  -*- octave -*-

## This is initiated by Antonio
## see http://i56578-swl.blogspot.com/search/label/T-207 for details
##
## usage: test_T207(fn)
##
## input:  fn ... filename containing a string of 0's and 1's (no space between)
## output: m  ... T-207 success rate [0..1]

function m=T207_test(fn)
  fid    = fopen(fn);
  fsk    = fgetl(fid);
  fclose(fid);
  bits   = fsk=='1';
  k      = 14; # frame size

  for start=1:14
    test_bits       = bits(start:end);
    n               = floor(length(test_bits)/k);
    frames{start}   = reshape(test_bits(1:n*k),k,n)';

    sum_bits        = sum(frames{start}(:,1:12), 2);
    checksum        = frames{start}(:,13:14)*[2; 1]; ## 00 -> 0; 10 -> 2; 11->3, 01 -> 1

    test_sum        = mod(sum_bits, 4); ## 2,6,10 -> 2, 3,7,11 -> 3, 0,4,8,12 -> 0, 1,5,9 -> 1

    success(start)  = 0;
    success(start) += sum(checksum(test_sum==0) == 3);
    success(start) += sum(checksum(test_sum==1) == 1);
    success(start) += sum(checksum(test_sum==2) == 2);
    success(start) += sum(checksum(test_sum==3) == 0);
    printf("processed rows:%i success(%2d):%i (%5.1f%%) \n", n, start, success(start), 100*success(start)/n);
    success(start) /= n;
  end

  [m,i] = max(success);

  clf;

  subplot(2,3,[1 2 4 5]);
  n = size(frames{i})(1);
  imagesc([1.5 14.5], [1 n], frames{i});
  title(fn);
  xlabel 'index';
  ylabel 'frame number';
  line([13 13 15 15 13], [0 n n 0 0], 'color', 'red', 'linewidth', 3);

  subplot(2,3,3);
  stairs(100*success);
  axis([1 14 0 100]);
  title(fn);
  xlabel 'start bit offset';
  ylabel 'success (%)';
  grid on

  subplot(2,3,6);
  text(-0.2, 0.7, sprintf('T-207 success:  %.1f%%', 100*m), 'fontweight', 'bold', 'fontsize', 13);
  axis('off');

endfunction







