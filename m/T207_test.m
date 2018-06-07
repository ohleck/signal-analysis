## Copyright (C) 2018  -*- octave -*-

## This is initiated by Antonio
## see http://i56578-swl.blogspot.com/search/label/T-207 for details
##
## usage: test_T207(fn)
##
## input:  fn ... filename containing a string of 0's and 1's (no space between)
## output: m  ... T-207 success rate [0..1]
##         a  ... aligned bit frames: crc in a(:,13:14)
##         p  ... permutation

function [m,a,p]=T207_test(fn)
  m      = 0;
  a      = [];
  fid    = fopen(fn);
  fsk    = fgetl(fid);
  fclose(fid);
  bits   = fsk=='1';
  k      = 14; # frame size

  modes = perms([0 1 2 3]);
  for start=1:14
    test_bits       = bits(start:end);
    n               = floor(length(test_bits)/k);
    frames{start}   = reshape(test_bits(1:n*k),k,n)';

    sum_bits        = sum(frames{start}(:,1:12), 2);
    checksum        = frames{start}(:,13:14)*[2; 1]; ## 00 -> 0; 10 -> 2; 11->3, 01 -> 1

    test_sum        = mod(sum_bits, 4); ## 2,6,10 -> 2, 3,7,11 -> 3, 0,4,8,12 -> 0, 1,5,9 -> 1

    for m=1:size(modes,1)
      success(start,m)  = 0;
      success(start,m) += sum(checksum(test_sum==0) == modes(m,1));
      success(start,m) += sum(checksum(test_sum==1) == modes(m,2));
      success(start,m) += sum(checksum(test_sum==2) == modes(m,3));
      success(start,m) += sum(checksum(test_sum==3) == modes(m,4));
      printf("mode:%i processed rows:%i success(%2d):%i (%5.1f%%) \n", m, n, start, success(start,m), 100*success(start,m)/n);
      success(start,m) /= n;
    end
  end

  [_m,i]   = max(success);
  [m,mode] = max(_m);
  p        = modes(mode,:);
  a        = frames{i(mode)};

  clf;

  subplot(2,3,[1 2 4 5]);
  n = size(a,1);
  imagesc([1.5 14.5], [1 n], a);
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
  text(-0.3, 0.79, sprintf('T-207(mode:%d) success: %.1f%% ', mode, 100*m), 'fontweight', 'bold', 'fontsize', 13);
  text(-0.3, 0.55, [sprintf('mode:%d [',mode) sprintf('%d', modes(mode,:)), ']'], 'fontweight', 'bold', 'fontsize', 13);
  axis('off');
endfunction
