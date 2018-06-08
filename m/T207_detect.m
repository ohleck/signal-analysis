## Copyright (C) 2018  -*- octave -*-

## This is initiated by Antonio
## see http://i56578-swl.blogspot.com/search/label/T-207 for details
##
## usage: T207_detect(fn)
##
## input:  fn ... filename containing a string of 0's and 1's (no space between)
## output: m  ... T-207 success rate [0..1]
##         a  ... aligned bit frames: crc in a(:,13:14)
##         p  ... permutation

function [m,a,p]=T207_detect(fn)
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

    if max(mean(frames{start}(:,13:14))) > 0.9
      success(start,1:length(modes)) = -1;
      continue;
    end
    sum_bits        = sum(frames{start}(:,1:12), 2);
    checksum        = frames{start}(:,13:14)*[2; 1]; ## 00 -> 0; 10 -> 2; 11->3, 01 -> 1

    test_sum        = mod(sum_bits, 4); ## 2,6,10 -> 2, 3,7,11 -> 3, 0,4,8,12 -> 0, 1,5,9 -> 1

    for m=1:size(modes,1)
      success(start,m)  = 0;
      success(start,m) += sum(checksum(test_sum==0) == modes(m,1));
      success(start,m) += sum(checksum(test_sum==1) == modes(m,2));
      success(start,m) += sum(checksum(test_sum==2) == modes(m,3));
      success(start,m) += sum(checksum(test_sum==3) == modes(m,4));
      if success(start,m)/n > 0.9
        printf("mode:%i processed rows:%i success(%2d):%i (%5.1f%%) \n", m, n, start, success(start,m), 100*success(start,m)/n);
      end
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
  imagesc([0 13], [1:24], 100*success', [0 100]);
  ylabel(colorbar('eastoutside'), 'success (%)', 'fontsize', 10);
  yticks([1:1:24])
  yt={};
  for i=1:length(modes)
    yt{i} = sprintf('%d%d%d%d', modes(i,1), modes(i,2), modes(i,3), modes(i,4));
  end
  yticklabels(yt);
  set(gca, 'fontsize', 7);
  set(get(gca, 'xlabel'), 'fontsize', 10)
  set(get(gca, 'ylabel'), 'fontsize', 10)
  title(fn);
  set(get(gca, 'title'), 'fontsize', 10)
  xlabel 'start bit offset';
  ylabel 'permutation';

  subplot(2,3,6);
  text(-0.3, 0.79, sprintf('T-207[%s] success: %.1f%% ', yt{mode}, 100*m), 'fontweight', 'bold', 'fontsize', 13);
  axis('off');
endfunction
