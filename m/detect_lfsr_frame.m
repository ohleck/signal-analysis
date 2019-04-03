## Copyright (C) 2018 Christoph Mayer -*- octave -*-

## usage: [frac,a,taps]=detect_lfsr_frame(fn, k, m)
##
## input:  fn   ... filename; file should contain only '0's and '1's
##         k    ... length of frame
##         m    ... max. LFSR length to test
## output: a    ... framed matrix of bits
##         frac ... success matrix [0,1]
##         taps ... cell matrix with taps

function [frac,a,taps]=detect_lfsr_frame(fn, k, m)
  if nargin != 3
    print_usage;
  endif

  fid = fopen(fn);
  s   = fgets(fid);
  fclose(fid);
  b   = s=='1';

  frac=[];
  taps={};
  n=floor(length(b)/k);
  a=reshape(b(1:n*k),k,n)';
  for i=1:k;
    for j=4:m;
      [frac(i,j),taps{i,j}]=lfsr_test(a(:,i)',j);
    end
  end

  [m,j]       = max(frac')
  [success,i] = max(m);
  j           = j(i);

  clf;

  subplot(2,3,[1 2 4 5]);
  imagesc(a);
  fn = strrep(fn, '_', '\_');
  title(fn);
  xlabel 'index';
  ylabel 'frame number';
  line(i+0.5*[-1 -1 1 1 -1], [0 n n 0 0], 'color', 'red', 'linewidth', 3);

  subplot(2,3,3);
  imagesc(frac', [0 1]);
  ylabel(colorbar('eastoutside'), 'success (%)', 'fontsize', 10);
  xlabel('start bit offset')
  ylabel('LFSR length')
  title(fn);
  set(get(gca, 'title'), 'fontsize', 10)

  subplot(2,3,6);
  text(-0.3, 1.1, sprintf('success: %.1f%% ', 100*success), 'fontweight', 'bold', 'fontsize', 13);

  taps{i,j}
  idx = j+1 - find(taps{i,j})(end:-1:1);
  s = 'S(i)=';
  for i=1:length(idx)
    s=[s sprintf('S(i-%d)', idx(i))];
    if mod(i,4)==0
      s=[s sprintf('\n      +')];
    else
      s=[s sprintf('+')];
    end
  end
  s
  text(-0.4, 0.62, s(1:end-1),  'fontsize', 11);
  axis('off');

end
