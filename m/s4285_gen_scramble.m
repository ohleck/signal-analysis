## -*- octave -*-

function x=s4285_gen_scramble
  state = [1 1 1 1 1 1 1 1 1];
  taps =  [0 0 0 0 1 0 0 0 1];
  x = [];
  for i=1:176
    x(i)  = sum(state(end-2:end).*[4 2 1]);
    for j=1:3
      state = [ mod(sum(state.*taps),2) state(1:end-1) ];
    end
  endfor
endfunction
