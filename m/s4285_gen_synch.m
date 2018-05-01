## -*- octave -*-

function x=s4285_gen_synch
  state = [1 1 0 1 0];
  taps  = [0 0 1 0 1];
  x = [];
  for i=1:80
    x(i)  = [ state(end) ];
    state = [ mod(sum(state.*taps),2) state(1:end-1) ];
  endfor
endfunction

