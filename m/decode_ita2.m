## Copyright (C) 2018 Christoph Mayer -*- octave -*-

## usage: s=decode_ita2(ita2, msg)
##
## assumes one start and two stop bits
##
## input: ita code table, message bits
## output: decoded messages

function s=decode_ita2(ita2, msg)

  idx         = cat(1, ita2.bits)*2.**[0:4]';
  letter      = 1;
  state       = 0; # unlocked
  msg_counter = 1;
  s = {};
  i = 0;
  while i < length(msg)-9
    frame     = msg(i+[1:8]);
    ## check for one start and for two stop bits
    new_state = sum(frame([1 7 8]) == [0 1 1]) == 3;
    if new_state == 0
      i += 1;
      if new_state != state
        msg_counter += 1;
      end
    else
      k = find(sum(frame([2:6]).*2.**[0:4]) == idx);
      if ita2(k).num == 29
        letter = 1;
      elseif ita2(k).num == 30
        letter = 2;
      elseif letter == 1
        s{msg_counter}(end+1) = ita2(k).letter(1);
      else
        s{msg_counter}(end+1) = ita2(k).figure(1);
      end
      i += 8;
    end
    state = new_state;
  end
endfunction
