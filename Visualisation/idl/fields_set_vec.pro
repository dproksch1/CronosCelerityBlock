FUNCTION fields_set_vec, name_quant, num

@cronos_commons.h

avail = 0
case num of
   0: begin
      nameloc='x-Component'
      if(name_quant eq 'Velocity') then begin
         var = size(sx)
         if(var(0) gt 0) then avail = 1
      endif
      if(name_quant eq 'mag. field') then begin
         var = size(Bx)
         if(var(0) gt 0) then avail = 1
      endif
      if(name_quant eq 'curr. dens') then avail = 1
   end
   1: begin
      nameloc='y-Component'
      if(name_quant eq 'Velocity') then begin
         var = size(sy)
         if(var(0) gt 0) then avail = 1
       endif
      if(name_quant eq 'mag. field') then begin
         var = size(By)
         if(var(0) gt 0) then avail = 1
      endif
      if(name_quant eq 'curr. dens') then avail = 1
   end
   2: begin
      nameloc='z-Component'
      if(name_quant eq 'Velocity') then begin
         var = size(sz)
         if(var(0) gt 0) then avail = 1
      endif
      if(name_quant eq 'mag. field') then begin
         var = size(Bz)
         if(var(0) gt 0) then avail = 1
      endif
      if(name_quant eq 'curr. dens') then avail = 1
   end
   3: begin
      nameloc='Abs-Value'
      case name_quant of
         'Velocity': begin
            var1 = size(sx)
            var2 = size(sy)
            var3 = size(sz)
            if(var1(0) gt 0 and var2(0) gt 0 and var3(0) gt 0) then avail = 1
         end
         'mag field': begin
            var1 = size(sx)
            var2 = size(sy)
            var3 = size(sz)
            if(var1(0) gt 0 and var2(0) gt 0 and var3(0) gt 0) then avail = 1
         end
         'curr. dens': begin
            avail = 1
         end
         else: begin
            avail = 0
         end
      endcase
   end
   4: begin
      nameloc='Divergence'
      case name_quant of
         'Velocity': begin
            var1 = size(sx)
            var2 = size(sy)
            var3 = size(sz)
            if(var1(0) gt 0 and var2(0) gt 0 and var3(0) gt 0) then avail = 1
         end
         'mag field': begin
            var1 = size(sx)
            var2 = size(sy)
            var3 = size(sz)
            if(var1(0) gt 0 and var2(0) gt 0 and var3(0) gt 0) then avail = 1
         end
         'curr. dens': begin
            avail = 1
         end
         else: begin
            avail = 0
         end
      endcase
   end
   5: begin
      nameloc='Curl'
      case name_quant of
         'Velocity': begin
            var1 = size(sx)
            var2 = size(sy)
            var3 = size(sz)
            if(var1(0) gt 0 and var2(0) gt 0 and var3(0) gt 0) then avail = 1
         end
         'mag field': begin
            var1 = size(sx)
            var2 = size(sy)
            var3 = size(sz)
            if(var1(0) gt 0 and var2(0) gt 0 and var3(0) gt 0) then avail = 1
         end
         'curr. dens': begin
            avail = 1
         end
         else: begin
            avail = 0
         end
      endcase
   end
 endcase

data_set = {name:'', set:0}
data_set.name = nameloc         ;
data_set.set = avail            ;
return, data_set
END
