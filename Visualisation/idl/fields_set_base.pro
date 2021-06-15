FUNCTION fields_set_base, num

@cronos_commons.h

avail = 0
vectorial = 0 ; State if vectorial field
all_comp = 0 ; State if all components are available
case num of
    0: begin
        nameloc='Density'
        var = size(rho)
        if(var(0) gt 0) then avail = 1
    end
    1: begin
        nameloc='Velocity'
        var1 = size(sx)
        var2 = size(sy)
        var3 = size(sz)
        if(var1(0) gt 0 or var2(0) gt 0 or var3(0) gt 0) then begin
           avail = 1
           vectorial = 1
        endif
        if(var1(0) gt 0 and var2(0) gt 0 and var3(0) gt 0) then all_comp = 1
    end
    2: begin
       nameloc='mag. field'
       var1 = size(Bx)
       var2 = size(By)
       var3 = size(Bz)
       if(var1(0) gt 0 or var2(0) gt 0 or var3(0) gt 0) then begin
          avail = 1
          vectorial = 1
       endif
       if(var1(0) gt 0 and var2(0) gt 0 and var3(0) gt 0) then all_comp = 1
    end
    3: begin
       nameloc='curr. dens'
       var1 = size(Bx)
       var2 = size(By)
       var3 = size(Bz)
       if(var1(0) gt 0 and var2(0) gt 0 and var3(0) gt 0) then begin
          avail = 1
          vectorial = 1
          all_comp = 1
       endif
    end
    4: begin
        nameloc='Etherm'
        var = size(eth)
        if(var(0) gt 0) then avail = 1
    end
    5: begin
        nameloc='Temp'
        var1 = size(rho)
        var2 = size(eth)
        if(var1(0) gt 0 and var2(0) gt 0) then avail = 1
    end
    6: begin
        nameloc='Ekin'
        var0 = size(rho)
        var1 = size(sx)
        var2 = size(sy)
        var3 = size(sz)
        if(var0(0) gt 0 and (var1(0) gt 0 or var2(0) gt 0 $
                             or var3(0) gt 0)) then avail = 1
    end
 endcase

data_set = {name:'', set:0, vectorial:0, all_comp:0}
data_set.name = nameloc         ;
data_set.set = avail            ;
data_set.vectorial = vectorial
data_set.all_comp = all_comp
return, data_set
END
