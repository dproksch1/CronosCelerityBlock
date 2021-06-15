PRO UNDEFINE, varname  
setstate = size(varname)
if(setstate(0) gt 1) then begin
    tempvar = SIZE(TEMPORARY(varname))
endif
END
