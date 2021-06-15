function value, filename, varname

  data = rd_tfile(filename,3)
  sizeData = size(data)
  n = sizeData(2)
  var = 0.

  for i =0l, n-1 do begin
     if(data(0,i) eq varname) then begin
        var = double(data(2,i))
     endif
  endfor
  return,var
end
