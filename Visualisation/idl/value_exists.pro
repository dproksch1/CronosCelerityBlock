function value_exists, filename, varname

  data = rd_tfile(filename,3)
  sizeData = size(data)
  n = sizeData(2)
  num = 0

  for i =0l, n-1 do begin
     if(data(0,i) eq varname) then begin
        num += 1
     endif
  endfor
  return,num
end
