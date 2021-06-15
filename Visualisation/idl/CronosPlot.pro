;*************************************************************
;* PROGRAM: CronosPlot
;*
;* AUTHOURS: Jens Kleimann <jk@tp4.rub.de>
;*           Ralf Kissmann <ralf.kissmann@uibk.ac.at>
;* FROM:     5-July-2011
;*
;* PURPOSE:
;*   Widget-based app for idl to read and plot hdf5-Data produced
;*   with the CRONOS code
;*
;* LICENSE:
;*   Use as you see fit, we don't care (but would still
;*   apprechiate acknowledgements and bug reports).
;*
; *************************************************************
;
; V0.1 - First build on 09-June-2011
; V0.8 - First semi-stable version (has still issues regarding
;        special outputs like, e.g., curl
; *************************************************************

;--------------------------------------------------------------
;  Base configuration
;--------------------------------------------------------------
PRO settings                    ; default settings
@cronos_commons.h
  
  rev_text = ' V0.8'
  maindir='../../data/'

  h5_names_iterare = 2          ; Iterate through file or read OmNames 
  frames = [0,0];
  rim = 0
  proj_nseq = -1                ; rank of current project in list if p.s
  proj_name = ''
  filebase  = ''                ; filename for PS output
  step_frame = -23;
  time_frame = -23.
  tog_lilo  = 0                 ; linear scaling
  tog_iLine = 0                 ; no line integrals
  tog_rangex = 0                ; use complete range (x)
  tog_rangey = 0                ; use complete range (y)
  tog_drange = 0                 ; Use exact data range
  tog_target = 0                ; output target (default: PostScript)
  tog_cbar  = 0                 ; Choice of color bar (0=right, 1=below, 2=none)
  tog_xyrescale = 0             ; Different scalings for x & y (default: no)
  tog_geom = 0                  ; Grid geometry (default is Cartesian)
  str_void = "---------"        ; const string for disabled menu items
  numEntries = 0;               ; number of fields in file
  frame_loaded = -23;           ; Number of loaded frame
  iquantity = 0;                ; Dafault is density
  iquan_Load = 0;               ; Load all fields by default
  icoltab = 51                  ; Default colortable
  imodifier = 0                 ; Default: no modification to data
  ivec_mod = 0                  ; Default: pick x-component
  cutlayer = 0                  ; {rad|theta|phi} layer at which to cut
  icutplane = 0                 ; cut at constant rad (0), theta (1) or phi (2)
  n_gc =       0                ; # of ghost cells to cut
  win_size = [800, 800]         ; window size (pixels)
  scale_y_to_x = 1.             ; Ratio of dy to dx for plot (default: unity)
  data_range = [-10. , 10. ]    ; min/max values in dataset
  plot_range = data_range       ; min/max values to plot
  xrange = [ 0, 1 ]             ; Range of x-axis
  yrange = [ 0, 1 ]             ; Range of y-axis
  xcoord = [ 0, 0]
  ycoord = [ 0, 0]
  zcoord = [ 0, 0]
  w_is_open = 0
  
  ; All available fields (so far)

  tab_fields = [ str_void ]
  tab_quant = [str_void]
  tab_quant_all = ['Density' , 'Velocity' , 'mag. field' , 'curr. dens' ,$
                   'Etherm' , 'Temp' , 'Ekin'] ; Gives all possible quantities
  tab_vec_modifier = [str_void]
  tab_modifier = [ ' none ' , 'abs. value' ]
  quant_vec = [ 0 ]
  quant_all = [ 0 ]
  tab_cutplane = [ '(x,y) with k', '(x,z) with j', '(y,z) with i' ]

END

;===============================================================
;  Translate tab-input into iquantity
;===============================================================

PRO get_iquantity
@cronos_commons.h

  imin = iquantity
  For iqua_loc = imin, 6 Do Begin
     If(tab_quant[imin] eq tab_quant_all[iqua_loc]) Then Begin
        iquantity = iqua_loc
        Break;
     EndIf
  Endfor
  alert," iquantity: ",iquantity
  alert," for quant: " + tab_quant_all[iquantity]
END


;===============================================================
;  Function to build name array
;===============================================================

FUNCTION GetName, varname
  name = varname
  case varname of
     'rho': begin
        name = "Density n"
     end
     'v_x': begin
        name = "Velocity v_x"
     end
     'v_y': begin
        name = "Velocity v_y"
     end
     'v_z': begin
        name = "Velocity v_z"
     end
     'B_x': begin
        name = "mag. field Bx"
     end
     'B_y': begin
        name = "mag. field By"
     end
     'B_z': begin
        name = "mag. field Bz"
     end
     'Etherm': begin
        name = "Etherm E_th"
     end
  endcase
  return, name
END


;===============================================================
;  Scan maindir folder for available Cronos projects
;===============================================================

PRO proj_scan
  @cronos_commons.h

  Print, 'Reading path ' + maindir
  dirlen = StrLen(maindir)

  path=STRCOMPRESS(maindir+'*float/')
  dir_list = file_search(path,count=N_proj)

  proj_NFrames = INTARR(N_proj);
  proj_names = STRARR(N_proj);

  count_zerof = 0 ; num of projects without frames
  ipro = 0
  for idir=0,N_proj-1 do begin
     dir_curr = dir_list(idir);
     print,idir,dir_list(idir)
     pname = STRMID(dir_curr,dirlen,STRLEN(dir_curr)-dirlen-6)

     ; Determine number of frames for each project and only use
     ; project with at least one frame:

     frame_patt = dir_list(idir) + '/*' + pname + '*.h5'
     frame_list = File_Search(frame_patt, Count=N_frames_all)

     N_frames = N_frames_all
     
     for ifra=0, N_frames-1 do begin
        parallel = STRMATCH(frame_list[ifra], '*_coord*') ;
        if(parallel eq 1) then N_frames -= 1;
     endfor

     If(N_frames gt 0) Then Begin
        proj_names[ipro] = pname;
        proj_NFrames[ipro] = N_frames
        ipro += 1
     Endif Else Begin
;        print," in: ",ipro,N_proj
;        dir_list[ipro:N_proj-2] = dir_list[ipro+1:N_proj-1]
;        proj_names[ipro:N_proj-2] = proj_names[ipro+1:N_proj-1]
;        proj_NFrames[ipro:N_proj-2] = proj_NFrames[ipro+1:N_proj-1]
        count_zerof = count_zerof+1
;        ipro -= 1
;        N_proj -= 1
;        print," Exit: ",pname
;        break;
     Endelse

  Endfor ; loop over projects
  
  ; Use only non-empty projs
  print,N_proj,count_zerof
  N_proj -= count_zerof
  proj_names = proj_names[0:N_proj-1]
  proj_NFrames = proj_NFrames[0:N_proj-1]


  proj_list = STRARR(N_proj);
  ; Build list entries:
  for ipro=0, N_proj-1 Do Begin
     frame_nstr = StrCompress(proj_NFrames[ipro], /Remove_all)
     proj_list[ipro] = StrMid('00000' + frame_nstr, 4, /Reverse_offset) $
                       + '  ' + proj_names[ipro]
     print, ipro, '  ', proj_list[ipro]
  endfor

  if (proj_nseq ne -1) Then Begin ; must update frame slider
     proj_nfra = proj_NFrames[proj_nseq]
     frames[0] = Min( [frames[0], proj_nfra-1] )
     frames[1] = Min( [frames[1], proj_nfra-1] )
  Endif


END ; proj_scan

;===============================================================
;  Scan for frames of specific project
;===============================================================
PRO scan_frames
@cronos_commons.h

  pname = proj_names(proj_nseq)
  path = STRCOMPRESS(maindir+pname+'_float/')
  catname = STRCOMPRESS(maindir+pname+'.cat')

  allFrames=file_search(path,'*',count=N_frames_all)
;  frame_files=file_search(path,'*',count=N_frames_all)
  frame_files=StrArr(N_frames_all)

  count = 0
  for ifra=0, N_frames_all-1 do begin
     parallel = STRMATCH(allFrames[ifra], '*_coord*') ;
     if(parallel eq 0) then begin
        frame_files(count) = allFrames(ifra);
        count += 1
     endif
  endfor
  N_frames = count
  frame_files = frame_files[0:N_frames-1];
  
  timesteps = IntArr(N_frames)
  times = DblArr(N_frames)
  for ifra=0, N_frames-1 do begin
;     print,ifra,n_frames,count
     file = frame_files(ifra);

     ;Open the file:
     file_id = H5F_OPEN(file)
     
     ; Open the group
     GroupName="/Data/"
     group_id    = H5G_OPEN(file_id, GroupName)

     ; Get time entry:
     AttrName="time"
     attrset_id = H5A_OPEN_NAME(group_id, AttrName)
     times(ifra) = H5A_READ(attrset_id)

     ; Get timestep:
;     AttrName="timestep"
     AttrName="itime"
     attrset_id = H5A_OPEN_NAME(group_id, AttrName)
     timesteps(ifra) = H5A_READ(attrset_id)

     H5A_CLOSE,attrset_id
     H5G_CLOSE,group_id
     H5F_CLOSE,file_id
     
  endfor

  frames_list = StrArr(N_frames)

  for ifra=0, N_frames-1 do begin
     tstep_nstr = StrCompress(timesteps[ifra], /Remove_all)
     times_nstr = StrCompress(times[ifra], /Remove_all)
     frames_list[ifra] = StrMid('00000' + tstep_nstr, 4, /Reverse_offset) $
                         + '  ' + times_nstr
  endfor
  
  alert,"Clearing all arrays"
   ; Clear all arrays:
  var = size(rho)
  if(var(0) gt 0) then undefine, rho
  var = size(sx)
  if(var(0) gt 0) then undefine, sx
  var = size(sy)
  if(var(0) gt 0) then undefine, sy
  var = size(sz)
  if(var(0) gt 0) then undefine, sz
  var = size(Bx)
  if(var(0) gt 0) then undefine, Bx
  var = size(By)
  if(var(0) gt 0) then undefine, By
  var = size(Bz)
  if(var(0) gt 0) then undefine, Bz
  var = size(Ax)
  if(var(0) gt 0) then undefine, Ax
  var = size(Ay)
  if(var(0) gt 0) then undefine, Ay
  var = size(Az)
  if(var(0) gt 0) then undefine, Az
  var = size(eth)
  if(var(0) gt 0) then undefine, eth

End


;===============================================================
;  Scan for fields for specific frame:
;===============================================================

PRO scan_fields
@cronos_commons.h

  ; Loading the file
  filename = frame_files[chosen_frame]
  file_id = H5F_OPEN(filename);

  ; Open data-group
  GroupName="/Data/"

  group_id = H5G_OPEN(file_id, GroupName)

  ; Get number of entries:
  AttrName="Entries"
  attrset_id = H5A_OPEN_NAME(group_id, AttrName)
  dump = H5A_READ(attrset_id)
  numEntries = dump[0]
  H5A_CLOSE,attrset_id


  Case h5_names_iterare of

      1: Begin
                                ; Alternative to get dataset names
        numloc = 0              ;
 ;    dump = H5G_Get_Nmembers (group_id, GroupName)
        print," Dump: ",dump
                                ;    numEntries = dump[0]
        numEntriesLoc = numEntries ;
        For imember = 0, numEntries-1 Do Begin
           namesloc = H5G_Get_Member_Name (file_id, GroupName, imember)
           Print, imember, " : ", namesloc
           dummy="om_user"
           if(STRCMP(namesloc,dummy,7, /FOLD_CASE) EQ 1) then begin
              numEntriesLoc -= 1
           endif else begin
              field_names[numloc] = H5A_READ(attrset_id)
              numloc += 1
           endelse 
           print,num,H5A_READ(attrset_id)
           H5A_CLOSE,attrset_id
        endfor
        numEntries = numEntriesLoc
        
     End
     0: Begin
     
                                ; Get all dataset names:
        numloc = 0              ;
        numEntriesLoc = numEntries ;
        field_names = strarr(NumEntries)
        for num=0, NumEntries-1 do begin
           if(num lt 10) then begin
              AttrName=STRCOMPRESS("Name_om0"+STRING(num), /REMOVE_ALL)
           endif else begin
              AttrName=STRCOMPRESS("Name_om"+STRING(num), /REMOVE_ALL)
           endelse
           attrset_id = H5A_OPEN_NAME(group_id, AttrName)
           namesloc = H5A_READ(attrset_id)
           dummy="om_user"
           if(STRCMP(namesloc,dummy,7, /FOLD_CASE) EQ 1) then begin
              numEntriesLoc -= 1
           endif else begin
              field_names[numloc] = H5A_READ(attrset_id)
           numloc += 1
        endelse 
           print,num,H5A_READ(attrset_id)
           H5A_CLOSE,attrset_id
        endfor
        numEntries = numEntriesLoc
        
     End
     2: Begin

        QuantName = [ "rho",                      $
                      "v_x", "v_y", "v_z",        $
                      "B_x", "B_y", "B_z", "Etherm" ]

        ; Get all dataset names (for broken files)
        numloc = 0
        numEntriesLoc = numEntries ;
        field_names = strarr(NumEntries)
        num_v = 0;
        num_b = 0 ;

        for num=0, NumEntries-1 do begin
           if(num lt 10) then begin
              AttrName=STRCOMPRESS("Name_om0"+STRING(num), /REMOVE_ALL)
           endif else begin
              AttrName=STRCOMPRESS("Name_om"+STRING(num), /REMOVE_ALL)
           endelse
           attrset_id = H5A_OPEN_NAME(group_id, AttrName)
           namesloc = H5A_READ(attrset_id)
           dummy="o"
           if(STRCMP(namesloc,dummy,1, /FOLD_CASE) EQ 1) then begin
              numEntriesLoc -= 1
           endif else begin
              if(STRCMP(namesloc,"r", 1, /FOLD_CASE) EQ 1) then begin
                 field_names[numloc] = 'rho'
              Endif

              if(STRCMP(namesloc,"v", 1) EQ 1 and num_v eq 2) then begin
                 field_names[numloc] = 'v_z'
                 num_v += 1
              Endif
              if(STRCMP(namesloc,"v", 1) EQ 1 and num_v eq 1) then begin
                 field_names[numloc] = 'v_y'
                 num_v += 1
              Endif
              if(STRCMP(namesloc,"v", 1) EQ 1 and num_v eq 0) then begin
                 field_names[numloc] = 'v_x'
                 num_v += 1
              Endif

              if(STRCMP(namesloc,"B", 1) EQ 1 and num_b eq 2) then begin
                 field_names[numloc] = 'B_z'
                 num_b += 1
              Endif
              if(STRCMP(namesloc,"B", 1) EQ 1 and num_b eq 1) then begin
                 field_names[numloc] = 'B_y'
                 num_b += 1
              Endif
              if(STRCMP(namesloc,"B", 1) EQ 1 and num_b eq 0) then begin
                 field_names[numloc] = 'B_x'
                 num_b += 1
              Endif
              
              if(STRCMP(namesloc,"e", 1, /FOLD_CASE) EQ 1) then begin
                 field_names[numloc] = 'Etherm'
              Endif
              numloc += 1
           endelse 
           print,num,H5A_READ(attrset_id)
           H5A_CLOSE,attrset_id
        endfor
        numEntries = numEntriesLoc
        
     End
     
  EndCase
  
  print,"Arrays:"
  tab_fields = StrArr(NumEntries+1)
  tab_fields[0] = ' All '
  for num=0, NumEntries-1 do begin
     print,num,"  ",field_names[num]," ",GetName(field_names[num])
     tab_fields[num+1] = GetName(field_names[num])
  endfor



  ; Determine extent of arrays:
  DataName = STRCOMPRESS(field_names[0])
  dataset_id = H5D_OPEN(group_id, DataName)
  dataspace_id = H5D_GET_SPACE(dataset_id) ; get dataspace
  
  ndims = H5S_GET_SIMPLE_EXTENT_NDIMS(dataspace_id) ;
  nx = H5S_GET_SIMPLE_EXTENT_DIMS(dataspace_id)    ;
  If (ndims ne 3) Then Begin
     alert, "ERROR: Dim <> 3!"
     Stop
  Endif
  alert, "Nx = ", nx

  attrset_id = H5A_Open_Name (dataset_id, 'origin')
  xb         = H5A_Read      (attrset_id)
  alert, "x0 = ", xb

  attrset_id = H5A_Open_Name (dataset_id, 'delta')
  dx         = H5A_Read      (attrset_id)
  alert, "dx = ", dx

  if(attribute_avail(group_id, 'rim') eq 1) then begin
     attrset_id = H5A_Open_Name (group_id, 'rim')
;  rim        = Fix(H5A_Read      (attrset_id))
     dump        = H5A_Read      (attrset_id)
     rim = Fix(dump[0])
     alert, "rim = ", rim
  Endif Else Begin
     rim = 0
     alert, "Attribute rim not found -> setting rim to zero "
  EndElse

  H5A_Close, attrset_id
  H5D_Close, dataset_id

  ; Close file
  H5G_CLOSE,group_id
  H5F_CLOSE,file_id

 ; Clear all arrays:
  undefine, rho
  undefine, sx
  undefine, sy
  undefine, sz
  undefine, Bx
  undefine, By
  undefine, Bz
  undefine, Ax
  undefine, Ay
  undefine, Az
  undefine, eth


  catname = STRCOMPRESS(maindir + proj_name + '.cat')
  mx    = intarr(ndims)
  ; Get mx and compare to expection:
  if(value_exists(catname, "Nx") and value_exists(catname, "Ny") and $
     value_exists(catname, "Nz")) then begin
     mx[0] = FIX(value(catname,"mx")) - 1;
     mx[1] = FIX(value(catname,"my")) - 1;
     mx[2] = FIX(value(catname,"mz")) - 1;
  endif else begin
     mx[0] = FIX(value(catname,"mx"));
     mx[1] = FIX(value(catname,"my"));
     mx[2] = FIX(value(catname,"mz"));
  endelse

  if((mx[0] + 1 + 2*rim ne Nx[0]) or $
     (mx[1] + 1 + 2*rim ne Nx[1]) or $
     (mx[2] + 1 + 2*rim ne Nx[2])) then begin
     alert, "Error: mx does not match cat-File description "
;     break
  endif
  alert, "mx =",mx



End ; of scan_fields--------------------------------------------

;===============================================================
;    Read dataset from file
;===============================================================

FUNCTION ReadH5, file_id, DataName
@cronos_commons.h
  dxloc   = dblarr(3);
  xbloc   = dblarr(3);

  dataloc = fltarr(nx);

  dataset_id1 = H5D_OPEN(file_id, DataName)
  AttrName="origin"
  attrset_id1 = H5A_OPEN_NAME(dataset_id1, AttrName)
  xbloc = H5A_READ(attrset_id1)
  AttrName="delta"
  attrset_id1 = H5A_OPEN_NAME(dataset_id1, AttrName)
  dxloc = H5A_READ(attrset_id1)
  dataloc = H5D_READ(dataset_id1)
  H5D_CLOSE,dataset_id1
  datacont = {data: fltarr(nx), dx: dblarr(3), xb: dblarr(3)}
  datacont.data = dataloc;
  datacont.dx = dxloc
  datacont.xb = xbloc
  return, datacont
END

;===============================================================
;   Load all chosen fields
;===============================================================

PRO load_frame
  @cronos_commons.h

  filename = frame_files[chosen_frame]

  ; Determine number of fields to be loaded
  if(iquan_Load eq 0) then begin
     loadnum = numEntries;
  endif else begin
     loadnum = 1
  endelse

  ; Open file
  file_id = H5F_OPEN(filename)

  if(iquan_Load eq 0) then begin
     chname = field_names[0]
     loadnum = NumEntries
  Endif Else Begin
     chname = field_names[iquan_Load-1]
     loadnum = 1
  Endelse


  for dnum=0, loadnum-1 do begin
   
     if(loadnum gt 1) then begin
        chname = field_names[dnum]
     endif

     alert,' Loading: ' + chname
     case chname of
        'rho': begin
           rho = fltarr(nx)                      ;
           data = ReadH5(file_id, "/Data/rho") ;
           rho = data.data                           ;
           xb  = data.xb                             ;
           dx  = data.dx                             ;
        end
        'v_x': begin
           sx = fltarr(nx)                       ;
           data = ReadH5(file_id, "/Data/v_x") ;
           sx = data.data                            ;
           xb  = data.xb                             ;
           dx  = data.dx                             ;
      end
        'v_y': begin
           sy = fltarr(nx)                       ;
           data = ReadH5(file_id, "/Data/v_y") ;
           sy = data.data                            ;
           xb  = data.xb                             ;
           dx  = data.dx                             ;
        end
        'v_z': begin
           sz = fltarr(nx)                       ;
           data = ReadH5(file_id, "/Data/v_z") ;
           sz = data.data                            ;
           xb  = data.xb                             ;
           dx  = data.dx                             ;
        end
        'B_x': begin
           Bx = fltarr(nx)                       ;
           data = ReadH5(file_id, "/Data/B_x") ;
           Bx = data.data                            ;
           xb  = data.xb                             ;
           dx  = data.dx                             ;
        end
        'B_y': begin
           By = fltarr(nx)                       ;
           data = ReadH5(file_id, "/Data/B_y") ;
           By = data.data                            ;
           xb  = data.xb                             ;
           dx  = data.dx                             ;
      end
        'B_z': begin
           Bz = fltarr(nx)                       ;
           data = ReadH5(file_id, "/Data/B_z") ;
           Bz = data.data                            ;
           xb  = data.xb                             ;
           dx  = data.dx                             ;
        end
        'Etherm': begin
           eth = fltarr(nx)                         ;
           data = ReadH5(file_id, "/Data/Etherm") ;
;      data = ReadH5(size, file_id, "/Data/Temp");
           eth = data.data      ;
           xb  = data.xb        ;
           dx  = data.dx        ;
        end
     endcase

  endfor ; loop over fields to load
  H5F_CLOSE,file_id

  alert, "...done!"
  frame_loaded = chosen_frame

  ; Display quantities available for plotting:
  

  tab_quant = StrArr(7)
  quant_vec = IntArr(7)
  n_avail = 0
  for iqua=0,6 do begin
     setstate = fields_set_base(iqua);
     if(setstate.set eq 1) then begin
        tab_quant[n_avail] = setstate.name
;        quant_vec[n_avail] = setstate.vectorial
        n_avail += 1            ;
        quant_vec[iqua] = setstate.vectorial
     endif Else Begin
        quant_vec[iqua] = 0
     EndElse
  Endfor
  tab_quant = tab_quant[0:n_avail-1]
;  quant_vec = quant_vec[0:n_avail-1]
     
;     If (rim eq -1) Then Begin
;        alert, "No rim/ngc entry found."
;        ngc = Max([ ngc, ngc_ran[0] ]) ; clamp within ngc_ran
;        ngc = Min([ ngc, ngc_ran[1] ])
;        ngc_max = Floor( (Min(nx)-1) / 2. ) ; max. extent of ghost layer
;        If (ngc gt ngc_max) Then Begin      ; grid smaller than ghost layers
;           warning = "WARNING: grid too narrow => changing NGC" + $
;                     StrCompress(ngc) + " ->" + StrCompress(ngc_max)
;           alert, warning
;           ngc = ngc_max
;           ngc_ran[1] = ngc_max
;        Endif
;     Endif
     
  prep_data_2D
;  Endelse
  
END                             ; o

;===============================================================
;    Get extreme values for subregion of a field
;===============================================================

FUNCTION GetExtremeValues, field, xx, yy, xmin, xmax, ymin, ymax, mx
  
  imin=0;
  imax=0;
  for ii=0, mx[0] do begin
     if(xx(ii) lt xmin) then imin++;
     if(xx(ii) lt xmax) then imax++;
  endfor

  jmin=0;
  jmax=0;
  for jj=0, mx[1] do begin
     if(yy(jj) lt ymin) then jmin++
     if(yy(jj) lt ymax) then jmax++
  endfor

;  print," Error: ",imin,imax,jmin,jmax

  iext = imax-imin+1;
  jext = jmax-jmin+1;

  fieldLoc = dblarr(iext, jext)
  fieldLoc = field(imin:imax, jmin:jmax);

;  print," mehr: ",fieldLoc(2,2),field(imin+2,jmin+2),min(fieldLoc)

  locmin = min(fieldLoc);
  locmax = max(fieldLoc);
  dump = MOMENT(fieldLoc);
  locave = dump[0];

;  print," results: ",locmin, locmax, locave

  extdata = {min: 0.d, max: 0.d, ave: 0.d}
  extdata.min = locmin;
  extdata.max = locmax;
  extdata.ave = locave;
  
  return, extdata
END

FUNCTION GetExtremeValuesCurvilinear, field, xpolar, ypolar, xmin, xmax, ymin, ymax, mx
  locmin = max(field)
  locmax = min(field)
  locave = 0.;
  numcells = 0.;
  for ii=0, mx[0] do begin
     for jj=0, mx[1] do begin
        if(xpolar(ii,jj) ge xmin and xpolar(ii,jj) le xmax) then begin
           if (ypolar(ii,jj) ge ymin and ypolar(ii,jj) le ymax) then begin
              if(field(ii,jj) lt locmin) then locmin = field(ii,jj);
              if(field(ii,jj) gt locmax) then locmax = field(ii,jj);
              locave += field(ii,jj)
              numcells += 1
           endif
        endif
     endfor
  endfor
  locave /= numcells;
  statdata = {min: 0.d, max: 0.d, ave: 0.d}
  statdata.min = locmin;
  statdata.max = locmax;
  statdata.ave = locave;

  return,statdata
END

;===============================================================
;   Set available field modifiers
;===============================================================
PRO get_vec_modifiers
@cronos_commons.h

  tab_vec_modifier = StrArr(6)
  n_avail = 0
  for iqua=0,5 do begin
     setstate = fields_set_vec(tab_quant[iquantity], iqua);
     if(setstate.set eq 1) then begin
        tab_vec_modifier[n_avail] = setstate.name
        n_avail += 1;
     endif
  Endfor
  tab_vec_modifier = tab_vec_modifier[0:n_avail-1]

END ; get_modifiers


;===============================================================
;   Get ranges of plot
;===============================================================

PRO get_xrange
@cronos_commons.h
  Case icutplane of
     0: Begin
        xrange = [Min(xcoord), Max(xcoord)]
     End
     1: Begin
        xrange = [Min(xcoord), Max(xcoord)]
     End
     2: Begin
        xrange = [Min(ycoord), Max(ycoord)]
     End
  Endcase
END

PRO get_yrange
@cronos_commons.h
  Case icutplane of
     0: Begin
        yrange = [Min(ycoord), Max(ycoord)]
     End
     1: Begin
        yrange = [Min(zcoord), Max(zcoord)]
     End
     2: Begin
        yrange = [Min(zcoord), Max(zcoord)]
     End
  Endcase
END

;===============================================================
;   Other geometries
;===============================================================
PRO set_to_geom
@cronos_commons.h
  xrange_polar = [9d99, -9d99]
  yrange_polar = xrange_polar

  xplot_polar = xplot;
  yplot_polar = yplot;
  
  Case icutplane of
     0: Begin
        xpolar = dblarr(mx[0]+1, mx[1]+1);
        ypolar = dblarr(mx[0]+1, mx[1]+1) ;
     End
     1: Begin
        xpolar = dblarr(mx[0]+1, mx[2]+1);
        ypolar = dblarr(mx[0]+1, mx[2]+1) ;
     End
     2: Begin
        xpolar = dblarr(mx[1]+1, mx[2]+1) ;
        ypolar = dblarr(mx[1]+1, mx[2]+1) ;
     End
  EndCase

  Case tog_geom of

     1: Begin

        Case icutplane of
           0: Begin             ; r,phi -> x,y
              For i=0, mx[0] Do Begin
                 For j=0, mx[1] Do Begin
                    xval = xplot(i)*cos(yplot(j))
                    yval = xplot(i)*sin(yplot(j))
                    xpolar(i,j) = xval ;
                    ypolar(i,j) = yval ;
                 EndFor
              Endfor
           End
           1: Begin             ; r,z -> r,z
              For i=0, mx[0] Do Begin
                 For k=0, mx[2] Do Begin
                    xpolar(i,k) = xplot(i) ;
                    ypolar(i,k) = yplot(k) ;
                 Endfor
              EndFor
           End
           2: Begin             ; phi,z -> r*phi,z
              rad = xplot(cutlayer)
              For j=0, mx[1] Do Begin
                 For k=0, mx[2] Do Begin
                    xval = rad*xplot(j)
                    yval = yplot(k)
                    xpolar(j,k) = xval ;
                    ypolar(j,k) = yval ;
                 EndFor
              Endfor
           End
        Endcase

     End
     2: Begin
        
        Case icutplane of
           0: Begin   ; r,theta -> to r_cyl,z
              For i=0, mx[0] Do Begin
                 For j=0, mx[1] Do Begin
                    z = xplot(i)*cos(yplot(j))
                    r_cyl = xplot(i)*sin(yplot(j))
                    xpolar(i,j) = r_cyl ;
                    ypolar(i,j) = z ;
                 EndFor
              Endfor
           End
           1: Begin             ; r,phi -> x,y
              theta = yplot(cutlayer)
              For i=0, mx[0] Do Begin
                 For k=0, mx[2] Do Begin
                    xval = xplot(i)*cos(yplot(k))*sin(theta)
                    yval = xplot(i)*sin(yplot(k))*sin(theta)
                    xpolar(i,k) = xval ;
                    ypolar(i,k) = yval ;
                 Endfor
              EndFor
           End
           2: Begin             ; theta,phi -> r*theta,r*phi
              rad = xplot(cutlayer)
              For j=0, mx[1] Do Begin
                 For k=0, mx[2] Do Begin
                    xval = rad*xplot(j)
                    yval = rad*yplot(k)
                    xpolar(j,k) = xval ;
                    ypolar(j,k) = yval ;
                 EndFor
              Endfor
           End
        Endcase

     End

  EndCase

  xplot = xpolar
  yplot = ypolar
     

END


;===============================================================
;   Prepare field to be plotted (cut out subregion, cut range,...)
;===============================================================
PRO prep_plot
@cronos_commons.h
  
  ; Set xrange and stuff...
  Case icutplane of
     0: Begin
        xplot = xcoord
        yplot = ycoord
     End
     1: Begin
        xplot = xcoord
        yplot = zcoord
     End
     2: Begin
        xplot = ycoord
        yplot = zcoord
     End
  Endcase

  plot_data = selc_data;

  if(tog_geom ne 0) then set_to_geom

  if(tab_modifier[imodifier] eq 'abs. value') then begin
     plot_data = abs(plot_data)
  EndIf

  data_range = [Min(plot_data), Max(plot_data)]
  flag_islog = 0
  If (data_range[0] gt 0.) Then Begin ; all values to plot > 0
     flag_posdata = 1
     If(tog_lilo eq 1) Then Begin
        data_takeLog            ;  Taking logarithm of the data
     Endif
  Endif Else Begin
     flag_posdata = 0
  EndElse

  if(tog_rangex eq 0) then xrange = [Min(xplot), Max(xplot)]
  if(tog_rangey eq 0) then yrange = [Min(xplot), Max(xplot)]

  ; Set ranges if tog is exact
  If(tog_rangex eq 1 or tog_rangey eq 1) Then Begin
     if(tog_geom eq 0) then begin
        extvals = GetExtremeValues(plot_data, xplot, yplot, $
                                   xrange[0], xrange[1],    $
                                   yrange[0], yrange[1], mx)
     Endif Else Begin
        extvals = GetExtremeValuesCurvilinear(plot_data, xplot, yplot, $
                                              xrange[0], xrange[1],    $
                                              yrange[0], yrange[1], mx)
     EndElse
     data_range[0] = extvals.min        ;
     data_range[1] = extvals.max        ;
  Endif


  Case tog_drange of
     0: Begin                  ; exact
        plot_range = data_range
     End
     1: Begin                  ; rounded (use first 2 non-zero digits)
        distance = plot_range[1] - plot_range[0]
        digshift = 10.^(Floor(Alog10(distance))-1)
        plot_range[0] = digshift*Floor( data_range[0]/digshift )
        plot_range[1] = digshift* Ceil( data_range[1]/digshift )
     End
     2: Begin                  ; custom -> no action needed
     End
  EndCase

END

;===============================================================
;    Put field into intermediate array for plotting
;===============================================================
PRO prep_data_2D
@cronos_commons.h

  ; Set up all 3 dimensions:
  xcoord = xb[0] + dx[0]*(rim + Findgen(mx[0]+1))
  ycoord = xb[1] + dx[1]*(rim + Findgen(mx[1]+1))
  zcoord = xb[2] + dx[2]*(rim + Findgen(mx[2]+1))

  get_xrange;
  get_yrange;

  Case icutplane of
     0: selc_data = FltArr( mx[0], mx[1] )
     1: selc_data = FltArr( mx[0], mx[2] )
     2: selc_data = FltArr( mx[1], mx[2] )
  Endcase
  alert," Plot ranges: "
  alert," xrange: ",xrange
  alert," yrange: ",yrange

  if( tab_quant_all[iquantity] eq 'Temp')  then begin
     TempNorm = value(catname,"Norm_Temp");
     gamma = value(catname, "Adiabatic_exponent");
  endif

  Case icutplane of
     0: Begin

        Case tab_quant_all[iquantity] of
           'Density': selc_data = rho(*,*,cutlayer)
           'Velocity':  Begin
              ; Vectorial modificators

              Case tab_vec_modifier[ivec_mod] of
                 'x-Component': selc_data = sx(*,*,cutlayer)
                 'y-Component': selc_data = sy(*,*,cutlayer)
                 'z-Component': selc_data = sz(*,*,cutlayer)
                 'Abs-Value': selc_data = sqrt(sx(*,*,cutlayer)^2 + $
                                              sy(*,*,cutlayer)^2 + $
                                              sz(*,*,cutlayer)^2)
                 'Divergence': get_DivVelo
                 'Curl': get_CurlVelo
              EndCase

           End
           'mag. field': Begin
              ; Vectorial modificators

              Case tab_vec_modifier[ivec_mod] of
                 'x-Component': selc_data = Bx(*,*,cutlayer)
                 'y-Component': selc_data = By(*,*,cutlayer)
                 'z-Component': selc_data = Bz(*,*,cutlayer)
                 'Abs-Value': selc_data = sqrt(Bx(*,*,cutlayer)^2 + $
                                              By(*,*,cutlayer)^2 + $
                                              Bz(*,*,cutlayer)^2)
              EndCase

           End
           'curr. dens': Begin
              ; Vectorial modificators

              Case tab_vec_modifier[ivec_mod] of
                 'x-Component': selc_data = get_CurrDens(0)
                 'y-Component': selc_data = get_CurrDens(1)
                 'z-Component': selc_data = get_CurrDens(2)
                 'Abs-Value': Begin
                    CurrX = get_CurrDens(0)
                    CurrY = get_CurrDens(1)
                    CurrZ = get_CurrDens(2)
                    selc_data = sqrt(CurrX^2 + CurrY^2 + CurrZ^2)
                 End
              EndCase
           End
           'Etherm': selc_data = eth(*,*,cutlayer)
           'Temp': selc_data = TempNorm*(gamma-1.)*eth(*,*,cutlayer)
           'Ekin': selc_data = 0.5*rho(*,*,layer)*(sx(*,*,cutlayer)^2 + $
                                                  sy(*,*,cutlayer)^2 + $
                                                  sz(*,*,cutlayer)^2)
        Endcase

     End

     1: Begin
        Case tab_quant_all[iquantity] of
           'Density': selc_data = rho(*,cutlayer,*)
           'Velocity':  Begin
              ; Vectorial modificators

              Case tab_vec_modifier[ivec_mod] of
                 'x-Component': selc_data = sx(*,cutlayer,*)
                 'y-Component': selc_data = sy(*,cutlayer,*)
                 'z-Component': selc_data = sz(*,cutlayer,*)
                 'Abs-Value': selc_data = sqrt(sx(*,cutlayer,*)^2 + $
                                              sy(*,cutlayer,*)^2 + $
                                              sz(*,cutlayer,*)^2)
                 'Divergence': get_DivVelo
                 'Curl': get_CurlVelo
              EndCase

           End
           'mag. field': Begin
              ; Vectorial modificators

              Case tab_vec_modifier[ivec_mod] of
                 'x-Component': selc_data = Bx(*,cutlayer,*)
                 'y-Component': selc_data = By(*,cutlayer,*)
                 'z-Component': selc_data = Bz(*,cutlayer,*)
                 'Abs-Value': selc_data = sqrt(Bx(*,cutlayer,*)^2 + $
                                              By(*,cutlayer,*)^2 + $
                                              Bz(*,cutlayer,*)^2)
              EndCase

           End
           'curr. dens': Begin
              ; Vectorial modificators

              Case tab_vec_modifier[ivec_mod] of
                 'x-Component': selc_data = get_CurrDens(0)
                 'y-Component': selc_data = get_CurrDens(1)
                 'z-Component': selc_data = get_CurrDens(2)
                 'Abs-Value': Begin
                    CurrX = get_CurrDens(0)
                    CurrY = get_CurrDens(1)
                    CurrZ = get_CurrDens(2)
                    selc_data = sqrt(CurrX^2 + CurrY^2 + CurrZ^2)
                 End
              EndCase
           End
           'Etherm': selc_data = eth(*,cutlayer,*)
           'Temp': selc_data = TempNorm*(gamma-1.)*eth(*,cutlayer,*)
           'Ekin': selc_data = 0.5*rho(*,*,layer)*(sx(*,cutlayer,*)^2 + $
                                                  sy(*,cutlayer,*)^2 + $
                                                  sz(*,cutlayer,*)^2)
        Endcase
        selc_data = REFORM(selc_data)

     End
     2: Begin
        Case tab_quant_all[iquantity] of
           'Density': selc_data = rho(cutlayer,*,*)
           'Velocity':  Begin
              ; Vectorial modificators

              Case tab_vec_modifier[ivec_mod] of
                 'x-Component': selc_data = sx(cutlayer,*,*)
                 'y-Component': selc_data = sy(cutlayer,*,*)
                 'z-Component': selc_data = sz(cutlayer,*,*)
                 'Abs-Value': selc_data = sqrt(sx(cutlayer,*,*)^2 + $
                                              sy(cutlayer,*,*)^2 + $
                                              sz(cutlayer,*,*)^2)
                 'Divergence': get_DivVelo
                 'Curl': get_CurlVelo
              EndCase

           End
           'mag. field': Begin
              ; Vectorial modificators

              Case tab_vec_modifier[ivec_mod] of
                 'x-Component': selc_data = Bx(cutlayer,*,*)
                 'y-Component': selc_data = By(cutlayer,*,*)
                 'z-Component': selc_data = Bz(cutlayer,*,*)
                 'Abs-Value': selc_data = sqrt(Bx(cutlayer,*,*)^2 + $
                                              By(cutlayer,*,*)^2 + $
                                              Bz(cutlayer,*,*)^2)
              EndCase

           End
           'curr. dens': Begin
              ; Vectorial modificators

              Case tab_vec_modifier[ivec_mod] of
                 'x-Component': selc_data = get_CurrDens(0)
                 'y-Component': selc_data = get_CurrDens(1)
                 'z-Component': selc_data = get_CurrDens(2)
                 'Abs-Value': Begin
                    CurrX = get_CurrDens(0)
                    CurrY = get_CurrDens(1)
                    CurrZ = get_CurrDens(2)
                    selc_data = sqrt(CurrX^2 + CurrY^2 + CurrZ^2)
                 End
              EndCase
           End
           'Etherm': selc_data = eth(cutlayer,*,*)
           'Temp': selc_data = TempNorm*(gamma-1.)*eth(cutlayer,*,*)
           'Ekin': selc_data = 0.5*rho(*,*,layer)*(sx(cutlayer,*,*)^2 + $
                                                  sy(cutlayer,*,*)^2 + $
                                                  sz(cutlayer,*,*)^2)
        Endcase
        selc_data = REFORM(selc_data)

     End
        
  Endcase

 prep_plot;   Plot relevant settings like x-range etc.

;  data_range = [Min(selc_data), Max(plotdata)]
;  flag_islog = 0
;  If (data_range[0] gt 0.) Then Begin ; all values to plot > 0
;     flag_posdata = 1
;     If(tog_lilo eq 1) Then Begin
;        data_takeLog            ;  Taking logarithm of the data
;     Endif
;  Endif Else Begin
;     flag_posdata = 0
;  EndElse


END                        ; of PRO prep_data ------------------


;===============================================================
;  Take logarithm of data (if necessary)
;===============================================================

PRO data_takeLog
@cronos_commons.h

  If(flag_islog eq 0) Then Begin ; Only transform if not done yet
     If(flag_posdata eq 1) Then Begin
        plot_data   = ALog10(plot_data)
        data_range = ALog10(data_range)
        flag_islog = 1
        alert, "Using log data "
     Endif Else Begin
        alert, "WARNING: Cannot use log scaling because Min[plot] <= 0."
        tog_lilo = 0
        flag_islog = 0
     Endelse
  EndIf

END

;===============================================================
;  Function to reset x-range
;===============================================================

PRO reset_xrange
@cronos_commons.h

  If(tog_rangex eq 1) Then Begin ; Set user defined range
     
  EndIf

END

;===============================================================
;  Function to compute current density
;===============================================================

FUNCTION get_CurrDens, comp
@cronos_commons.h

  imin = IntArr(3);
  imax = IntArr(3);

  imin = [1, 1, 1]   ;
  imax = [mx[0]-1, mx[1]-1, mx[2]-1];

  Case icutplane of
     0: Begin
        imin[2] = cutlayer;
        imax[2] = cutlayer;
        CurrDens = dblarr(mx[0]+1,mx[1]+1)    ;
     End
     1: Begin
        imin[1] = cutlayer;
        imax[1] = cutlayer;
        CurrDens = dblarr(mx[0]+1,mx[2]+1)    ;
     End
     2: Begin
        imin[2] = cutlayer;
        imax[2] = cutlayer;
        CurrDens = dblarr(mx[1]+1,mx[2]+1)    ;
     End
  EndCase
  

  If (comp eq 0) Then Begin
     for i=imin[0], imax[0] do begin
        for j=imin[1], imax[1] do begin
           for k=imin[2], imax[2] do begin

              dyBz = 0.;
              If(mx[1] gt 1 and mx[2] gt 1) then begin
                 dyBz = (Bz(i,j+1,k  ) - Bz(i,j-1,k  ) +$
                         Bz(i,j+1,k-1) - Bz(i,j-1,k-1))/(4.*dx[1]) ;
              Endif
              dzBy = 0.;
              If(mx[1] gt 1 and mx[2] gt 1) Then Begin
                 dzBy = (By(i,j  ,k+1) - By(i,j  ,k+1) +$
                         By(i,j-1,k-1) - By(i,j-1,k-1))/(4.*dx[2]) ;
              Endif

              Result = dyBz - dzBy;

              Case icutplane of
                 0: CurrDens(i,j) = Result;
                 1: CurrDens(i,k) = Result;
                 2: CurrDens(j,k) = Result ;
              EndCase
           
           Endfor
        Endfor
     Endfor
  Endif
  return,CurrDens

END

;===============================================================
;  Function to compute divergence of velocity
;===============================================================

PRO get_DivVelo
@cronos_commons.h
  
  imin = IntArr(3);
  imax = IntArr(3);

  imin = [1, 1, 1]   ;
  imax = [mx[0]-1, mx[1]-1, mx[2]-1];

  Case icutplane of
     0: Begin
        imin[2] = cutlayer;
        imax[2] = cutlayer;
     End
     1: Begin
        imin[1] = cutlayer;
        imax[1] = cutlayer;
     End
     2: Begin
        imin[2] = cutlayer;
        imax[2] = cutlayer;
     End
  EndCase

  for i=imin[0], imax[0] do begin
     for j=imin[1], imax[1] do begin
        for k=imin[2], imax[2] do begin

           dxVx = 0.;
           if(mx[0] gt 1) then begin
              dxVx = (sx(i+1,j,k) - sx(i-1,j,k))/(2.*dx[0]) ;
           Endif

           dyVy = 0.;
           if(mx[1] gt 1) then begin
              dyVy = (sy(i,j+1,k) - sy(i,j-1,k))/(2.*dx[1]) ;
           Endif

           dzVz = 0.;           
           if(mx[2] gt 1) then begin
              dzVz = (sz(i,j,k+1) - sz(i,j,k-1))/(2.*dx[2]) ;
           EndIf

           Case icutplane of
              0: selc_data(i,j) = dxVx + dyVy + dzVz ;
              1: selc_data(i,k) = dxVx + dyVy + dzVz ;
              2: selc_data(j,k) = dxVx + dyVy + dzVz ;
           Endcase

        Endfor
     EndFor
  EndFor


END


;===============================================================
;  Function to compute curl of velocity
;===============================================================

PRO get_CurlVelo
@cronos_commons.h
  
  imin = IntArr(3);
  imax = IntArr(3);

  imin = [1, 1, 1]   ;
  imax = [mx[0]-1, mx[1]-1, mx[2]-1];

  Case icutplane of
     0: Begin
        imin[2] = cutlayer;
        imax[2] = cutlayer;
     End
     1: Begin
        imin[1] = cutlayer;
        imax[1] = cutlayer;
     End
     2: Begin
        imin[2] = cutlayer;
        imax[2] = cutlayer;
     End
  EndCase

  for i=imin[0], imax[0] do begin
     for j=imin[1], imax[1] do begin
        for k=imin[2], imax[2] do begin

           dxVy = 0.;
           dxVz = 0.;
           if(mx[0] gt 1) then begin
              dxVy = (sy(i+1,j,k) - sy(i-1,j,k))/(2.*dx[0]) ;
              dxVz = (sz(i+1,j,k) - sz(i-1,j,k))/(2.*dx[0]) ;
           Endif

           dyVx = 0.;
           dyVz = 0.;
           if(mx[1] gt 1) then begin
              dyVx = (sx(i,j+1,k) - sx(i,j-1,k))/(2.*dx[1]) ;
              dyVz = (sz(i,j+1,k) - sz(i,j-1,k))/(2.*dx[1]) ;
           Endif

           dzVx = 0.;           
           dzVy = 0.;           
           if(mx[2] gt 1) then begin
              dzVx = (sx(i,j,k+1) - sx(i,j,k-1))/(2.*dx[2]) ;
              dzVy = (sy(i,j,k+1) - sy(i,j,k-1))/(2.*dx[2]) ;
           EndIf

           curl = sqrt((dyVz - dzVy)^2 + $
                       (dzVx - dxVz)^2 + $
                       (dxVy - dyVx)^2);
           Case icutplane of
              0: selc_data(i,j) = curl ;
              1: selc_data(i,k) = curl ;
              2: selc_data(j,k) = curl ;
           Endcase

        Endfor
     EndFor
  EndFor


END

;===============================================================
;   Plotting routine
;===============================================================
PRO do_plot
@cronos_commons.h

; --- plot layout constants -----------------------------------------
  xBo   = 0.1                  ; relative border thickness (x dir.)
  wcb   = xBo*0.4               ; relative width of color bar (x d.)
  wan   = 0.35                  ; relative width of annotation box
  han   = 0.15                  ; relative height of annotation box
  nlevels = 24
  nColors = 254.               
  thickn  = [1  , 3]            ; field line thickness for {Xwin|PS}
  fillco  = [0  , 255]          ;    circle fill color for    -"-
  charsz  = [1.2, 2.5]          ;  annotation charsize for    -"-
  ndig    = 5                   ; # of frame digits in PS filename
  coltab  = 39                  ; color table (39 = rainbow + white)
  lev0    = 0                   ; level of zero line
; -------------------------------------------------------------------

  title = tab_quant[iquantity]

  If (tog_lilo eq 1) Then title = 'log( ' + title + ' )'


  ;size ratio of contour plot
  size_ratio = (xrange[1]-xrange[0]) / (yrange[1]-yrange[0]) 
  size_ratio /= scale_y_to_x
  win_ratio = (1.* win_size[0]) / (1.* win_size[1]) ; window aspect ratio
  yBo = xBo * 0.8                              ; border thickness (y dir.)
  hcb = wcb * win_ratio                        ; color bar height (y d.)

  If (tog_cbar ne 1) Then Begin                ; check available space for...
     avail_sq = [1-4*xBo-wcb-wan, 1-2*yBo    ] ;  (sq)uare layout
     avail_lo = [1-2.4*xBo-wcb  , 1-3*yBo-han] ;  (lo)ng   layoout
  Endif Else Begin                             ; same for {color bar below}
     avail_sq = [1-3*xBo-wan, 1-3*yBo-hcb    ]
     avail_lo = [1-2*xBo    , 1-3*yBo-hcb-han]
  Endelse

  ; relative plot positions:
  pos_tx = [xBo, 1-yBo-han, xBo+wan, 1-yBo] ; text field
  zoom_sq = Min(avail_sq * [win_ratio/size_ratio, 1] )
  zoom_lo = Min(avail_lo*win_size/win_size[1] / [size_ratio, 1] )
  If (zoom_sq ge zoom_lo) Then Begin ; -> choose square layout
     wcp = zoom_sq / (win_ratio/size_ratio)     ; rel. width & height of
     hcp = zoom_sq                   ;  contour plot w/o color bar
     print,"square layout"
     If (tog_cbar ne 1) Then Begin ; Square + Right
        pos_lp = [    xBo        , yBo, 1-3*xBo-wcb-wcp, 1-2*yBo-han]
        pos_mp = [1-2*xBo-wcb-wcp, yBo, 1-2*xBo-wcb    ,     yBo+hcp]
        pos_cb = [1-  xBo-wcb    , yBo, 1-  xBo        ,     yBo+hcp]
     Endif Else Begin           ;    Square + Below
        pos_lp = [  xBo    ,   yBo    , 1-2*xBo-wcp, 1-2*yBo-han    ]
        pos_mp = [1-xBo-wcp, 2*yBo+hcb, 1-  xBo    ,   2*yBo+hcb+hcp]
        pos_cb = [1-xBo-wcp,   yBo    , 1-  xBo    ,     yBo+hcb    ]
     Endelse

  Endif Else Begin              ; -> choose long layout
     print,"long layout"
     wcp = zoom_lo / (win_ratio/size_ratio)
     hcp = zoom_lo
     If (tog_cbar ne 1) Then Begin ; Long + Right
        print,"right"
        pos_lp = [2*xBo+wan, 2*yBo+hcp,  1-xBo       , 1-yBo    ]
        pos_mp = [  xBo    ,   yBo    ,    xBo+wcp   ,   yBo+hcp]
        pos_cb = [1.4*xBo+wcp,   yBo    , 1.4*xBo+wcp+wcb,   yBo+hcp]
     Endif Else Begin           ; Long + Below
        print,"below"
        pos_lp = [2*xBo+wan, 3*yBo+hcb+hcp, 1-xBo, 1-yBo        ]
        pos_mp = [1-xBo-wcp, 2*yBo+hcb    , 1-xBo, 2*yBo+hcb+hcp]
        pos_cb = [  xBo    ,   yBo        , 1-xBo,   yBo+hcb    ]
     Endelse
  Endelse

  nframes = 1  

  ; Prepare and load colortables.
  ommax = plot_range[1]
  ommin = plot_range[0]

  loadct,0
   @colortabs

   !p.background=16777215L
   !p.color=0L


                                ; Prepare Plot window:
   If (tog_target eq 1) Then Begin ; output to Xwindow
      Set_Plot, 'X'
      Window, 1, Retain=2, Xsize=win_size[0], Ysize=win_size[1]
      Device, Decompose=0
      w_is_open = (0 eq 0)      ; window_open <- true
      alert, "Plotting frame #" + StrCompress(frames[0], /Remove_all)
      
      !P.Font=0
                                ; Here we use specfic -readable-
                                ; xwindow font. Available fonts can be
                                ; found from foo = xfont()
      DEVICE, SET_FONT='-adobe-courier-bold-r-normal--14-100-100-100-m-90-iso8859-1'
      !P.THICK=2.0
      !X.THICK=2.0
      !Y.THICK=2.0
      !P.CHARTHICK=2.0
      !P.CHARSIZE=2.
      erase,255
      
   Endif Else Begin             ; output to PS file

      progress_message = "Plotting frame "
      alert, progress_message
      
      Set_Plot, 'PS'

      filename = filebase + ".ps"
      alert, "Output file name: " + filename

      La4 = [29.7 , 21.0]                ; DIN A4 format in cm
      dpc = Max( win_size / La4 )        ; px/cm (pad where needed)
      ps_size = win_size / dpc           ; size of PS plot in cm
      ps_offs = (La4 - ps_size) / 2.0    ; padding in cm
      Device, Filename=filename, /Color, Bits=24, /Landscape, /helvetica,$
              Xsize=ps_size[0], Xoffset=       ps_offs[1], $
              Ysize=ps_size[1], Yoffset=La4[0]-ps_offs[0]

      !P.Font=0
      !P.THICK=4.0
      !X.THICK=6.0
      !Y.THICK=6.0
      !P.CHARTHICK=5.0
      !P.CHARSIZE=1.2
   Endelse
   !P.COLOR=0
   !P.BACKGROUND=255


   If (data_range[1] eq data_range[0]) Then Begin
      alert, "WARNING: field is constant."
   Endif

   rangediff =     plot_range[1]  -     plot_range[0]
   rangeextr = Abs(plot_range[1]) + Abs(plot_range[0])
   If (rangediff le (1e-6)*rangeextr) Then Begin
      alert, "ERROR: Cannot plot min=max range."
   Endif Else Begin

; BEGIN { annotation box }
      If (tog_target eq 1 or tog_target eq 0) Then Begin
         coordname = [ "z", "y", "x" ]
         coordval  = xb[icutplane] + (cutlayer) * dx[icutplane]
         If (Abs(coordval/dx[icutplane]) lt 1e-6) Then coordval = 0.
         output = [ "Project:  " + proj_name, $
                    "Grid:    [" +            $
                    StrCompress(nx[0]-2*rim) + " x" + $
                    StrCompress(nx[1]-2*rim) + " x" + $
                    StrCompress(nx[2]-2*rim) + " ] ", $
;                    cat_string[0] + cat_vstr[0], $
;                    cat_string[1] + cat_vstr[1], $
;                    cat_string[2] + cat_vstr[2], $
                    "Quantity: "        + title, $
                    "Time:   "          + Strcompress(time_frame) + $
                    " at "+coordname[icutplane] + " =" + $
                    StrCompress(coordval) ]
         For row = 0, N_Elements(output)-1 Do Begin
            Xyouts, pos_tx[0], -0.02+pos_tx[3]-0.04*row, $
                    output[row], Charsize=charsz[tog_target], /Normal
         Endfor
      EndIf
; END { annotation box }

; BEGIN { colorbar } DELETE??
      ; Annotation format for colorbar
      if(ommax gt 100.) then begin
        Digits = 1+FIX(ALOG10(ommax))
        formatString = STRCOMPRESS('(I'+STRING(Digits)+')',/REMOVE_ALL)
     endif else if (ommax lt 0.01 and ommax ge 0.0001) then begin
        Digits = 2+ABS(FIX(ALOG10(ommax)))
        DigVar = Digits+3
        formatString = STRCOMPRESS('(F'+STRING(DigVar)+'.'+ $
                                   STRING(DIGITS)+')',/REMOVE_ALL)
     endif else if (ommax lt 0.0001 and ommin gt -0.0001) then begin
        formatString = STRCOMPRESS('(e8.1)') ;
     endif else begin
        formatString = STRCOMPRESS('(F6.3)') ;
     endelse

     if (tog_cbar ne 2) Then Begin
        If (tog_cbar eq 0) Then Begin

           COLORBAR, NCOLORS=nColors,BOTTOM=bottom    $
                     ,POSITION=pos_cb  $
                     ,/vertical,/right,MAXRANGE=userLevels[nlevels-1]  $
                     ,MINRANGE=userLevels[0],FORMAT=formatstring
           
        Endif Else Begin
           
           COLORBAR, NCOLORS=nColors,BOTTOM=bottom    $
                     ,POSITION=pos_cb  $
                     ,/right,MAXRANGE=userLevels[nlevels-1]  $
                     ,MINRANGE=userLevels[0],FORMAT=formatstring
           
        Endelse
     EndIf

; END { colorbar }
        
; BEGIN { contour plot }

     xtitle = ['x' , 'x', 'y']
     ytitle = ['y' , 'z', 'z']

     Contour, plot_data, xplot, yplot, Position=pos_mp, $
              background=255,color=0, $
              Levels=userLevels, C_colors=userColors, /Fill, $
              Xrange=xrange, Xstyle=1, $
              Yrange=yrange, Ystyle=1, $
              Zrange=[xrange[0],yrange[1]],$
              xtitle=xtitle[icutplane],ytitle=ytitle[icutplane],$
              /Noerase

; END { contour plot }
        
     If (tog_target eq 0) Then Begin
        Device, /Close
        Set_Plot, 'X'
     Endif
     
  Endelse
  alert, "...finished."

END

;===============================================================
;  Respond to widget action
;===============================================================
PRO ptest_event, event          ; respond to user action
@cronos_commons.h

  Widget_Control, event.id, Get_Uvalue=uvalue
  Print, "ev.id = ",event.id    ; to see if widget is still alive

  Case uvalue of
     'selc_proj' : Begin
        proj_nseq = event.index
 ;       menu_item = proj_list[proj_nseq]
        proj_nfra = proj_nframes[proj_nseq]
        proj_name = proj_names[proj_nseq]
;        If (proj_nfra gt 1) Then Begin
;           frames[0] = Min([ frames[0], proj_nfra-1 ])
;           frames[1] = Min([ frames[1], proj_nfra-1 ])
;        Endif Else Begin
;           tog_fram = 2
;           frames = [0, 0]
;        Endelse
     End
     'butt_scan' : Begin
        alert, "Scanning folder, please wait..."
        proj_scan
        alert, "done."
     End
     'butt_scan_frames' : Begin
        alert, "Scanning frames, please wait..."
        scan_frames
        alert, "done."
     End
     'selc_frame' : Begin
        chosen_frame = event.index
        step_frame = timesteps[chosen_frame]
        time_frame = times[chosen_frame]
        alert, "Scanning fields, please wait..."
        scan_fields
        alert, "done."
     End
;     'butt_scan_fields' : Begin
;        alert, "Scanning fields, please wait..."
;        scan_fields
;        alert, "done."
;     End
     'selc_quanLoad' : Begin
        cur_text = Widget_Info (id_selc_quanLoad, /Combobox_Gettext)
;        If (cur_text ne str_void) Then Begin
        iquan_Load = event.index
;           prep_data
;        Endif
     End
     'butt_load' : Begin
        tog_geom = 0            ; Reset to Cartesian grid
        tog_rangex = 0          ; Reset to exact data range
        tog_rangey = 0          ; Reset to exact data range
        load_frame
       End

     ; right-handed widgets:
     'selc_quan' : Begin
        cur_text = Widget_Info (id_selc_quan, /Combobox_Gettext)
        If (cur_text ne str_void) Then Begin
           iquantity = event.index
           get_iquantity;
           if(quant_vec[iquantity] eq 1) then get_vec_modifiers
           prep_data_2D
        Endif
     End

     'rbut_lilo' : Begin
        Widget_Control, id_rbut_lilo, Get_Value=tog_lilo
        tog_drange = 0            ; reset range to 'exact'
        prep_plot
     End

     'selc_vec_modi' : Begin
        cur_text = Widget_Info (id_selc_vec_modi, /Combobox_Gettext)
        If (cur_text ne str_void) Then Begin
           ivec_mod = event.index
           prep_data_2D
        Endif
     End

     'selc_modi' : Begin
        cur_text = Widget_Info (id_selc_modi, /Combobox_Gettext)
        If (cur_text ne str_void) Then Begin
           imodifier = event.index
;           prep_data
           prep_plot
        Endif
     End

     'rbut_drang' : Begin
        Widget_Control, id_rbut_drang, Get_Value=tog_drange
;        prep_data
     End
     'ffld_dran1' : Begin
        If (event.value lt plot_range[1]) Then Begin
           plot_range[0] = event.value
        Endif
        prep_plot
;        prep_data
     End
     'ffld_dran2' : Begin
        If (event.value gt plot_range[0]) Then Begin
           plot_range[1] = event.value
        Endif
        prep_plot
;        prep_data
     End

     'rbut_rangex' : Begin
        Widget_Control, id_rbut_rangex, Get_Value=tog_rangex
        if(tog_rangex eq 0) then begin
           get_xrange
        endif else Begin
        Endelse
        prep_plot
     End
     'ffld_ran1x' : Begin
        If (event.value lt xrange[1]) Then Begin
           xrange[0] = event.value
        Endif
        prep_plot
     End
     'ffld_ran2x' : Begin
        If (event.value gt xrange[0]) Then Begin
           xrange[1] = event.value
        Endif
        prep_plot
     End

     'rbut_rangey' : Begin
        Widget_Control, id_rbut_rangey, Get_Value=tog_rangey
        if(tog_rangey eq 0) Then Begin
           get_yrange
        Endif Else Begin
        Endelse
        prep_plot
     End
     'ffld_ran1y' : Begin
        If (event.value lt yrange[1]) Then Begin
           yrange[0] = event.value
        Endif
        print,event.value,yrange[1],yrange[0]
        prep_plot
        print,"Again: ",event.value,yrange[1],yrange[0]
     End
     'ffld_ran2y' : Begin
        If (event.value gt yrange[0]) Then Begin
           yrange[1] = event.value
        Endif
        prep_plot
     End

     'selc_cutp' : Begin ; Select place to cut
        icutplane = event.index
        prep_data_2D
     End
     'isld_cutp' : Begin
        Widget_Control, id_isld_cutp, Get_Value=cutlayer
        prep_data_2D
     End
     'rbut_geom' : Begin
        Widget_Control, id_rbut_geom, Get_Value=tog_geom
        prep_plot
     End
     'rbut_iLine' : Begin
        Widget_Control, id_rbut_iLine, Get_Value=tog_iLine
     End

     ; Plot settings:
     'selc_targ' : Begin
        tog_target = event.index
;        If (tog_target eq 1) Then Begin
;           tog_fram = 0         ; only single frames for Xwindow output
;        Endif
     End
     'text_fnam' : Begin
        Widget_Control, id_text_fnam, Get_Value=filebase
     End
     'npix_winx' : Begin
        Widget_Control, id_npix_winx, Get_Value=val
        win_size[0] = Min([2000, Max( [200, val ] ) ])
        If (val ne win_size[0]) Then alert, 'Xsize clipped to range [200,2000].'
     End
     'npix_winy' : Begin
        Widget_Control, id_npix_winy, Get_Value=val
        win_size[1] = Min([2000, Max( [200, val ] ) ])
        If (val ne win_size[1]) Then alert, 'Ysize clipped to range [200,2000].'
     End
     ; Change display ratio of x and y
     'rbut_xyrescale' : Begin
        Widget_Control, id_rbut_xyrescale, Get_Value=tog_xyrescale
        if(tog_xyrescale eq 0) then scale_y_to_x = 1.;
     End
     'ival_xyrescale' : Begin
        if(event.value gt 0.) then scale_y_to_x = event.value
     End
     ; Choose color-table:
     'isld_ctab' : Begin
        icoltab = event.value;
     End
     'ival_ctab' : Begin
        icoltab = Max([ event.value, 0            ]) ; not < 0
        icoltab = Min([ icoltab    , 54 ])   ; not > # of frames
     End
     ; Location / Suppression of colorbar
     'rbut_cbar' : Begin
        Widget_Control, id_rbut_cbar, Get_Value=tog_cbar
     End

     ; Final do-something buttons
     'butt_plot' : Begin
        Print, "Plotting quant:", iquantity
        do_plot
     End
     'butt_vset' : Begin
        alert, ""
        alert, "--Sectorplot (w) JK/RK ----"
        alert, "revision: " + rev_text
        alert, "data dir: " + maindir
        alert, "-----------------------"
        alert, ""
     End
     'butt_stop' : Begin
        alert, "Aborting for debugging."
        Stop
     End
     'butt_exit' : Begin
        If (w_is_open) Then Wdelete, 1
        Widget_Control, /Destroy, event.top
     End

  Endcase

  If (uvalue ne 'butt_plot') and (uvalue ne 'text_fnam') Then Begin
     filebase = proj_name + '_'        + $
                tab_quant[iquantity] ; restore output filename
     If (tog_lilo eq 1) Then filebase = filebase + '_Log'
  Endif

  widget_update
END


;===============================================================
;  Print info to message window
;===============================================================

PRO alert, text, valarray       ; print messages to text window
  @cronos_commons.h                ;  1st arg: plain tex
  If N_params() eq 2 Then Begin ;  2nd arg: array, gets printed inside []s
     val_string = ''
     For imem = 0, N_elements(valarray)-1 Do Begin
        val_string = val_string + StrCompress(valarray[imem])
     Endfor
     If N_elements(valarray) gt 1 Then Begin
        val_string = '[' + val_string + ']'
     Endif
     text = text + val_string
  Endif
  Widget_control, id_messages, Set_Value=text, /Append
END        

;===============================================================
;   Dependencies of widgets
;===============================================================
PRO widget_update
@cronos_commons.h

  sen_p = (proj_name    ne '')       ; flag "project has been selected"
  sen_f = (step_frame   ne -23)      ; flag "frame has been selected"
  sen_s = (tab_fields[0] ne str_void ) ; flag: fields determined
  sen_q = (tab_quant[0] ne str_void )  ; flag: quantities determined
  sen_v = (quant_vec[iquantity] eq 1)  ; Flag: is or is not vectorial
  sen_ps = (tog_target eq 0)           ; Flag: ps plot chosen
  sen_l = (frame_loaded ne -23)      ; flag "frame has been loaded"
  
; Projects box:
  Widget_Control, id_list_proj, Set_Value=proj_list
  Widget_Control, id_scan_frames, Sensitive=sen_p
  if (sen_p) Then Widget_Control, id_list_proj, Set_List_Select=proj_nseq 
  
; Frames box:
  Widget_Control, id_list_frames, Set_Value=frames_list
;  Widget_Control, id_scan_fields, Sensitive=sen_f
  if (sen_f) Then Widget_Control, id_list_frames, Set_List_Select=chosen_frame
  
; Load quantities box:
  tab_fieldsmenu = tab_fields
  Widget_Control, w3QLoad, Sensitive=sen_s
  Widget_Control, id_selc_quanLoad, Set_Value=tab_fieldsmenu
  Widget_Control, id_selc_quanLoad, Set_Combobox_Select=iquan_Load

; Choose quantity box:
  tab_quantmenu = tab_quant
  Widget_Control, w2Quant, Sensitive=sen_q
  Widget_Control, id_selc_quan, Set_Value=tab_quantmenu
  Widget_Control, id_selc_quan, Set_Combobox_Select=iquantity
  Widget_Control, id_rbut_lilo, Sensitive=flag_posdata, Set_Value=tog_lilo

; Modifier for vectorial quantities box:
  tab_vec_modifiermenu = tab_vec_modifier
  Widget_Control, w3Qvmod, Sensitive=sen_v
  Widget_Control, id_selc_vec_modi, Set_Value=tab_vec_modifiermenu
  Widget_Control, id_selc_vec_modi, Set_Combobox_Select=ivec_mod

; Data Range box
  Widget_Control, w2DRange, Sensitive=sen_f
  dranlabel = ['Data range:', 'Log data range:']
  Widget_Control, id_label_dran, Set_Value=dranlabel[tog_lilo]
  Widget_Control, id_rbut_drang, Set_Value=tog_drange
  Widget_Control, id_ffld_dran1, Set_Value=plot_range[0]
  Widget_Control, id_ffld_dran2, Set_Value=plot_range[1]
  If (tog_drange eq 0 or tog_drange eq 1) Then Begin
     Widget_Control, id_ffld_dran1, Sensitive=0
     Widget_Control, id_ffld_dran2, Sensitive=0
  Endif Else Begin
     Widget_Control, id_ffld_dran1, Sensitive=1
     Widget_Control, id_ffld_dran2, Sensitive=1
  Endelse
  
; Plot Range box
  Widget_Control, w2PRange, Sensitive=sen_f
  Widget_Control, id_rbut_rangex, Set_Value=tog_rangex
  Widget_Control, id_ffld_ran1x, Set_Value=xrange[0]
  Widget_Control, id_ffld_ran2x, Set_Value=xrange[1]
  Widget_Control, id_ffld_ran1x, Sensitive=tog_rangex
  Widget_Control, id_ffld_ran2x, Sensitive=tog_rangex

  Widget_Control, id_rbut_rangey, Set_Value=tog_rangey
  Widget_Control, id_ffld_ran1y, Set_Value=yrange[0]
  Widget_Control, id_ffld_ran2y, Set_Value=yrange[1]
  Widget_Control, id_ffld_ran1y, Sensitive=tog_rangey
  Widget_Control, id_ffld_ran2y, Sensitive=tog_rangey
  
  
; Select view box:
  Widget_Control, w2Views, Sensitive=sen_l
  Widget_Control, id_selc_cutp, Set_Combobox_Select=icutplane
  If (sen_s eq 1) Then Begin
     Widget_Control, id_isld_cutp, Set_Slider_Max=nx[2-icutplane]-2*n_gc-1
  Endif
  Widget_Control, id_rbut_geom, Set_Value=tog_geom
  Widget_Control, id_rbut_iLine, Sensitive=sen_s, Set_Value=tog_iLine
  
; Plot settings box:
  Widget_Control, w2Outpt, Sensitive=sen_l
  Widget_Control, id_selc_targ, Sensitive=sen_s
  Widget_Control, id_text_fnam, Sensitive=sen_ps
  Widget_Control, id_text_fnam, Set_Value=filebase
  Widget_Control, id_npix_winx, Set_Value=win_size[0]
  Widget_Control, id_npix_winy, Set_Value=win_size[1]
  Widget_Control, id_rbut_xyrescale, Set_Value=tog_xyrescale
  Widget_Control, id_ival_xyrescale, Set_Value=scale_y_to_x
  Widget_Control, id_ival_xyrescale, Sensitive=tog_xyrescale
  Widget_Control, id_isld_ctab, Sensitive=sen_s
  Widget_Control, id_isld_ctab, Set_Value=icoltab
  Widget_Control, id_isld_ctab, Set_Slider_Max=54
  Widget_Control, id_ival_ctab, Sensitive=sen_s
  Widget_Control, id_rbut_cbar, Set_Value=tog_cbar
  Widget_Control, id_rbut_cbar, Sensitive=sen_s

; do buttons
  Widget_Control, id_butt_plot, Sensitive=sen_l
End ; widget_update

;===============================================================
;  Main Program:
;===============================================================

@cronos_commons.h
settings
proj_scan

;---------------------------------------------------------------
; Build widgets
;---------------------------------------------------------------

; Declare position for widgets
w0Base  = Widget_Base (  /Row, Title='Cronos-Plot widget'+rev_text)
; Build two rows of sub-widgets
w1Lcol  = Widget_Base (w0Base, /Column)
w1Rcol  = Widget_Base (w0Base, /Column)

; Left-handed widgets
w2Proj  = Widget_Base (w1Lcol, /Column, /Frame)
w2Load  = Widget_Base (w1Lcol, /Row, /Frame)
w2Quant = Widget_Base (w1Lcol, /Column, /Frame)

; Left-handed sub-widgets
w3Frame = Widget_Base (w2Load, /Column, /Frame)
w3QLoad = Widget_Base (w2Load, /Column, /Frame)
; Quantity sub-widgets
w3Qvar = Widget_Base (w2Quant, /Row)
w3Qvmod = Widget_Base (w2Quant, /Row)
w3Qmod = Widget_Base (w2Quant, /Row)
;w3F_mod = Widget_Base (w2Frame, /Row)
;w3F_fr1 = Widget_Base (w2Frame, /Row)


; Right-handed widgets
w2DRange = Widget_Base (w1Rcol, /Column, /Frame)
w2PRange = Widget_Base (w1Rcol, /Row, /Frame)
w2Views = Widget_Base (w1Rcol, /Column, /Frame)
w2Outpt = Widget_Base (w1Rcol, /Column, /Frame) ; Output frame
w2DoBut = Widget_Base (w1Rcol, /Row   , /Frame) ; do plot frame

; Right-handed sub-widgets

; View sub-widgets
w3Views_cut = Widget_Base (w2Views, /Column)
w3V_cut = Widget_Base (w3Views_cut, /Row)
w3V_lay = Widget_Base (w3Views_cut, /Row)
w3V_geom = Widget_Base (w2Views, /Column)
w3Views_int = Widget_Base (w2Views, /Column)
w3V_int = Widget_Base (w3Views_int, /Row)
;w3V_ngc = Widget_Base (w3Views, /Row)

; Plot range subwidgets
w2PRangeX = Widget_Base (w2PRange, /Column, /Frame)
w2PRangeY = Widget_Base (w2PRange, /Column, /Frame)

; Plot setting subwidgets:
w3O_target = Widget_Base (w2Outpt, /Row)
w3O_fba = Widget_Base (w2Outpt, /Row)
w3O_wsize = Widget_Base (w2Outpt, /Row)
w3O_rescale = Widget_Base (w2Outpt, /Row)
w3O_coltit = Widget_Base (w2Outpt, /Row)
w3O_ctab = Widget_Base (w2Outpt, /Row)
w3O_cbar = Widget_Base (w2Outpt, /Row)


; Fill the left-handed widgets: (data handling and message window)

; Choose project widget (L)
Label_Proj  = Widget_Label   (w2Proj, Value='Projects:')
id_list_proj = Widget_List   (w2Proj, Value=proj_list, Ysize=8, $
                              Uvalue='selc_proj')
id_butt_scan = Widget_Button (w2Proj, Value='Re-scan folder', $
                              Uvalue='butt_scan')
id_scan_frames = Widget_Button (w2Proj, Value='Scan frames', $
                                Uvalue='butt_scan_frames')
;id_butt_load = Widget_Button (w2Proj, Value='Load (first) frame', $
;                              Uvalue='butt_load')

; Choose frame widget (L)
Label_Frame  = Widget_Label (w3Frame, value='Frame(s): timestep / time')
id_list_frames = Widget_List   (w3Frame, Value=frames_list, Ysize=8, $
                                Uvalue='selc_frame')
;id_scan_fields = Widget_Button (w3Frame, Value='Show fields', $
;                                Uvalue='butt_scan_fields')

; Choose frame widget (L)
;Label_Frame  = Widget_Label (w3F_mod, value='Frame(s):')
;id_isld_fra1 = Widget_Slider (w3F_fr1, Title='first', Scr_Xsize=200, $
;                              Uvalue='isld_fra1')
;id_ival_fra1 = CW_Field (w3F_fr1, Title='', Xsize=6,  /Return_Events, $
;                         Uvalue='ival_fra1')


; Choose quantity to load widget (L)
Label_QuantLoad = Widget_Label (w3QLoad, value='Pick field(s):')
id_selc_quanLoad  = Widget_Combobox (w3QLoad, Value=tab_quantLoad, $
                                 Uvalue='selc_quanLoad')
id_butt_load = Widget_Button (w3QLoad, Value='Load field(s)', $
                              Uvalue='butt_load')

; message widget (L)
id_messages = Widget_Text (w1Lcol, Xsize=40, Ysize=22, /Wrap, /Scroll)



; Quantity widget (R)
Label_Quant   = Widget_Label    (w3Qvar, Value='Quantity:')
id_selc_quan  = Widget_Combobox (w3Qvar, Value=tab_quant, $
                                 Uvalue='selc_quan')
Label_Modif  = Widget_Label     (w3Qvmod, Value='Mod vec field:')
id_selc_vec_modi = Widget_Combobox  (w3Qvmod, Value=tab_vec_modifier, $
                                 Uvalue='selc_vec_modi')
Label_Modif  = Widget_Label     (w3Qmod, Value='Modifier:')
id_selc_modi = Widget_Combobox  (w3Qmod, Value=tab_modifier, $
                                 Uvalue='selc_modi')
id_rbut_lilo = CW_Bgroup (w2Quant, ['linear','log'], Column=2, $
                          Label_left='scaling:', Uvalue = 'rbut_lilo', $
                          /Exclusive)

; Data Range widget (R)
id_label_dran = Widget_Label (w2DRange, Xsize=100, Value='Data range:')
id_rbut_drang = CW_BGroup    (w2DRange, ['exact','rounded','custom'], $
                             Column=3, Uvalue='rbut_drang', /Exclusive)
id_ffld_dran1 = CW_Field     (w2DRange, Title='min', Uvalue='ffld_dran1', $
                             /Floating, /Return_Events)
id_ffld_dran2 = CW_Field     (w2DRange, Title='max', Uvalue='ffld_dran2', $
                             /Floating, /Return_Events)

; Plot Range widget (R)
id_label_rangex = Widget_Label (w2PRangeX, Xsize=100, Value='x range:')
id_rbut_rangex = CW_BGroup    (w2PRangeX, ['exact','custom'], $
                               Column=2, Uvalue='rbut_rangex', /Exclusive)
id_ffld_ran1x = CW_Field     (w2PRangeX, Title='min', Uvalue='ffld_ran1x', $
                               /Floating, /Return_Events, Xsize=12)
id_ffld_ran2x = CW_Field     (w2PRangeX, Title='max', Uvalue='ffld_ran2x', $
                               /Floating, /Return_Events, Xsize=12)

id_label_rangey = Widget_Label (w2PRangeY, Xsize=100, Value='y range:')
id_rbut_rangey = CW_BGroup    (w2PRangeY, ['exact','custom'], $
                             Column=2, Uvalue='rbut_rangey', /Exclusive)
id_ffld_ran1y = CW_Field     (w2PRangeY, Title='min', Uvalue='ffld_ran1y', $
                             /Floating, /Return_Events, Xsize=12)
id_ffld_ran2y = CW_Field     (w2PRangeY, Title='max', Uvalue='ffld_ran2y', $
                             /Floating, /Return_Events, Xsize=12)


; View widget (R)
Label_Views  = Widget_Label    (w3V_cut, Value='Slice at')
id_selc_cutp = Widget_Combobox (w3V_cut, Value=tab_cutplane, $
                                Uvalue='selc_cutp')
Label_ViewE  = Widget_Label    (w3V_cut, Value='=')
id_isld_cutp = Widget_Slider   (w3V_cut, Value=cutlayer, $
                                Uvalue='isld_cutp')
id_rbut_geom = CW_BGroup        (w3V_geom, ['Cart.','Cyl.','Sph.'], Column=3, $
                                 Label_left='Data geometry:', $
                                 Uvalue='rbut_geom', /Exclusive)
Label_LineInt = Widget_Label    (w3V_int, Value='Line integral instead')
id_rbut_iLine  = CW_Bgroup (w3V_int, ['Off','On'], Column=2, $
                            Label_left='Switch:', Uvalue = 'rbut_iLine', $
                            /Exclusive)


; Plot settings widget (R)
Label_Target = Widget_LABEL    (w3O_target, Value='Output to')
id_selc_targ = Widget_Combobox (w3O_target, Value=['PostScript','Xwindow'], $
                                Uvalue='selc_targ')
id_text_fnam = CW_Field        (w3O_fba, Title='Filebase: ', $
                                /Return_Events, /String, $
                                Uvalue='text_fnam' )
id_npix_winx = CW_Field         (w3O_wsize, Title='Image size [px]:', $
                                 Uvalue='npix_winx', /Return_Events, $
                                 Xsize=8, /Integer)
id_npix_winy = CW_Field         (w3O_wsize, Title='x', $
                                 Uvalue='npix_winy', /Return_Events, $
                                 Xsize=8, /Integer)

; Scaling y to x
id_rbut_xyrescale = CW_Bgroup (w3O_rescale, ['Exact','Rescale'], Column=2, $
                               Label_left='dy/dx', Uvalue = 'rbut_xyrescale', $
                               /Exclusive)
id_ival_xyrescale = CW_Field (w3O_rescale, Title='Ratio', Xsize=6, $
                              /Return_Events, Uvalue='ival_xyrescale')

; Colortable choice: w3O_ctab
Label_Frame  = Widget_Label (w3O_coltit, value='Select Color-table:')
id_isld_ctab = Widget_Slider (w3O_ctab, Title='iCol', Scr_Xsize=200, $
                             Uvalue='isld_ctab')
id_ival_ctab = CW_Field (w3O_ctab, Title='', Xsize=6,  /Return_Events, $
                         Uvalue='ival_ctab')
id_rbut_cbar = CW_BGroup        (w3O_cbar, ['right','below','none'], Column=3, $
                                 Label_left='Type of color bar:', $
                                 Uvalue='rbut_cbar', /Exclusive)

; Final box: plot buttons and so on... (R)
id_butt_plot = Widget_Button (w2DoBut, Value='PLOT'         , $
                              Uvalue='butt_plot')
id_butt_vset = Widget_Button (w2DoBut, Value='View settings', $
                              Uvalue='butt_vset')
id_butt_stop = Widget_Button (w2DoBut, Value='Abort'        , $
                              Uvalue='butt_stop')
id_butt_exit = Widget_Button (w2DoBut, Value='Exit'         , $
                              Uvalue='butt_exit')


; Display the widget:
Widget_Control, w0Base, /Realize, /Append

widget_update

Xmanager, 'ptest', w0Base, /No_Block

END ; of Main

;      !P.Font=1
;      DEVICE, SET_FONT='Times', /TT_FONT
;      DEVICE, SET_FONT='-adobe-courier-bold-o-normal--11-80-100-100-m-60-iso8859-1'
;      DEVICE,
;      SET_FONT='-adobe-courier-medium-r-normal--14-100-100-100-m-90-iso8859-1'
