;@cplot_coms
COMMON COLORS, R_orig, G_orig, B_orig, R_curr, G_curr, B_curr 


coltab = icoltab

nlevels     = 100                ;
bottom     = 100                ;
nlevels = 253
nColors = 253
bottom = 1
if(coltab eq 24) then begin
    nlevels = 60                 ;
    nColors = 60
    bottom = 127                ;
endif
if(coltab eq 37) then begin
    nlevels = 253
    nColors = 253
    bottom = 1
endif
if(coltab eq 41) then begin
    nlevels = 135                ;
    bottom = 120                ;
endif
if(coltab eq 42 or coltab eq 43 or coltab eq 44 or coltab eq 45 or coltab eq 46 or coltab eq 47 or coltab eq 48 or coltab eq 50 or coltab eq 51) then begin
    nlevels = 253                ;
    bottom = 1                  ;
endif

if(coltab eq 52) then begin
	nlevels = 253
	bottom = 1
endif


;        bottom     = 155;127;100        ;
dcol       = (ommax - ommin)/nlevels ;
userLevels = IndGen(nlevels)*dcol + ommin
userColors = Indgen(nlevels)+bottom

; Here bottom is the first color-index to use -
; in this case we assign nlevels 99-199 with new colors
; from the table -  the others are kept as bw.
; i. e. we only use 101 colors from coltab
if(coltab le 40) then begin
    loadct,coltab,ncolors=nlevels,bottom=bottom
endif
if(coltab eq 16) then begin
;    loadct,coltab,ncolors=nlevels+1,bottom=bottom-1
;    loadct,coltab,ncolors=nlevels+1,bottom=bottom+1
    loadct,coltab,ncolors=nlevels+1,bottom=bottom-1
    r_curr(1) = 0
    g_curr(1) = 0
    b_curr(1) = 0
    r_curr(0) = 0
    g_curr(0) = 0
    b_curr(0) = 0
    tvlct,r_curr,g_curr,b_curr
endif
if(coltab eq 24) then begin
    loadct,coltab,ncolors=253,bottom=1
endif
if(coltab eq 41) then begin
    loadct,4,ncolors=253,bottom=1
endif
if(coltab eq 42) then begin
    restore,filename='temp5coltab.sav'
    tvlct,red,green,blue
endif
if(coltab eq 43) then begin
    restore,filename='densitycoltab5.sav'
    tvlct,red,green,blue
endif

if(coltab eq 44) then begin
    restore,filename='colsternmod.sav'
    tvlct,red,green,blue
endif

if(coltab eq 45) then begin
    restore,filename='whiteblack.sav'
    tvlct,red,green,blue
endif

if(coltab eq 46) then begin
    restore,filename='whiteblackalt.sav'
    tvlct,red,green,blue
endif

if(coltab eq 47) then begin
;    restore,filename='coltabtempniceinv.sav'
    restore,filename='whiteblackalt.sav'
    tvlct,red,green,blue
endif

if(coltab eq 48) then begin
    restore,filename='bluetobrown02.sav'
    tvlct,red,green,blue
endif

if(coltab eq 50) then begin
    restore,filename='bluetobrownWithLine.sav'
    tvlct,red,green,blue
endif

if(coltab eq 51) then begin
    restore,filename='fromredtogreen.sav'
    tvlct,red,green,blue
endif

if(coltab eq 52) then begin
    restore,filename='fromredtogreenThreeLines.sav'
    tvlct,red,green,blue
endif

if(coltab eq 53) then begin
    restore,filename='testtableBunt.sav'
    tvlct,red,green,blue
 endif

if(coltab eq 54) then begin
    restore,filename='testtableBunt2.sav'
    tvlct,red,green,blue
endif

nColors = nlevels
