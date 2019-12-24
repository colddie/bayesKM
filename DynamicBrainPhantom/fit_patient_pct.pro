; turn off delay estimation temporally
test  = 0
debug = 0


restore,filename = '/home/tsun/data/CTperfusion_sydney.sav'
workimg = regimgs      ;frame
nframe = (size(workimg))[4];
nplane = (size(workimg))[3];
ncol   = (size(workimg))[2];
nrow   = (size(workimg))[1];
baseline   = NIimgmu2hounsfield(workimg[*,*,*,0],/h2m) * 10
tacall     = workimg
for iframe=0,nframe-1 do begin   
    tacall[*,*,*,iframe] = NIimgmu2hounsfield(tacall[*,*,*,iframe],/h2m) *10
    tacall[*,*,*,iframe] -= baseline
    tacall[*,*,*,iframe] = tacall[*,*,*,iframe]  > 0
endfor

dosmooth = 1
if dosmooth then  print, 'test 3d aif smooth...'
cbfall    = fltarr(ncol,nrow,nplane)
mttall    = cbfall*0
delayall  = cbfall*0
cbfall_sd = cbfall*0
mttall_sd = cbfall*0
delayall_sd = cbfall*0



;;;;;;;;;;;;;;;;;;;
;; ------------------------------
; read in all time stamps
a=fltarr(512,512,726)
tstamps = []
for i=1,726  do begin
    a(*,*,i-1)= niread_dicom2('/home/tsun/data/Siemens_recons/0610697/'+'IM'+string(i,format='(I05)'),/reset,/singlefile, header=header1)
    tstamp = float(header1(1,where(header1(0,*) eq 'FRAME_START_TIME_MS'))) / 1000.
    if i mod 22 eq 0 then  tstamps = [tstamps, [tstamp]]
endfor

tstamps -= tstamps[0]
tdiff = (tstamps-shift(tstamps,1)) 
tdiff[0] = 0.0



dt       = mean(tdiff)  ;;60./33    ;?
addskull = 1
tstart   = min(tstamps)    ;0 ; 
tstop    = max(tstamps)    ;32   ; ?
frameNr = (tstop-tstart)/dt
t = tstamps   ;indgen(frameNr)*dt;     ;0:1:49
to = indgen(25)*2 +1  ;1:2:49;
aifs = [0, 0, 0, 0, 25, 105, 220, 350, 440, 485, 430, 300, 180, 110, 104, 108, 115, 125, 115, 108, 98, 90, 98, 108, 112];
to = to[where(to le (tstop-tstart))]
aifs = aifs[where(to le (tstop-tstart))]
sample = 1.
ts = indgen((tstop-tstart)*sample+1)/sample   ;0:0.1:49;             ?????????????????????????????????????????
aif = interpol(float(aifs),to,ts)    ;, /spline);



; goto, line1
stop

;;;;;;;;;;;;;;;;;;;
;; ------------------------------
isweight = 1
def_pmin = [-100.0,0.0, 0.0]    ;[0.0,0.00001,0.0]     ; cbf, mtt; cbv and ttp are calculated 
def_pmax = [100.0,100.,100.0]    ;[100.0,100.0,0.0]  
doSD = 1
doCL = 0
bootstrapIter = 200  ; has to be larger than 100!


for iplane=11,12 do begin   ;0, nplane-1 do begin
    print, 'Processing plane number '+nistring(iplane)
    ; for irow=0, nrow-1 do begin
    ; for irow=0, 99 do begin  
    ; for irow=100, 199 do begin
    ; for irow=200, 299 do begin
    ; for irow=300, 399 do begin  
    for irow=400, 499 do begin   
    print, 'ROW ' + nistring(irow) + '...'
        for icol=0, ncol-1 do begin

        if test then begin
            icol = 256 ;79
            irow = 256 ;126
            iplane = 11
        endif

        tac = reform(double(tacall[icol,irow,iplane,*]))     ;-baseline[icol,irow,iplane]))
        if total(abs(tac)) lt 1e-5 or tac[10] eq tac[20] then continue
        ; tac1 = dblarr(25)
        ; for j=0,24 do $
        ;     tac1[j] = tac[j*2+1]
        ; too = indgen(50) +1 
        tac1 =  interpol(tac,t,ts)             ;, /spline);

        lib = '/home/tsun/bin/tpcclib-master/build/bin/libmtga_idl.so'
        matrix = double(fltarr(bootstrapIter*n_elements(def_pmin))) ; change with num_param
        output = double(fltarr(8))
        weights = fltarr(n_elements(ts))+1.0
        success = call_external(lib,'pCT_idl',long(n_elements(ts)),double(ts),tac1, $
                        double(aif),output,long(debug),$
                        long(isweight),double(weights),double(def_pmin),double(def_pmax), $
                        long(doSD),long(doCL), $
                        long(bootstrapIter),matrix) 

        ; coarse sampling
        if 0 then begin
            weights = fltarr(n_elements(t))+1.0
            aifss = interpol(float(aifs),to,t, /spline);
            to[0] = 0. ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;??????????????
            success = call_external(lib,'pCT_idl',long(n_elements(t)),double(t),tac, $
                        double(aifss),output,long(debug),$
                        long(isweight),double(weights),double(def_pmin),double(def_pmax), $
                        long(doSD),long(doCL), $
                        long(bootstrapIter),matrix) 
        end

        ; stop
        cbfall[icol,irow,iplane] = output[0]   ;cbf
        mttall[icol,irow,iplane] = output[1]   ;mtt
        delayall[icol,irow,iplane] = output[2]   ;delay
        cbfall_sd[icol,irow,iplane] = output[5]   ;cbf sd
        mttall_sd[icol,irow,iplane] = output[6]   ;mtt sd 
        delayall_sd[icol,irow,iplane] = output[7]   ;delay sd       
        ; cbv = cbf*mtt/60;
        ; ttp = where(tac eq max(tac));

        if test then  stop
        endfor
    endfor
    cbvall = cbfall * mttall/60;
    

endfor



stop




line1:
; ---------- compare with block-svd ------------
cbfall1    = fltarr(ncol,nrow,nplane)
mttall1    = cbfall*0
delayall1  = cbfall*0
inmap     = tacall *0
for iframe = 0, frameNr-1 do $
     inmap[*,*,*,iframe] = tacall[*,*,*,iframe]     ;-baseline

lambda = 0.2   ; truncation
mpad   = 2
mask   = 0
dt    /= 1    ; /=10.
first  = 0
last   = 32    ; 200     ; depends ont mtt
; tac = congrid(tac,10)
; aif = congrid(aif,10)


;--- plane based compromise between speed and memory ----
inmap = inmap[*,*,11:12,*]
inmap1 = fltarr(ncol,nrow,2,n_elements(ts))
for iplane=0,1 do begin   ;0, nplane-1 do begin
    for irow=0, nrow-1 do begin
        for icol=0, ncol-1 do begin
        tac = inmap[icol,irow,iplane,*]
        tac1 =  interpol(tac,t,ts, /spline);
        inmap1[icol,irow,iplane,*] = tac1
        endfor
    endfor
endfor
; tac = inmap1[icol,irow,iplane,*]
; if total(tac) lt 1e-5 or tac[10] eq tac[20] then continue
pct_bsvd, inmap1,aif,dt,lambda,mpad,mask, $
    cbf=cbfmap,cbv=cbvmap,mtt=mttmap,delay=delaymap, k=k, $
    first=first,last=last
mttmap /= 10



stop



line2:
; ---------- calculate permeability map -----------

; PCT Parameters
rho = 1.05;   %Average brain tissue density
PRE_bbbp = 20; %First frame in BBBP calculation
POST_bbbp = 49; %Last frame in BBBP calculation
POST_bbbp = min([POST_bbbp,frameNr]); %Last frame for BBBP cannot exceed the actual last frame
dt = 1; % time interval between samples in seconds


aifss = interpol(float(aifs),to,t, /spline);
inmap  = tacall *0
cbvmap = cbfall_ref*mttall_ref/60
for iframe = 0, frameNr-1 do $
     inmap[*,*,*,iframe] = tacall[*,*,*,iframe]-baselineall

; inmap = inmap[*,*,*,153:157]
pct_bbbp, inmap, cbvmap, aifss, dt, rho, PRE_bbbp, POST_bbbp, mask,   bbbmap, x, ymap, R 
; BBBP    - A map of brain permeability in mL/100g/min [Y x X]
; X       - The independent variable of the patlak plot [T x 1]
; MAP    - A map of dependent variables of the patlak plot [T x Y x X]
; R       - A map of coefficients of determination (R^2)









stop
End