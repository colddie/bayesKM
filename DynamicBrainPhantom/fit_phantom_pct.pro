; turn off delay estimation temporally
test              = 0
debug             = 0
presmooth         = 1
usemotion         = 0
simulate_interval = 1
simulate_AIF      = 0
has_delay         = 0
noisefree         = 0     ;0-noisy
use_svd           = 0

; restore, filename='tmp.sav'
; restore,filename='tmp_delay.sav'
; if usemotion then restore, filename='tmp_delay_transform0.sav'

if has_delay eq 0 then begin
    if noisefree eq 1 then  restore,filename='simulatedphantom_nodelay.sav'
    if noisefree eq 0 then  restore,filename='simulatephantom_motionfree_nodelay.sav'
    if usemotion eq 1 then  restore,filename='simulatephantom_moving_nodelay.sav'
endif else begin
    if noisefree eq 1 then  restore, filename='simulatedphantom_delay.sav'
    if noisefree eq 0 then  restore, filename='simulatedphantom_motionfree_delay.sav'
    if usemotion eq 1 then  restore, filename='simulatedphantom_moving_delay.sav'
endelse

if presmooth then begin
  for iframe = 0, (size(reconall))[4]-1 do begin 
    for iplane = 0, (size(reconall))[3]-1 do begin 
        reconall[*,*,iplane,iframe]= niconvolgauss(reconall[*,*,iplane,iframe],fwhm=2.0,dimensions=[0,1])
    endfor
  endfor
endif

; restore,filename='tmp_delay_new.sav'
; reconall =tacall 


; cbfall_ref = cbfall
; mttall_ref = mttall
; delayall_ref = delayall
baselineall = reconall[*,*,*,0]

dt       = 1    ; %s
addskull = 1
tstart   = 0 ; %s
tstop    = 49 ; %s
upsample = 1
frameNr  = (tstop-tstart)/dt+1
t        = indgen(frameNr)*dt;     ;0:1:49
to       = indgen(25)*2 +1  ;1:2:49;
aifs     = [0, 0, 0, 0, 25, 105, 220, 350, 440, 485, 430, 300, 180, 110, 104, 108, 115, 125, 115, 108, 98, 90, 98, 108, 112];
ts       = indgen(tstop*10+1)*0.1   ;0:0.1:49;
ts       = congrid(ts,n_elements(ts)/upsample)
aif      = interpol(float(aifs),to,ts, /spline);
nplane   = 256;
ncol     = 256;
nrow     = 256;

; if simulate_AIF then begin
;     aif = shift(aif,20)
;     ; aif = aif + randomu(seed,n_elements(aif))*mean(aif)/2
; endif

; goto,line0
if use_svd then goto, line1
; if use_pmap then goto, line2


line0:

isweight      = 1
def_pmin      = [0.0,0.00001,0.0]     ; cbf, mtt; cbv and delay are calculated 
def_pmax      = [100.0,100.0,0.0]  
doSD          = 1
doCL          = 0
bootstrapIter = 200  ; has to be larger than 100!
if has_delay eq 1 then def_pmax[2] = 12.0

nplane      = 256;
ncol        = 256;
nrow        = 256;
cbfall      = fltarr(ncol,nrow,nplane)
mttall      = cbfall*0
delayall    = cbfall*0
cbfall_sd   = cbfall*0
mttall_sd   = cbfall*0
delayall_sd = cbfall*0

for iplane = 151,151 do begin   ;153,157 0, nplane-1 do begin
    print, 'Processing plane number '+nistring(iplane)
    for irow = 0, nrow-1 do begin
        print, 'Processing row number '+nistring(irow)
        for icol = 0, ncol-1 do begin

        if test then begin
            icol   = 156 ;79
            irow   = 135 ;126
            iplane = 154
            icol   = 128 ;79
            irow   = 173 ;126
            iplane = 154
        endif
        
        tmp = reform(double(reconall[icol,irow,iplane,*]))
        tac = reform(double(reconall[icol,irow,iplane,*]-baselineall[icol,irow,iplane]))
        if total(tmp) lt 0 or total(tac) lt 1e-5 or tac[10] eq tac[20] then continue
        ; tac1 = dblarr(25)
        ; for j=0,24 do $
        ;     tac1[j] = tac[j*2+1]
        ; too = indgen(50) +1 

        if simulate_interval then begin
            undersample = 2   ;4
            t     = indgen(frameNr)*dt;  
            index = indgen(n_elements(tac)/undersample)*undersample
            tac   = tac[index]
            t     = t[index]      
        endif

        if simulate_AIF then begin
            truncate = 23   ;30
            t     = indgen(frameNr)*dt;  
            index = where(t lt truncate)
            tac   = tac[index]
            t     = t[index]   
            index = where(ts lt truncate)
            aif   = aif[index]
            ts    = ts[index]
        endif

        tac1 =  interpol(tac,t,ts, /spline);

        lib = '/home/tsun/bin/tpcclib-master/build/bin/libmtga_idl.so'
        matrix = double(fltarr(bootstrapIter*n_elements(def_pmin))) ; change with num_param
        output = double(fltarr(8))
        weights = fltarr(n_elements(ts))+1.0
        to[0] = 0. ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;??????????????
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
        cbfall[icol,irow,iplane]      = output[0]   ;cbf
        mttall[icol,irow,iplane]      = output[1]   ;mtt
        delayall[icol,irow,iplane]    = output[2]   ;delay
        cbfall_sd[icol,irow,iplane]   = output[5]   ;cbf sd
        mttall_sd[icol,irow,iplane]   = output[6]   ;mtt sd 
        delayall_sd[icol,irow,iplane] = output[7]   ;delay sd       
        ; cbv = cbf*mtt/60;
        ; ttp = where(tac eq max(tac));

        if test then  stop
        endfor
    endfor
    cbvall = cbfall * mttall/60;
    

endfor


if usemotion eq 1 then $
save,filename='results_simulation_moving.sav', cbfall,mttall,cbvall,delayall,cbfall_sd,mttall_sd,delayall_sd,aif,ts
stop
if simulate_interval eq 1 then $
save,filename='results_simulation_interval.sav', cbfall,mttall,cbvall,delayall,cbfall_sd,mttall_sd,delayall_sd,aif,ts

if simulate_AIF eq 1 then $
save,filename='results_simulation_aif.sav', cbfall,mttall,cbvall,delayall,cbfall_sd,mttall_sd,delayall_sd,aif,ts

if simulate_interval eq 0 and simulate_AIF eq 0 then $
save,filename='results_simulation.sav', cbfall,mttall,cbvall,delayall,cbfall_sd,mttall_sd,delayall_sd,aif,ts





line1:
; ---------- compare with block-svd ------------
cbfall1    = fltarr(ncol,nrow,nplane)
mttall1    = cbfall1*0
delayall1  = cbfall1*0
inmap      = reconall *0
for iframe = 0, frameNr-1 do $
     inmap[*,*,*,iframe] = reconall[*,*,*,iframe]-baselineall

lambda = 0.2   ; truncation
mpad   = 2
mask   = 0
dt    /= 1.    ; /=10.
first  = 0
last   = 49    ; 200     ; depends ont mtt
; tac = congrid(tac,10)
; aif = congrid(aif,10)


;--- plane based compromise between speed and memory ----
inmap  = inmap[*,*,151:152,*]    ;153:157
nframe = (size(inmap))[4]
index  = ts*0
if simulate_interval then begin
    undersample = 2   ;4
    t     = indgen(frameNr)*dt;  
    index = indgen(frameNr/undersample)*undersample 
    t     = t[index]    
    index = indgen(n_elements(ts)/undersample)*undersample 
    ts    = ts[index]
    aif   = aif[index]
endif
if simulate_AIF then begin
    truncate = 23   ;30
    t     = indgen(frameNr)*dt;  
    index = where(t lt truncate)
    t     = t[index]   
    index = where(ts lt truncate)
    ts    = ts[index]
    aif   = aif[index]
endif
inmap1 = fltarr(ncol,nrow,(size(inmap))[3],n_elements(index)) 

for iplane=0,0 do begin   ;0, nplane-1 do begin
    for irow=0, nrow-1 do begin
        for icol=0, ncol-1 do begin
        tac = inmap[icol,irow,iplane,*]
        
        if simulate_interval then begin 
            index = indgen(n_elements(tac)/undersample)*undersample
            tac   = tac[index]  
        endif
        if simulate_AIF then begin
            index = where(t lt truncate)
            tac   = tac[index] 
        endif

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
cbvmap /= 10

;---- volume based, fast but with large memory ------
; aifss = interpol(float(aifs),to,t, /spline);
; pct_bsvd, inmap,aifss,dt,lambda,mpad,mask, $
;         cbf=cbfmap,cbv=cbvmap,mtt=mttmap,delay=delaymap, k=k, $
;         first=first,last=last


;----- voxel based, slow with little memory
; for iplane=153,157 do begin   ;0, nplane-1 do begin
;     print, 'Processing Plane '+ nistring(iplane)+' uisng b-SVD...'
;     for irow=0, nrow-1 do begin
;         for icol=0, ncol-1 do begin
;         tac = inmap[icol,irow,iplane,*]
;         tac1 =  interpol(tac,t,ts, /spline);
;         if total(tac) lt 1e-5 or tac[10] eq tac[20] then continue
;     pct_bsvd, tac1,aif,dt,lambda,mpad,mask, $
;                 cbf=cbfmap,cbv=cbvmap,mtt=mttmap,delay=delaymap, k=k, $
;                 first=first,last=last
; ; cbvmap = cbfmap * mttmap;
;         cbfall1[icol,irow,iplane] = cbfmap
;         mttall1[icol,irow,iplane] = mttmap/10
;         cbvall1[icol,irow,iplane] = cbvmap
;         delayall1[icol,irow,iplane] = delaymap
;         ; print, cbfmap, mttmap/10., cbvmap, delaymap, ttp/10.
;         endfor
;     endfor
; endfor

; not sure why cbf 10 times larger, so is cbv, mtt seems to be fine, fixed
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
inmap  = reconall *0
cbvmap = cbfall_ref*mttall_ref/60
for iframe = 0, frameNr-1 do $
     inmap[*,*,*,iframe] = reconall[*,*,*,iframe]-baselineall

; inmap = inmap[*,*,*,153:157]
pct_bbbp, inmap, cbvmap, aifss, dt, rho, PRE_bbbp, POST_bbbp, mask,   bbbmap, x, ymap, R 
; BBBP    - A map of brain permeability in mL/100g/min [Y x X]
; X       - The independent variable of the patlak plot [T x 1]
; MAP    - A map of dependent variables of the patlak plot [T x Y x X]
; R       - A map of coefficients of determination (R^2)


;
stop

End