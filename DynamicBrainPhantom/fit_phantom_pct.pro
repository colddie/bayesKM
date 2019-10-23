; turn off delay estimation temporally
test = 0
debug = 0
usemotion = 1


; restore, filename='tmp.sav'
restore,filename='tmp_delay.sav'
if usemotion then restore, filename='tmp_delay_transform0.sav'
cbfall_ref = cbfall
mttall_ref = mttall
delayall_ref = delayall

dt       = 1    ; %s
addskull = 1
tstart   = 0 ; %s
tstop    = 49 ; %s
frameNr = (tstop-tstart)/dt+1
t = indgen(frameNr)*dt;     ;0:1:49
to = indgen(25)*2 +1  ;1:2:49;
aifs = [0, 0, 0, 0, 25, 105, 220, 350, 440, 485, 430, 300, 180, 110, 104, 108, 115, 125, 115, 108, 98, 90, 98, 108, 112];
ts = indgen(tstop*10+1)*0.1   ;0:0.1:49;
aif = interpol(float(aifs),to,ts, /spline);


isweight = 1
def_pmin = [0.0,0.00001,0.0]     ; cbf, mtt; cbv and ttp are calculated 
def_pmax = [100.0,100.0,0.0]  
doSD = 1
doCL = 0
bootstrapIter = 200  ; has to be larger than 100!

nplane = 256;
ncol   = 256;
nrow   = 256;
cbfall    = fltarr(ncol,nrow,nplane)
mttall    = cbfall*0
delayall  = cbfall*0
cbfall_sd = cbfall*0
mttall_sd = cbfall*0
delayall_sd = cbfall*0
for iplane=153,157 do begin   ;0, nplane-1 do begin
    print, 'Processing plane number '+nistring(iplane)
    for irow=0, nrow-1 do begin
        for icol=0, ncol-1 do begin

        if test then begin
            icol = 156 ;79
            irow = 135 ;126
            iplane = 154
            icol = 128 ;79
            irow = 173 ;126
            iplane = 154
        endif

        tac = reform(double(tacall[icol,irow,iplane,*]-baselineall[icol,irow,iplane]))
        if total(tac) lt 1e-5 or tac[10] eq tac[20] then continue
        ; tac1 = dblarr(25)
        ; for j=0,24 do $
        ;     tac1[j] = tac[j*2+1]
        ; too = indgen(50) +1 
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
    
    

endfor


stop


End