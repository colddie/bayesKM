; create aif 0..49 s in 0.1 sec sampling interval
; to = 1:2:49;
; aifs = [0 0 0 0 25 105 220 350 440 485 430 300 180 110 104 108 115 125 115 108 98 90 98 108 112];
; ts = 0:0.1:49;
; aif = interp1(to,aifs,ts,'spline');



; phantom from external 
; simple model patlak, Ki,Vp
; dump = call_external()



; 2-compartment model, K1, k2, Dv
; nplane = 60
; ncol = 80
; nrow = 80

; phantomK1 = fltarr(ncol,nrow,nplane)
; phantomk2 = fltarr(ncol,nrow,nplane)
; phantomk3 = fltarr(ncol,nrow,nplane)
; phantomk4 = fltarr(ncol,nrow,nplane)
; 0.6,  1.4,  0.06, 0.002, 0.15

fromphantom = 1
; img=readraw('/storage0/home/tsun/code/idl/mgh/fabber_pet/fdg_0_1.bin',256,256,190,'float')
img=readraw('/home/tsun/work/fabber_pet/fdg_0_1.bin',256,256,190,'float')

size    = size(img)
ncol    = size[1]
nrow    = size[2]
nplane  = size[3]
phantom = img
phantom[where(img gt 0. and img lt 50.)]=1.
phantom[where(img ge 50. and img lt 100.)]=2.
phantom[where(img ge 150. and img lt 250.)]=3.
phantom[where(img ge 100. and img lt 150.)]=3.
phantom[where(img ge 150. and img lt 250.)]=4.

phantomK1 = phantom / (max(phantom)/0.6)
phantomk2 = phantom / (max(phantom)/1.4)
phantomk3 = phantom / (max(phantom)/0.06)
phantomk4 = phantom / (max(phantom)/0.002)
phantomm  = phantom / (max(phantom)/0.15)*0

; genereate dynamic frames

;  fdg_3comp_basic, 'init', a, b, c, d
;  fdg_3comp_basic, 'start_parms', a, b, c, d, parms
;  fdg_3comp_basic, 'calc', A, B, C, D, parms

  ; Ask which model to be used.
  ;------------------
  modeldef = 'fdg_3comp_basic'

  plasma_t = $
  [ 0.166667,    0.283333,    0.366667,   0.466667,    0.550000,    0.616667,$ 
    0.716667,    0.833333,    0.950000,    1.06667,     1.18333,     1.31667,$
    1.46667,     1.60000,     1.78333,     1.95000,     3.70000,     6.76667,$
    11.8833,     20.1167,     32.0333 ,    45.5833,     57.8333 ]

  if 1 then begin
   plasma_c = $
   [9.99229e-06, 1.10063e-05, 0.000347761,   0.0495032,  0.343662,  0.486564, $
    0.604856,   0.428336,    0.105233,  0.0904635,   0.0633207,   0.0558369, $
    0.0506417,   0.0459184,   0.04085060,  0.0381968,   0.030804,   0.0342774,$
    0.0305163 ,  0.0296231,   0.0243044,   0.0203284,   0.0173419] - 0.012 > 0
   endif else begin
    plasma_c = $
    [9.99229e-06, 1.10063e-05, 0.000347761,  0.0495032,  0.343662,  0.486564, $
     0.214856,   0.0968336,    0.105233,  0.0904635,   0.0733207,  0.0698369, $
     0.0686417,   0.0679184,   0.0685060, 0.0681968,   0.0546804,   0.0442774,$
     0.0365163 ,  0.0296231,   0.0243044, 0.0203284,   0.0173419]
   endelse
 
  ; Decide, based on the model, what to do with the spill-over curves
  ;------------------
  spill_t = plasma_t * 0
  spill_c = plasma_c * 0
  if total(modeldef eq ['fdg_3comp_spill', 'acetate_2comp_metab', $
                        'acetate_2comp_efmetab', 'nh3_3comp_basic', $
                        'nh3_3comp_kul', 'nh3_2comp_kul']) gt 0 $
    then begin
     spill1_t = plasma_t
     spill1_c = plasma_c
     spill2_t = spill_t
     spill2_c = spill_c
  endif else begin
     spill1_t = plasma_t
     spill1_c = plasma_c
  endelse

  ; Generate times for scan duration and mid scan times.
  ;--------------------
  if 0 then begin
    scan_dur  = float([.1, .1, .1, .1, .1, .1, .15, .15, .15, .15, $
                       .5, .5, .5, .5, 2,2,2,2,2, 5,5,5,5,5,5,5,5])
  endif else begin
    scan_dur = (plasma_t[1] - plasma_t[0]) + 0 * plasma_t
  endelse
  scan_time = scan_dur * 0
  for i = 1, n_elements(scan_time)-1 do begin
    scan_time(i) = scan_time(i-1) + float(scan_dur(i-1))
  end
  mid_time = scan_time + ( 0.5 * scan_dur )  ; [min]

  totaldur = string(max(mid_time))
;   read, totaldur, prompt='Total scan duration in min [' + totaldur + ']: '
  if totaldur ne '' then begin
    totaldur = float(totaldur)
    mid_time = mid_time(where(mid_time le totaldur))
    scan_dur = scan_dur(where(mid_time le totaldur))
  endif
  tissue_t = mid_time



; -----------------------------------------
 tissue_t = plasma_t
; -----------------------------------------


imgall = fltarr(ncol, nrow, nplane, n_elements(plasma_t))
for iplane = 100, 100 do begin;nplane-1 do begin
    print, 'plane number ' + nistring(iplane)
    for irow = 0,  nrow-1 do begin
        for icol = 0,  ncol-1 do begin
    ; for irow = 100, 100 do begin ; nrow-1 do begin
    ;     for icol = 100, 100 do begin ; ncol-1 do begin

        ; Generate default parms and ask the user if he wants to 
        ; modify them
        ;-------------------
        if fromphantom then begin
            parms = [phantomK1[icol,irow,iplane],phantomk2[icol,irow,iplane],phantomk3[icol,irow,iplane],phantomk4[icol,irow,iplane],phantomm[icol,irow,iplane]]
        endif else $
            call_procedure, modeldef, 'start_parms', parms

        if total(parms) eq 0. then continue

        call_procedure, modeldef, 'parms', descrip
        ; for i = 0, n_elements(descrip)-1 do $
            ; print, descrip(i), parms(i), format = "(A, ' = ', F8.4)"
        ;   keuze =''
        ;   read, keuze, prompt="Do you want to change k's [y/N]: "
        ;   if keuze eq 'y' then begin
        ;     for i = 0, n_elements(descrip)-1 do begin
        ;       prompt = string(descrip(i), parms(i), $
        ;                  format = "(A, ' = ', F8.4, ': ')")
        ;       nieuw = ''
        ;       read, nieuw, prompt = prompt
        ;       if nieuw ne '' then begin
        ;         parms(i) = float(nieuw)
        ;       endif
        ;     endfor
        ;   print, ''
        ;   for i = 0, n_elements(descrip)-1 do $
        ;     print, descrip(i), parms(i), format = "(A, ' = ', F8.4)"
        ;   endif

        ; Compute the tissue counts
        ;--------------------------
        tissue_c = NIcalc_outcurve(modeldef, parms, plasma_t, plasma_c, tissue_t, $
                                spill1_t = spill1_t, spill1_c = spill1_c, $
                                spill2_t = spill2_t, spill2_c = spill2_c)

        
        tissue_cc = double(tissue_c *0)
        dummy = call_external('libtpccm.so', 'simC2_idl', double(plasma_t), double(plasma_c), fix(n_elements(plasma_t)), $
                             double(parms[0]),double(parms[1]),double(parms[2]),double(parms[3]), tissue_cc,return_type=2,/verbose)
          ; ret=simC2(input.x, input.c[0].y, input.sampleNr, 0.4, 0.5, 0.4, 0.3, 
          ;   tissue.c[0].y, NULL, NULL);
        ; stop

        imgall[icol,irow,iplane,*] = tissue_cc

        endfor
    endfor
endfor


stop

; validate






; Ask which model to be used.
fitphantomK1 = fltarr(ncol,nrow,nplane)
fitphantomk2 =  fitphantomK1 *0
fitphantomk3 =  fitphantomK1 *0
fitphantomk4 =  fitphantomK1 *0
 
  ; Decide, based on the model, what to do with the spill-over curves
  ;------------------
  if total(modeldef eq ['fdg_3comp_spill', 'acetate_2comp_metab', $
                        'acetate_2comp_efmetab', $
                        'nh3_3comp_basic', 'nh3_3comp_kul', 'nh3_2comp_kul']) gt 0 $
    then begin
     spill1_t = plasma_t
     spill1_c = plasma_c
     spill2_t = spill_t
     spill2_c = spill_c
  endif else begin
     spill1_t = spill_t
     spill1_c = spill_c
  endelse

  ; Read startt and stopt if requested
  ;---
    startt = 0.0
    stopt  = max(tissue_t)
    startstring = NIstring(startt)
    stopstring  = NIstring(stopt)


  ; Start params?
  ;---
    parms   = NImodeler(modeldef, 'start_parms')
    parms[4] = 0

  ; Ask which parameters to fix
  ;----------------
    descrip = NImodeler(modeldef, 'parms')
    fix_parms = intarr(n_elements(descrip))
    ; fix_parms[n_elements(descrip)-1] = 1


  ; Apply the model
  ;----------------
for iplane = 100, 100 do begin;nplane-1 do begin
    print, 'plane number ' + nistring(iplane)
    for irow = 0,  nrow-1 do begin
        for icol = 0,  ncol-1 do begin
    ; for irow = 100, 100 do begin ; nrow-1 do begin
    ;     for icol = 100, 100 do begin ; ncol-1 do begin

        tissue_c = imgall[icol,irow,iplane,*]
        
        if total(tissue_c) lt 1e-3 then continue
        fitparms = NImodeler(modeldef, plasma_t, plasma_c, tissue_t, tissue_c, $
                    sigmaa = sigma, start_parms = parms, $     ;/plotcurve,  
                    weight=tissue_t / sqrt(tissue_c > max(tissue_c / 5)), $
                    spill1_t = spill1_t, spill1_c = spill1_c, $
                    spill2_t = spill2_t, spill2_c = spill2_c, $
                    startt = startt, stopt = stopt, fix_parms = fix_parms)

        descrip = NImodeler(modeldef, 'parms')
        for i = 0, n_elements(descrip)-1 do $
            ; print, descrip(i), fitparms(i), sigma(i), $
            ;     format = "(A, ' = ', F8.4, ' +- ', F8.4)"


            fitphantomK1[icol,irow,iplane] = fitparms[0]
            fitphantomk2[icol,irow,iplane] = fitparms[1]
            fitphantomk3[icol,irow,iplane] = fitparms[2]
            fitphantomk4[icol,irow,iplane] = fitparms[3]

        endfor
    endfor
endfor


stop

End