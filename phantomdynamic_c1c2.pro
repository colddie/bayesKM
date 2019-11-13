; create aif 0..49 s in 0.1 sec sampling interval
; to = 1:2:49;
; aifs = [0 0 0 0 25 105 220 350 440 485 430 300 180 110 104 108 115 125 115 108 98 90 98 108 112];
; ts = 0:0.1:49;
; aif = interp1(to,aifs,ts,'spline');


pad = '/home/tsun/bin/fsl/install/src/fabber_core/fabber_pet_c1/'
fromphantom = 1
modeldef = 'fdg_3comp_basic'
tracer = 'fdg'    ; way, recropide, fdg, fdopa, fluropirdize
isotope = 'F-18'     ;'C-11'   ; 'F-18'
compartment = 'C2'    ; 'C1', 'C2', 'srtm', 'rtcm'
bound = 0.05*0             ; high the high variance in Ks
tacnoiselevel = 0.005     ; higher the high noise
nframe =  158   ; by default, deprecated
rebin_fac = 4.0; 0.5
fwhm = 5.0

case tracer of 
  'fluropirdize' : begin
    ; img=readraw('/storage0/home/tsun/code/idl/mgh/fabber_pet/fdg_0_1.bin',256,256,190,'float')
    img=readraw('/home/tsun/bin/fsl/install/src/fabber_core/fabber_pet_c1/fdg_1.bin',256,256,190,'float')
    size    = size(img)
    ncol    = size[1]
    nrow    = size[2]
    nplane  = size[3]
    phantom = img
    phantom[where(img gt 0. and img lt 50.)]=1.
    phantom[where(img ge 50. and img lt 100.)]=2.
    phantom[where(img ge 100. and img lt 150.)]=3.
    phantom[where(img ge 150. and img lt 250.)]=4.
    phantomK1 = phantom * 0
    phantomk2 = phantom * 0
    phantomk3 = phantom * 0
    phantomk4 = phantom * 0
    phantomm = phantom * 0

    ; According to 
    phantomK1 = phantom / (max(phantom)/0.12)
    phantomk2 = phantom / (max(phantom)/0.2)

    ; genereate dynamic frames

      ; Ask which model to be used.
      ;------------------
      ; tracer = 'F-18'
      addnoise = 0    ; simple noise
  end

  'fdg' : begin
  ; img=readraw('/storage0/home/tsun/code/idl/mgh/fabber_pet/fdg_0_1.bin',256,256,190,'float')
  img=readraw('/home/tsun/bin/fsl/install/src/fabber_core/fabber_pet_c1/fdg_1.bin',256,256,190,'float')
  size    = size(img)
  ncol    = size[1]
  nrow    = size[2]
  nplane  = size[3]
  phantom = img
  phantom[where(img gt 0. and img lt 50.)]   =1.
  phantom[where(img ge 50. and img lt 100.)] =2.
  phantom[where(img ge 100. and img lt 150.)]=3.
  phantom[where(img ge 150. and img lt 250.)]=4.
  phantomK1 = phantom * 0
  phantomk2 = phantom * 0
  phantomk3 = phantom * 0
  phantomk4 = phantom * 0
  phantomm  = phantom * 0

  ; According to "Dynamic Positron Emission Tomography Image Restoration via a Kinetics-Induced Bilateral Filter"
  phantomK1 = phantom / (max(phantom)/0.12)
  phantomk2 = phantom / (max(phantom)/0.2)
  phantomk3 = phantom / (max(phantom)/0.1)
  phantomk4 = phantom / (max(phantom)/0.01)
  phantomm  = phantom / (max(phantom)/0.15)*0

  ; genereate dynamic frames

  ;  fdg_3comp_basic, 'init', a, b, c, d
  ;  fdg_3comp_basic, 'start_parms', a, b, c, d, parms
  ;  fdg_3comp_basic, 'calc', A, B, C, D, parms

    ; Ask which model to be used.
    ;------------------
    ; tracer = 'F-18'
    addnoise = 0    ; simple noise
  end

  'way' : begin
    img=readraw('/home/tsun/bin/fsl/install/src/fabber_core/fabber_pet_c1/wayc11_1.bin',256,256,190,'float')
    size    = size(img)
    ncol    = size[1]
    nrow    = size[2]
    nplane  = size[3]
    phantom = img
    phantom[where(img gt 0. and img le 50.)]=1.
    phantom[where(img gt 50. and img lt 100.)]=2.
    phantom[where(img ge 100. and img lt 150.)]=3.
    phantom[where(img ge 150. and img lt 250.)]=4.
    phantomK1 = phantom * 0
    phantomk2 = phantom * 0
    phantomk3 = phantom * 0
    phantomk4 = phantom * 0
    phantomm = phantom * 0

    ; According to The PET Radioligand [carbonyl-11C]Desmethyl-WAY-100635 Binds to 5-HT1A Receptors, JNM2002
    if 0 then begin
        ; frontal, ant, temporal, insular lobes
        phantomK1[where(phantom eq 4)] = 0.14
        phantomk2[where(phantom eq 4)] = 0.25
        phantomk3[where(phantom eq 4)] = 0.18
        phantomk4[where(phantom eq 4)] = 0.04
        phantomm[where(phantom eq 4)] = 0.0
        ; hippocampus
        phantomK1[where(phantom eq 3)] = 0.11
        phantomk2[where(phantom eq 3)] = 0.15
        phantomk3[where(phantom eq 3)] = 0.21
        phantomk4[where(phantom eq 3)] = 0.04
        phantomm[where(phantom eq 3)] = 0.0
        ; raphe nucleus
        phantomK1[where(phantom eq 2)] = 0.12
        phantomk2[where(phantom eq 2)] = 0.54
        phantomk3[where(phantom eq 2)] = 0.24
        phantomk4[where(phantom eq 2)] = 0.05
        phantomm[where(phantom eq 2)] = 0.0
        ; cerebellum
        phantomK1[where(phantom eq 1)] = 0.16
        phantomk2[where(phantom eq 1)] = 0.33
        phantomk3[where(phantom eq 1)] = 0.03 *0         ;set zero to use simplified model
        phantomk4[where(phantom eq 1)] = 0.04 *0          
        phantomm[where(phantom eq 1)] = 0.0
    endif else begin
        ; bound = 0.8   ; 50% variance
        ; if you want to use reference model, two assumptions should be made
        ; 1.reference region has no specific binding 2. K1/k2 is same everywhere
        meanC = [ [0.14,0.29,0.18,0.05,0.0], $      ; 2tcm
                [0.11,0.225,0.21,0.05,0.0], $
                [0.22,0.45,0.24,0.02,0.0], $
                [0.16,0.33,0.0,0.0,0.0] ]

        if compartment eq 'srtm' then begin        ; srtm, do not generate same tacs as full comaprtment!
        ; k2a = k2/(1+BPND)
        ;
        meanC = [ [0.875,0.063,0.27,0.0,0.0], $
                [0.6875,0.043,0.15,0.0,0.0], $
                [1.375,0.034,2.0,0.0,0.0], $
                [1.,0.1,0.0,0.0,0.0] ]     ;was 0.33
        endif
        
        if compartment eq 'rtcm' then begin         ; rtcm
        meanC = [ [0.875,0.29,0.18,0.05,0.0], $
                [0.6875,0.225,0.21,0.05,0.0], $
                [1.375,0.45,0.24,0.02,0.0], $
                [1.,0.33,0.0,0.0,0.0] ]
        endif

        num = [n_elements(where(phantom eq 4)),n_elements(where(phantom eq 3)), $
                n_elements(where(phantom eq 2)),n_elements(where(phantom eq 1))]
        idx = [4,3,2,1]
        for i = 0,n_elements(idx)-1 do begin
          ; lb+(ub-lb)*randomu(seed,100) 
          phantomK1[where(phantom eq idx[i])] = meanC[0,i]*(1-bound+bound*2*randomu(seed,num[i]) ) 
          phantomk2[where(phantom eq idx[i])] = meanC[1,i]*(1-bound+bound*2*randomu(seed,num[i]) )   
          phantomk3[where(phantom eq idx[i])] = meanC[2,i]*(1-bound+bound*2*randomu(seed,num[i]) )   
          phantomk4[where(phantom eq idx[i])] = meanC[3,i]*(1-bound+bound*2*randomu(seed,num[i]) )
        endfor         
    endelse
    ; tracer = 'C-11'
  end

  'fdopa' : print, 'not implemented yet>!'

  'recropide' : begin
    img=readraw('/home/tsun/bin/fsl/install/src/fabber_core/fabber_pet_c1/recropide_1.bin',256,256,190,'float')
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
    phantomK1 = phantom * 0
    phantomk2 = phantom * 0
    phantomk3 = phantom * 0
    phantomk4 = phantom * 0
    phantomm = phantom * 0

    ; According to Quantitative PET Analysis of the Dopamine D2 Receptor Agonist Radioligand, JNM2009
    if 0 then begin
        ; putamen
        phantomK1[where(phantom eq 4)] = 0.44
        phantomk2[where(phantom eq 4)] = 0.07
        phantomk3[where(phantom eq 4)] = 0.15
        phantomk4[where(phantom eq 4)] = 0.19
        phantomm[where(phantom eq 4)] = 0.0
        ; caudate nucleus
        phantomK1[where(phantom eq 3)] = 0.39
        phantomk2[where(phantom eq 3)] = 0.06
        phantomk3[where(phantom eq 3)] = 0.11
        phantomk4[where(phantom eq 3)] = 0.20
        phantomm[where(phantom eq 3)] = 0.0
        ; brain stem /thalamus
        phantomK1[where(phantom eq 2)] = 0.43
        phantomk2[where(phantom eq 2)] = 0.07
        phantomk3[where(phantom eq 2)] = 0.03
        phantomk4[where(phantom eq 2)] = 0.13
        phantomm[where(phantom eq 2)] = 0.0
        ; cerebellum, grey, white matters
        phantomK1[where(phantom eq 1)] = 0.41
        phantomk2[where(phantom eq 1)] = 0.06
        phantomk3[where(phantom eq 1)] = 0.0
        phantomk4[where(phantom eq 1)] = 0.0
        phantomm[where(phantom eq 1)] = 0.0
    endif else begin
        ; bound = 0.5   ; 50% variance
        meanC =  [[0.44,0.07,0.15,0.19,0.0], $
                [0.39,0.06,0.11,0.20,0.0], $
                [0.43,0.07,0.03,0.13,0.0], $
                [0.41,0.06,0.0,0.0,0.0]]
        num = [n_elements(where(phantom eq 4)),n_elements(where(phantom eq 3)), $
                n_elements(where(phantom eq 2)),n_elements(where(phantom eq 1))]
        idx = [4,3,2,1]
        for i = 0,n_elements(idx)-1 do begin
          ; lb+(ub-lb)*randomu(seed,100) 
          phantomK1[where(phantom eq idx[i])] = meanC[0,i]*(1-bound+bound*2*randomu(seed,num[i]) ) 
          phantomk2[where(phantom eq idx[i])] = meanC[1,i]*(1-bound+bound*2*randomu(seed,num[i]) )   
          phantomk3[where(phantom eq idx[i])] = meanC[2,i]*(1-bound+bound*2*randomu(seed,num[i]) )   
          phantomk4[where(phantom eq idx[i])] = meanC[3,i]*(1-bound+bound*2*randomu(seed,num[i]) )
        endfor         
    endelse
    ; compartment = 'C2'
    ; tracer = 'C-11'
  end

endcase

; stop


  ; plasma_t = $
  ; [ 0.166667,    0.283333,    0.366667,   0.466667,    0.550000,    0.616667,$ 
  ;   0.716667,    0.833333,    0.950000,    1.06667,     1.18333,     1.31667,$
  ;   1.46667,     1.60000,     1.78333,     1.95000,     3.70000,     6.76667,$
  ;   11.8833,     20.1167,     32.0333 ,    45.5833,     57.8333 ]

  ; if 0 then begin
  ;  plasma_c = $
  ;  [9.99229e-06, 1.10063e-05, 0.000347761,   0.0495032,  0.343662,  0.486564, $
  ;   0.604856,   0.428336,    0.105233,  0.0904635,   0.0633207,   0.0558369, $
  ;   0.0506417,   0.0459184,   0.04085060,  0.0381968,   0.030804,   0.0342774,$
  ;   0.0305163 ,  0.0296231,   0.0243044,   0.0203284,   0.0173419] - 0.012 > 0
  ;  endif else begin
  ;   plasma_c = $
  ;   [9.99229e-06, 1.10063e-05, 0.000347761,  0.0495032,  0.343662,  0.486564, $
  ;    0.214856,   0.0968336,    0.105233,  0.0904635,   0.0733207,  0.0698369, $
  ;    0.0686417,   0.0679184,   0.0685060, 0.0681968,   0.0546804,   0.0442774,$
  ;    0.0365163 ,  0.0296231,   0.0243044, 0.0203284,   0.0173419]
  ;  endelse

  ; Read input function from file
  ;------------------
  inputfile = pad+'input.dat'
  dataStruct = { plasma_t:0.0, plasma_c:0.0}
  nrows = File_Lines(inputfile)
  data = Replicate(dataStruct, nrows)
  OpenR, lun, inputfile, /GET_LUN
  ReadF, lun, data
  Free_Lun, lun
  plasma_t = (data.plasma_t)   ;[0:nframe-1]   ; truncate equilibrium
  plasma_c = (data.plasma_c)   ;[0:nframe-1]
  plasma_t = congrid(plasma_t,n_elements(plasma_t)/rebin_fac,/center,/interp)
    ; rebin(plasma_t,n_elements(plasma_t)/rebin_fac)
  plasma_c = congrid(plasma_c,n_elements(plasma_c)/rebin_fac,/center,/interp)
    ; rebin(plasma_c,n_elements(plasma_c)/rebin_fac)

tosave=pad+tracer+compartment+'rebin'+nistring(rebin_fac)
file_mkdir, tosave
cd, tosave


  if compartment eq 'srtm' or compartment eq 'rtcm' then begin
    tmp = read_ascii('ref_tissuec.txt')
    plasma_c = reform(tmp.FIELD1,n_elements(tmp.FIELD1)) 
    plasma_c = plasma_c[0:n_elements(plasma_t)-1]
    ; plasma_c[0] = 0.
  endif
  

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

;   ; Generate times for scan duration and mid scan times.
;   ;--------------------
;   if 0 then begin
;     scan_dur  = float([.1, .1, .1, .1, .1, .1, .15, .15, .15, .15, $
;                        .5, .5, .5, .5, 2,2,2,2,2, 5,5,5,5,5,5,5,5])
;   endif else begin
;     scan_dur = (plasma_t[1] - plasma_t[0]) + 0 * plasma_t
;   endelse
;   scan_time = scan_dur * 0
;   for i = 1, n_elements(scan_time)-1 do begin
;     scan_time(i) = scan_time(i-1) + float(scan_dur(i-1))
;   end
;   mid_time = scan_time + ( 0.5 * scan_dur )  ; [min]

;   totaldur = string(max(mid_time))
; ;   read, totaldur, prompt='Total scan duration in min [' + totaldur + ']: '
;   if totaldur ne '' then begin
;     totaldur = float(totaldur)
;     mid_time = mid_time(where(mid_time le totaldur))
;     scan_dur = scan_dur(where(mid_time le totaldur))
;   endif
;   tissue_t = mid_time


; stop

; -----------------------------------------
 tissue_t = plasma_t
; -----------------------------------------


; model the resolution in parametric space
if fwhm gt 0 then begin
      phantomK1 = NIconvolgauss(phantomK1,fwhm=fwhm)     ;,dimensions=3
      phantomk2 = NIconvolgauss(phantomk2,fwhm=fwhm)     ;,dimensions=3
      phantomk3 = NIconvolgauss(phantomk3,fwhm=fwhm)     ;,dimensions=3
      phantomk4 = NIconvolgauss(phantomk4,fwhm=fwhm)     ;,dimensions=3
end

imgall = fltarr(ncol, nrow, nplane, n_elements(plasma_t))
for iplane = 82, 83 do begin;nplane-1 do begin
    print, 'plane number ' + nistring(iplane)
    for irow = 0,  nrow-1 do begin
        for icol = 0,  ncol-1 do begin
    ; for irow = 100, 100 do begin ; nrow-1 do begin
    ;     for icol = 100, 100 do begin ; ncol-1 do begin

        ; Generate default parms and ask the user if he wants to 
        ; modify them
        ;-------------------
        if fromphantom then begin
            parms = [phantomK1[icol,irow,iplane],phantomk2[icol,irow,iplane], $
                     phantomk3[icol,irow,iplane],phantomk4[icol,irow,iplane], $
                     phantomm[icol,irow,iplane]]
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
        if 0 then begin
          tissue_c = NIcalc_outcurve(modeldef, parms, plasma_t, plasma_c, tissue_t, $
                                  spill1_t = spill1_t, spill1_c = spill1_c, $
                                  spill2_t = spill2_t, spill2_c = spill2_c)
        endif
        
        tissue_cc = double(plasma_c *0)
        if compartment eq 'C1' then $
          dummy = call_external('libtpccm.so', 'simC1_idl', double(plasma_t), double(plasma_c), fix(n_elements(plasma_t)), $
                              double(parms[0]),double(parms[1]), tissue_cc,return_type=2,/verbose)

        if compartment eq 'C2' then $
          dummy = call_external('libtpccm.so', 'simC2_idl', double(plasma_t), double(plasma_c), fix(n_elements(plasma_t)), $
                              double(parms[0]),double(parms[1]),double(parms[2]),double(parms[3]), tissue_cc,return_type=2,/verbose)
          ; ret=simC2(input.x, input.c[0].y, input.sampleNr, 0.4, 0.5, 0.4, 0.3, 
          ;   tissue.c[0].y, NULL, NULL);

        if compartment eq 'srtm' then $
          dummy = call_external('libtpccm.so', 'simSRTM_idl', double(plasma_t), double(plasma_c), fix(n_elements(plasma_t)), $
                              double(parms[0]),double(parms[1]),double(parms[2]), tissue_cc,return_type=2,/verbose)

        if compartment eq 'rtcm' then $
          dummy = call_external('libtpccm.so', 'simRTCM_idl', double(plasma_t), double(plasma_c), fix(n_elements(plasma_t)), $
                              double(parms[0]),double(parms[1]),double(parms[2]),double(parms[3]), tissue_cc,return_type=2,/verbose)

        ;  addnoise = 1
        ; if addnoise then begin
        ;   rand = randomn(SEED,n_elements(plasma_t))
        ;   scale = mean(tissue_cc)/mean(rand) *rand /20
        ;   tissue_cc += scale
        ; endif
      

        imgall[icol,irow,iplane,*] = tissue_cc

        endfor
    endfor
endfor

; if fwhm gt 0 then begin
;     for iframe=0,nframe/2-1 do begin
;       tmp = NIconvolgauss(imgall[*,*,*,iframe],fwhm=fwhm)     ;,dimensions=3
;       imgall[*,*,*,iframe] = tmp
;     end 
; end
tmpimg = imgall[*,*,82:83,*]
mask = niread_nifti(pad+'mask.nii')
for i=0,n_elements(plasma_c)-1 do tmpimg[*,*,*,i] *= mask
padname = tosave+'/'+'actimg_'+ tracer+compartment
save, filename= padname +'.sav', tmpimg
help,niwrite_nii(tmpimg, padname +'.nii')


; write the plasma_t.txt and plasma_c.txt
if compartment eq 'C1' or compartment eq 'C2' then begin
  fname= tosave+'/'+'plasma_t.txt'
  OPENW,1,fname 
  PRINTF,1,plasma_t
  CLOSE,1
  fname= tosave+'/'+'plasma_c.txt'
  OPENW,1,fname 
  PRINTF,1,plasma_c
  CLOSE,1

  fname = padname +'.nii' ; pad+'_'+nistring(bound)+'_'+nistring(tacnoiselevel)+'_noise.nii'
  mask = phantom *0
  mask[where(phantom eq 1)] = 1.0
  mask=mask[*,*,82:83]
  img = niread_nii(fname,orientation='RAS') 
  tissue_c = fltarr(n_elements(plasma_t))
  for ip = 0, n_elements(plasma_t)-1 do begin
    tmp = img[*,*,*,ip] * mask
    tissue_c[ip] = mean(tmp[where(tmp gt 0.)])
  endfor
  fname = tosave+'/'+'ref_tissuec.txt'
  OPENW,1,fname 
  PRINTF,1,tissue_c
  CLOSE,1

  ; ; write the interpolated plasma_t.txt and plasma_c.txt
  ; plasma_tt = congrid(plasma_t,2*n_elements(plasma_t))
  ; plasma_cc = interpol(float(plasma_c),plasma_t,plasma_tt, /spline);

  ; fname= pad+'plasma_tt.txt'
  ; OPENW,1,fname 
  ; PRINTF,1,plasma_tt
  ; CLOSE,1
  ; fname= pad+'plasma_cc.txt'
  ; OPENW,1,fname 
  ; PRINTF,1,plasma_cc
  ; CLOSE,1

endif





; nifti header
plasma_t_sif = [0]
nsample = n_elements(tissue_cc)
; for i = 0, nsample-1 do $
plasma_t[0] = 0.00001
plasma_t_sif = [plasma_t_sif, plasma_t]

plasma_t_sif= [ [plasma_t_sif[0:nsample-1]], [plasma_t_sif[1:nsample]], $
               [fltarr(nsample)], [fltarr(nsample)] ]

fname= padname + '.sif' 
OPENW,1,fname 
PRINTF,1, '01/01/1970 00:00:00 ' + nistring(n_elements(plasma_t)) + ' 4 1 . '+isotope
PRINTF,1,transpose(plasma_t_sif),FORMAT='(F,1X,F,1X,I9,1X,F7.2)' 
CLOSE,1
stop

; insert counts based on images and then add noise
; imgweigh actimg_way.nii actimg_way.sif
; fvar4img -i=C-11 actimg_way.nii 0.005 actimg_way_noise.nii
spawn, 'imgweigh '+padname+'.nii '+padname+'.sif' 
spawn, 'cp ' +padname+'.sif ' +padname+ 'tmp.sif'  

inputfile = padname+'tmp.sif' 
dataStruct = { plasma_t:0.0, plasma_c:0.0, weight:0L, unkown:0.0}
nrows = File_Lines(inputfile)-1
data = Replicate(dataStruct, nrows)
OpenR, lun, inputfile, /GET_LUN
SKIP_LUN, lun, 1, /LINES
ReadF, lun, data
Free_Lun, lun
plasma_t_sif[*,2] = data.weight
plasma_t_sif[0,1] = 0.0000001
plasma_t_sif[1,0] = 0.0000001
fname= padname + '.sif' 
OPENW,1,fname 
PRINTF,1, '01/01/1970 00:00:00 ' + nistring(n_elements(plasma_t)) + ' 4 1 . '+isotope
PRINTF,1,transpose(plasma_t_sif),FORMAT='(F,1X,F,1X,I9,1X,F7.2)' 
CLOSE,1
; check if the start and end time-stamps are identical 
spawn, 'fvar4img -i=' +isotope+' '+padname+'.nii '+nistring(tacnoiselevel)+' '+padname+'_'+nistring(bound)+'_'+nistring(tacnoiselevel)+'_noise.nii'






; plot time-activity-curve
tac_cere=tmpimg[128,128,0,*] 
tac_matter=tmpimg[116,76,0,*]   
tac_nuc=tmpimg[125,172,0,*]  
window,0
plot,plasma_t,plasma_c  
oplot,plasma_t,tac_cere,col=1 
oplot,plasma_t,tac_matter,col=2   ; green
oplot,plasma_t,tac_nuc,col=3    
window, 1
plot, plasma_t, tissue_c, col=4  
stop


















;;;;;;;;;;; validate



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
for iplane = 82,83 do begin;nplane-1 do begin
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