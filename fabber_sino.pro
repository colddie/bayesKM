fit_mode = 1;   ;0 - fabber, 1 - regular nlls fitting
t0 = systime(1)

; read dynamic FDG phantom generated from phantomdynamic_c1c2.pro
pad = '/home/tsun/bin/fsl/install/src/fabber_core/fabber_pet_c1/'
; tmpimg = niread_nii(pad+'actimg_fdgC2_0.00000_0.00500000_noise.nii')
; tmpimg = niread_nii(pad+'actimg_fdgC2.nii')
restore, filename = pad+'actimg_fdgC2.sav'

stop


; ===================================================================
; ================generate simulated  noise-free projections==============
; ===================================================================
; 
; modified from mmr_ocl_rayt_test.pro
;
 test_fwd  = 0
 test_back = 0 
 use_mltr  = 0 
 fast_math = 1
 force_cpu = 1
 rigmotion = 0
 addnoise  = 0
 usempi    = 1

; pixelsizemm = 1.0 ;1.5 ;3.0; 1.36719
; planesepmm  = 1.0 ;1.5 ;3.0; 3.0
oversample = 0
; raytracing = 1
; maxmu_mm = 0.05

sumdet = 2
nrcols = 172*2 /  sumdet                ;344;
nrrows = 172*2 /  sumdet               ;344;
nrplanes = 2
img_dim = [nrcols, nrrows, nrplanes]
pixelsizemm =  [2.0445, 2.0445, 2.03125] ;[2.08625, 2.08625, 2.03125] ;[3.06675, 3.06675, 2.03125];  recon pixel size
binsize = [2.0445, 2.0445, 2.03125*2]     ; if span eq 1 then planesep is 4.0625 mm
relreconsize = pixelsizemm[0] / binsize[0]
relreconplanesize = pixelsizemm[2] / binsize[2]
volumefwhm =  1.0    ;    4.5 / [2.0445, 2.03125]
; projd_s1  = nidef_proj(/pet3d_mmr, span = 1, $
;                                         relreconsize = relreconsize, $
;                                         relreconplanesize = relreconplanesize , $
;                                         nrcols = img_dim[0])
projd0 = nidef_proj(/pet3d_mmr, span = 11,$
                                        relreconsize = relreconsize, $
                                        relreconplanesize = relreconplanesize *2,$
                                        nrcols = nrcols, nrplanes=nrplanes,$
                                        volumefwhm = volumefwhm, sumdet=sumdet) 
                   


nframes = (size(tmpimg))[4]

; ===================================================================
goto, line2     ; skip the rest and go directly to fitting
; ===================================================================


tmpsino = fltarr(172,252,837,nframes)

  for iframe = 0, nframes-1 do begin 
    print, 'forward projecting frame ' + nistring(iframe)
    img = congrid(tmpimg[*,*,*,iframe], nrcols,nrrows,nrplanes)
    sino = 0
    NIproj, img, sino, projd = projd0
    if addnoise then $                 ;add Poisson noise
      sino  = nipoisson(seed,sino)  
    tmpsino[*,*,*,iframe] = sino
    ; stop
  endfor 
print, 'TIme elapsed '+ nistring(systime(1)-t0)+' seconds...'
stop

; ===================================================================
; ================ image Reconstruction ==============
; ===================================================================

line1:
restore,filename='filetmpsino.sav'
print, '-- Start Reconstruction --'
  time0 = systime(1)
if test_fwd eq 1 then begin
  ; a test for projection
  sino_fwd = 0
  NIproj, img, sino_fwd, projd = projd0;, subset = 0
  
endif else if test_back eq 1 then begin
  ; a test for back-projection
  img_back = 0
  NIproj, img_back, sino, projd = projd0, /back;, subset = 110

endif else begin
  ; start reconstruction
  nriter = 1
  nrsub  = 32

  tmprecon = fltarr(nrcols,nrrows,nrplanes,nframes)
  for iframe = 0, nframes-1 do begin 
    print, 'MLEM w/o motion at frame '+nistring(iframe)
    recon = 0
    sino = tmpsino[*,*,*,iframe]
    recon = NImaposem(sino, projd = projd0, nriter=nriter, nrsub=nrsub,$
        /norm_per_sub) 
    tmprecon[*,*,*,iframe] = recon
  endfor
endelse
print, 'Recon time (min): ', (systime(1)-time0)/60.


stop


; ===================================================================
; ================ image NLLS fitting ==============
; ===================================================================
;compare to the truth Ki=K1*k3/(k2+k3)
line2:
restore,filename='filetmprecon.sav'
; tmprecon = fltarr(nrcols,nrrows,nrplanes,nframes)
; for iframe=0, nframes-1 do $
;   tmprecon[*,*,*,iframe]=congrid(tmpimg[*,*,*,iframe],nrcols,nrrows,nrplanes)

fitmode = 'nlls' ;'mcmc'



nframe=158
inputfile = pad+'input.dat'
dataStruct = { plasma_t:0.0, plasma_c:0.0}
nrows = File_Lines(inputfile)
data = Replicate(dataStruct, nrows)
OpenR, lun, inputfile, /GET_LUN
ReadF, lun, data
Free_Lun, lun
plasma_t = (data.plasma_t)[0:nframe-1]   ; truncate equilibrium
plasma_c = (data.plasma_c)[0:nframe-1]
plasma_t = rebin(plasma_t,n_elements(plasma_t)/2)
plasma_c = rebin(plasma_c,n_elements(plasma_c)/2)
plasma_t0 = [0,plasma_t[0:n_elements(plasma_t)-2]]

frameNr= n_elements(plasma_t0)
tstart = 70
tstop  = 120
output = double(fltarr(5)) 
debug  = 0
llsq_model = 0
isweight   = 0;
patlog     = '/home/tsun/bin/tpcclib-master/build/bin/libmtga_idl.so'
tmprecon   = reform(tmprecon,long(nrcols)*nrrows*nrplanes,nframes)
Kiimg  = fltarr(nrcols,nrrows,nrplanes)
Vbimg  = Kiimg*0
Kiimgsd= Kiimg*0
Vbimgsd= Kiimg*0
SWSSimg= Kiimg*0
for jvoxel = 0, long(nrcols)*nrrows*nrplanes-1 do begin
    tac = reform(tmprecon[jvoxel, *])
    if total(tac) lt 80. then continue ;;n_elements(where(tac ne 0)) le 1 then continue
    if (jvoxel mod 1000 eq 0) then print,'processed '+nistring(jvoxel)+' voxels...'

    if fitmode eq 'nlls' then begin
    ; do fitting
      success = call_external(patlog, 'patlak_idl', long(frameNr), double(plasma_t0), $
                  double(plasma_t), double(tac), $
                  double(plasma_c),double(tstart),double(tstop),double(output),long(debug),long(llsq_model),$
                  long(isweight)) 
      Kiimg[jvoxel]   = output[0]
      Kiimgsd[jvoxel] = output[1]
      Vbimg[jvoxel]   = output[2]
      Vbimgsd[jvoxel] = output[3]
      SWSSimg[jvoxel] = output[4]
    endif
; stop
    if fitmode eq 'mcmc' then begin
      ; do sampling
      mcmc  = 0
      useprior = 0
      debug = 0
      model = 5  
      initialK = [6e-3,2e-7,0.,0.]
      lb = [0.,0.,0.,0.]
      ub = [1e-2, 1e-6, 0.,0.]
      rwmh_par_scale = double(0.016) 
      hmc_step_size  = double(0.001)
      rwmh_n_burnin  = 10000L *8          ; ?
      rwmh_n_draws   = 10000L *2
      weight = (plasma_t -shift(plasma_t,1)) / plasma_t[n_elements(plasma_t)-1]
      weight[0] = 0.
      output = fltarr(rwmh_n_draws, n_elements(initialK))
      prior = dblarr(n_elements(initialK))
      mcmc_tac = '/home/tsun/bin/mcmc-master/tests/example/librwmh_tac.so'
      success = call_external(mcmc_tac, 'rwmh_tac_2tpc', long(frameNr), double(tac), double(plasma_t), double(plasma_c), $
                      double(weight), double(prior), output, $
                      double(initialK), double(lb), double(ub), $
                      double(rwmh_par_scale),double(hmc_step_size),long(rwmh_n_burnin),long(rwmh_n_draws),long(model),$
                      long(debug),long(mcmc),long(useprior),$
                      double(plasma_t0), double(tstart), double(tstop)) 

      Kiimg[jvoxel] = mean(output(*,0))
      Vbimg[jvoxel] = mean(output(*,1)) 
      Kiimgsd[jvoxel] = stddev(output(*,0))
      Vbimgsd[jvoxel] = stddev(output(*,1))
      print,Kiimg[jvoxel], Vbimg[jvoxel], Kiimgsd[jvoxel], Vbimgsd[jvoxel]
    endif
    
    if fitmode eq 'vi' then begin
      ; do varaitional inference
    endif


  ; if jvoxel gt 20000  then stop
endfor


; print, 0.6*0.02/(1.4+0.02),max(Kiimg)      ; compare to the truth!
stop




; ===================================================================
; ================ projection NLLS fitting ==============
; ===================================================================
line3:
restore,filename='filetmpsino.sav'
help, tmpsino
nframe=158
inputfile = pad+'input.dat'
dataStruct = { plasma_t:0.0, plasma_c:0.0}
nrows = File_Lines(inputfile)
data = Replicate(dataStruct, nrows)
OpenR, lun, inputfile, /GET_LUN
ReadF, lun, data
Free_Lun, lun
plasma_t = (data.plasma_t)[0:nframe-1]   ; truncate equilibrium
plasma_c = (data.plasma_c)[0:nframe-1]
plasma_t = rebin(plasma_t,n_elements(plasma_t)/2)
plasma_c = rebin(plasma_c,n_elements(plasma_c)/2)
plasma_t0 = [0,plasma_t[0:n_elements(plasma_t)-2]]



; ; pre-compute the input function sinogram
; tmpsino1 = fltarr(172,252,837,nframes)
; for iframe = 0, nframes-1 do begin 
;   print, 'forward projecting frame ' + nistring(iframe)
;   img = fltarr(nrcols,nrrows,nrplanes) + plasma_c[iframe]
;   sino = 0
;   NIproj, img, sino, projd = projd0
;   tmpsino1[*,*,*,iframe] = sino
; endfor 
restore,filename='inputsino.sav'
stop

frameNr= n_elements(plasma_t0)
tstart = 70
tstop  = 120
output = double(fltarr(5)) 
debug  = 0
llsq_model = 0
isweight   = 0;
patlog     = '/home/tsun/bin/tpcclib-master/build/bin/libmtga_idl.so'
n1=(size(tmpsino))[1]
n2=(size(tmpsino))[2]
n3=(size(tmpsino))[3]
tmpsino = reform(tmpsino,long(n1)*n2*n3,nframes)
tmpsino1 = reform(tmpsino1,long(n1)*n2*n3,nframes)
Kisino  = fltarr(n1,n2,n3)
Vbsino  = Kisino*0
Kisinosd= Kisino*0
Vbsinosd= Kisino*0
SWSSsino= Kisino*0
for jvoxel = 0, long(n1)*n2*n3-1 do begin
    tac = tmpsino[jvoxel,*]
    ; plasma_cc = tmpsino1[jvoxel,*]
    plasma_cc = plasma_c
    if n_elements(where(tac ne 0)) le 1 then continue
    if (jvoxel mod 1000 eq 0) then print,'processed '+nistring(jvoxel)+' voxels...'
    ; do fitting
    success = call_external(patlog, 'patlak_idl', long(frameNr), double(plasma_t0), $
                double(plasma_t), double(reform(tac)), $
                double(plasma_cc),double(tstart),double(tstop),double(output),long(debug),long(llsq_model),$
                long(isweight)) 
    Kisino[jvoxel]   = output[0]
    Kisinosd[jvoxel] = output[1]
    Vbsino[jvoxel]   = output[2]
    Vbsinosd[jvoxel] = output[3]
    SWSSsino[jvoxel] = output[4]   
    ; if output[4] lt 0.9 then stop       ;just check the goodness of fitting 
  ; if jvoxel gt 20000  then stop
endfor

stop





; ===================================================================
; ================ projection Reconstucting (mean and variance) ==============
; ===================================================================

print, '-- Start Reconstruction --'
  time0 = systime(1)
  Kiimg   = 0
  Kiimgsd = 0
  Vbimg   = 0
  Vbimgsd = 0
  Kisinosd = Kisinosd*Kisinosd
  Vbsinosd = Vbsinosd*Vbsinosd
if test_back eq 1 then begin
  ; a test for back-projection
  NIproj, Kiimg, Kisino, projd = projd0, /back;, subset = 110
  NIproj, Kiimgsd, Kisinosd, projd = projd0, /back;, subset = 110
  NIproj, Vbimg, Vbsino, projd = projd0, /back;, subset = 110
  NIproj, Vbimgsd, Vbsinosd, projd = projd0, /back;, subset = 110
endif else begin
  ; start reconstruction
  nriter = 3     ; need more iterations than image-based approach?
  nrsub  = 32
  print, 'MLEM w/o motion paramereic images... '

  Kiimg = NImaposem(Kisino, projd = projd0, nriter=nriter, nrsub=nrsub,$
      /norm_per_sub,/shows) 
  Kiimgsd = NImaposem(Kisinosd, projd = projd0, nriter=nriter, nrsub=nrsub,$
      /norm_per_sub) 
  Vbimg = NImaposem(Vbsino, projd = projd0, nriter=nriter, nrsub=nrsub,$
      /norm_per_sub) 
  Vbimgsd = NImaposem(Vbsinosd, projd = projd0, nriter=nriter, nrsub=nrsub,$
      /norm_per_sub) 
endelse
print, 'Recon time (min): ', (systime(1)-time0)/60.

Kiimgsd = sqrt(Kiimgsd)
Vbimgsd = sqrt(Vbimgsd)

stop













if 1 then begin
; save the projections to fabber




; read parametric projections from fabber





; perform reconstructions for distributions



endif




if 0 then begin
; perfrom dynamic reconstructions





; save the reconstructed images to fabber




; read parametric images from fabber


endif










End