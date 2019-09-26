fit_mode = 1;   ;0 - fabber, 1 - regular nlls fitting


t0 = systime(1)

; read dynamic FDG phantom generated from phantomdynamic_c1c2.pro
pad = '/home/tsun/bin/fsl/install/src/fabber_core/fabber_pet_c1/'
; tmpimg = niread_nii(pad+'actimg_fdgC2_0.00000_0.00500000_noise.nii')
; tmpimg = niread_nii(pad+'actimg_fdgC2.nii')
restore, filename = pad+'actimg_fdgC2.sav'

; stop


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
; goto, line3     ; skip the rest and go directly to fitting
tmpsino = fltarr(172,252,837,nframes)

; if usempi eq 0 then begin 
  ; for iframe = 0, nframes-1 do begin 
  ;   print, 'forward projecting frame ' + nistring(iframe)
  ;   img = congrid(tmpimg[*,*,*,iframe], nrcols,nrrows,nrplanes)
  ;   sino = 0
  ;   NIproj, img, sino, projd = projd0
  ;   if addnoise then $                 ;add Poisson noise
  ;     sino  = nipoisson(seed,sino)  
  ;   tmpsino[*,*,*,iframe] = sino
  ;   ; stop
  ; endfor 
; endif else begin
  mpi_rank, rank
  start = 0
  stop  = nframes
  ;; divide task among MPI ranks
  mpi_helper, start, stop, i_start, i_stop

  ;; main loop
  for iframe = i_start, i_stop-1 do begin $ 
    print, 'forward projecting frame_' + nistring(iframe)+'_at rank' +nistring(rank) &$
    img = congrid(tmpimg[*,*,*,iframe], nrcols,nrrows,nrplanes) &$
    sino = 0 &$
    NIproj, img, sino, projd = projd0  &$
    if addnoise then sino  = nipoisson(seed,sino)  &$                 ;add Poisson noise
    tmpsino[*,*,*,iframe] = sino
    ; stop
;   endfor 
  save, filename='sinogram trunck'+'_at rank' +nistring(rank), rank
; endelse

print, 'TIme elapsed '+ nistring(systime(1)-t0)+' seconds...'
exit

; End