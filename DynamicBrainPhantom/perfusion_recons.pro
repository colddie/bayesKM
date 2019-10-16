
fdk      = 0
use_mltr = 0
use_mlem = 0

noise = 1S
nriter = 2                 ; niter for true_recon
nrsub = 40                 ; nrsub for all reconstruction
det_rebin = 1
angle_rebin = 1
raytracing = 0
   
maxmu_mm = 0.05
nrcols   = 256 
nrrows   = 256 
nrplanes = 5
nrcoldet = 400                ; 
nrrowdet = 400 /40
detcolsize = 0.2 *5   
detrowsize = 0.2 *5 
detdist = 216.5
focusheight = 575.
pixelsizemm = 1.0  ;/5
planesepmm  = 1.0 ;/5
coloffset   = 0.0
rowoffset   = -40.0 * 0      ; 0.0 if full scan
planeoffset = 0.0
imgcoloffset = 200.0 *0     ; offset scan

reconsize   = [1.0, 1.0, 1.0]  /5                                    ; add /2 if full scan                   
reconoffset = [coloffset, rowoffset+40, planeoffset] *0    

totangles = 180 *2
angles = findgen(totangles) * (2*!Pi) / totangles ;+ 40*!pi/180   ; 2* !pi


nframe=50
reconimgs = fltarr(nrcols,nrrows,nrplanes,nframe)
for iframe = 0, nframe-1 do begin
    print, 'processing frame... ' + nistring(iframe)
    framename = nistring(iframe+1)
    tmp = readraw(framename,nrcols,nrrows,256,'float')
    img = tmp[*,*,138:142]
    img = niimgmu2hounsfield(img, /h2m) /10


; simulate projections
    projd0  = nidef_conebeam(nrcols, nrrows, nrplanes, $
                            nrcoldet, nrrowdet, $
                            detcolsize  = detcolsize, $
					                  detrowsize  = detrowsize, $
                            detdist     = detdist, $
					                  focusheight = focusheight, $
                            angles       = angles, $
                            pixelsizemm  = pixelsizemm,  $
                            planesepmm   = planesepmm,   $
                            detcolcenter = detcolcenter, $  ;
                            detrowcenter = detrowcenter, $  ;
                            imgcoloffset = imgcoloffset, $  ;COR
                            coloffset    = coloffset, $
                            rowoffset    = rowoffset, $
                            planeoffset  = planeoffset, $
                            focuscoloffset = focuscoloffset, $
                            focusplaneoffset = focusplaneoffset, $
                            fwhm0 = fwhm0, raytracing = raytracing) 
    pix = [pixelsizemm, pixelsizemm, planesepmm]

  sinor = 0
  for i = 0, totangles-1 do begin    
    NIproj, img, sinor, subset = i, projd = projd0
  endfor

  ; --------------------------------- Add noise ----------------------------------------------------
  if noise eq 1 then begin
    blank = sinor * 0 + 5e6                                            ; less counts more noise, 1e6,1e5,1e4,1e3,5e2,3e2,1e2
    trans = blank * exp(-sinor)
    trans = nipoisson(seed, trans)
    sinor = -alog(trans/blank)  ; add back noise to sinogram
    if n_elements(where(~finite(sinor))) gt 1 then  sinor[where(~finite(sinor))] = 0.0
  endif else begin
    blank = sinor * 0 + 1e6
    trans = blank * exp(-sinor)
  endelse

  ; stop
   if noise eq 0 then begin
     save, filename='projs/projs' + nistring(iframe), sinor
   endif else begin
     save, filename='projs/projs_noise' + nistring(iframe), sinor
   endelse
; reconstruct frames and save them
;  projd_recon = projd0
    projd_recon  = nidef_conebeam(nrcols, nrrows, nrplanes, $
                            nrcoldet, nrrowdet, $
                            detcolsize  = detcolsize, $
					                  detrowsize  = detrowsize, $
                            detdist     = detdist, $
					                  focusheight = focusheight, $
                            angles       = angles, $
                            pixelsizemm  = pixelsizemm,  $
                            planesepmm   = planesepmm,   $
                            detcolcenter = detcolcenter, $  ;
                            detrowcenter = detrowcenter, $  ;
                            imgcoloffset = imgcoloffset, $  ;COR
                            coloffset    = coloffset, $
                            rowoffset    = rowoffset, $
                            planeoffset  = planeoffset, $
                            focuscoloffset = focuscoloffset, $
                            focusplaneoffset = focusplaneoffset, $
                            fwhm0 = fwhm0, raytracing = raytracing) 

 if fdk eq 1 then $
    niconebeamfdk, recon, sinor, projd = projd_recon, extrapolate = 100, /parker, /postscale

 
  if use_mltr eq 1 then $
     recon = NImapostr2(trans, blank, projd = projd_recon, nriter=nriter, $
      nrsub=nrsub, /shows)

  if use_mlem eq 1 then $
     recon = NImaposem(sinor, projd = projd_recon, nriter=nriter, $
        nrsub=nrsub,recon=recon, /norm_per_sub, /shows)	
    

  if fdk eq 1 or use_mltr eq 1 or use_mlem eq 1 then reconimgs[*,*,*,iframe] = recon


endfor
stop
   save, filename='recons' + nistring(iframe), reconimgs

End