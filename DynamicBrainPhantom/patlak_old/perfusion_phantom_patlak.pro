

; Function getTac, imsize, icol, irow, iplane

;     ; read phantom images in one voxel
;     tac = []
;     for iframe = 0, 49 do begin
;         framename = nistring(iframe+1)   ;+'crop'
;         tmp = readraw(framename,imsize[0],imsize[1],imsize[2],'float')
;         framename2 = nistring(iframe+1)+'crop'
;         tmp2 = tmp[*,*,138:142]
;         writeraw, framename2, tmp2
;         tac = [tac, tmp[icol,irow,iplane]]
;     endfor
;     stop
;     return, tac

; End

Function getTac, imgs, imsize, icol, irow, iplane

    ; read phantom images in one voxel
    tac = imgs[icol,irow,iplane,*]
    ; tac = []
    ; for iframe = 0, 49 do begin
        ; framename = nistring(iframe+1)  +'crop'
        ; tmp = readraw(framename,imsize[0],imsize[1],imsize[2],'float')
        ; tac = [tac, tmp[icol,irow,iplane]]
    ; endfor
    return, tac
End

Pro modeling, tissue_t, tissue_c, plasma_t, plasma_c, parms, sigma, $
    plotcurve=plotcurve, startt=startt, stopt=stopt, xpatlak=xpatlak, ypatlak=ypatlak
  if n_elements(tissue_c) eq 0 then begin
    print, 'Please generate a tissue curve first'
    return
  endif

  modeldef = 'patlak'    ; subject to change, splii over

  ; Start params?
  ;---
  ;descrip = NImodeler(modeldef, 'parms')
  ;parms   = NImodeler(modeldef, 'start_parms')
  ;parms[i] = 1

  ; Ask which parameters to fix
  ;----------------
  ;descrip = NImodeler(modeldef, 'parms')
  ;fix_parms = intarr(n_elements(descrip))
  ;fix_parms[i] = 1


  ; Apply the model
  ;----------------
  if keyword_set(plotcurve) then begin
    niwin, 2
    wset, 2
    wshow, 2
  endif

  ; here we go for patlak directly
  if 1 then begin
    sigma = [-1,-1,-1]
   if n_elements(startt) ne 0 or n_elements(stopt) ne 0 then begin
       if n_elements(startt) eq 0 then startt = min(tissue_t)
       if n_elements(stopt)  eq 0 then stopt  = max(tissue_t)
     daar = where(tissue_t le stopt)
     if daar(0) eq -1 then stop
     tissue2_t = tissue_t(daar)
     tissue2_c = tissue_c(daar)
     daar = where(tissue_t ge startt)
     if daar(0) eq -1 then stop
     nrpoints = n_elements(daar)
     ;print, nrpoints
     
     parms = patlak(plasma_t, plasma_c, tissue2_t, tissue2_c, $
                      nrpoints = nrpoints, plotcurve=plotcurve, xpatlak=xpatlak, ypatlak=ypatlak)                     
   endif else begin   
     parms = patlak(plasma_t, plasma_c, tissue_t, tissue_c, $
                      nrpoints = nrpoints, plotcurve=plotcurve, xpatlak=xpatlak, ypatlak=ypatlak)
   endelse
  endif else begin
      parms = NImodeler(modeldef, plasma_t, plasma_c, tissue_t, tissue_c, $
            plotcurve=plotcurve, sigmaa = sigma, start_parms = parms, $
            weight=tissue_t / sqrt(tissue_c > max(tissue_c / 5)), $
            spill1_t = spill1_t, spill1_c = spill1_c, $
            spill2_t = spill2_t, spill2_c = spill2_c, $
            startt = startt, stopt = stopt, fix_parms = fix_parms)
            stop
  endelse

  descrip = NImodeler(modeldef, 'parms')
  ;for i = 0, n_elements(descrip)-1 do $
  ;   print, descrip(i), parms(i), sigma(i), $
  ;          format = "(A, ' = ', F8.4, ' +- ', F8.4)"

End


Pro perfusion_phantom_patlak

; get AIF curve
nframe = 50
tt = indgen(25)*2  ; 1:2:49
aifs = [0, 0, 0, 0, 25, 105, 220, 350, 440, 485, 430, 300, 180, 110, 104, 108, 115, 125, 115, 108, 98, 90, 98, 108, 112]
; ts = indgen(nframe)*1
ts = indgen(50)*1 ; 0:0.1:49
aif = interpol(aifs,tt,ts,/spline)

    nframe = 20
    tt = indgen(25)*2  ; 1:2:49
    aifs = [0, 0, 0, 0, 25, 105, 220, 350, 440, 485, 430, 300, 180, 110, 104, 108, 115, 125, 115, 108, 98, 90, 98, 108, 112]
    ; ts = indgen(nframe)*1
    ts = indgen(50)*1 ; 0:0.1:49
    aif = interpol(aifs,tt,ts,/spline)
    ts  = ts[0:nframe-1]
    aif = aif[0:nframe-1]

; load images
sz = [256,256,5]   ; [256,256,256]
imgs = fltarr(sz[0],sz[1],sz[2],nframe)
for iframe = 0, nframe-1 do begin
    framename = nistring(iframe+1)  +'crop'
    tmp = readraw(framename,sz[0],sz[1],sz[2],'float')
    tmp = niimgmu2hounsfield(tmp,/h2m) /10 ;*1e3        ; otherwise too small
    if iframe eq 0 then baseline = tmp
    imgs[*,*,*,iframe] = tmp - baseline          ; subtract baseline
endfor
          

; compute patlak for each voxel
startt = ts[10]
; stopt  = ts[40]
slopeimg     = fltarr(sz[0],sz[1],sz[2])
interceptimg = fltarr(sz[0],sz[1],sz[2])
; for k = 1, 1 do begin
for k = 0, sz[2]-1 do begin
print, 'processing plane number ' + nistring(k)
  for j = 0, sz[1]-1 do begin
    for i = 0, sz[0]-1 do begin
        tac = getTac(imgs, sz, i, j, k)
        

        ; stop

        ; patlak
        plasma_c = aif
        tissue_c = reform(tac)                     ; subject to change, e7imgs
        ; if total(tissue_c) lt 1e3 then continue                ;ignore background
        if tac[nframe-1] gt 0.0003*1 then begin          ; plot the curve of myocardium
            modeling, ts, tissue_c, ts, plasma_c, parms, sigma, startt=startt, stopt=stopt, /plotcurve 
            ; stop 
        endif else begin
            modeling, ts, tissue_c, ts, plasma_c, parms, sigma, startt=startt, stopt=stopt
        endelse
        slopeimg[i,j,k]     = parms[1]
        interceptimg[i,j,k] = parms[0]

    endfor
  endfor
endfor


stop

; validate the difference
; ct = ktrans*integral(aif)+CBV*aif
imgs_integral = imgs*0
for iframe=0,nframe-1 do begin
  if iframe eq 0 then imgs_integral[*,*,*,iframe] = slopeimg*float(aif[iframe]) + interceptimg*aif[iframe]
  if iframe gt 0 then imgs_integral[*,*,*,iframe] = slopeimg*int_tabulated(float(ts[0:iframe]),float(aif[0:iframe])) + interceptimg*aif[iframe]

;   imgs[*,*,*,iframe] += baseline
;   imgs_integral[*,*,*,iframe] += baseline
endfor

stop



; simulate projections







; reconstruct frames and save them







; derive AIF from images (optional)






; compute patlak for each voxel








End