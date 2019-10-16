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


Pro perfusion_fitrecon

pad = 'projs/projs_noise'    ;'projs/projs'
tvv = 0
fit = 1
tv  = 0
test = 1
tips = 0
; tips_in = 1

; restore projections
nrcols    = 256 
nrrows    = 256 
nrplanes  = 5
baseline  = fltarr(nrcols,nrrows,nrplanes)
nrsub     = 20
nriter    = 2
sz        = [nrcols, nrrows, nrplanes]

det_rebin = 1
angle_rebin = 1
raytracing = 0
   
maxmu_mm = 0.05
nrcoldet = 400                ; 
nrrowdet = 400 /40
detcolsize = 0.2 *5   
detrowsize = 0.2 *5 
detdist    = 216.5
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

if test ne 1 then startf = 30 ;25
if test eq 1 then startf = 10

; get AIF curve
if test eq 0 then begin
    nframe = 50
    tt = indgen(25)*2  ; 1:2:49
    aifs = [0, 0, 0, 0, 25, 105, 220, 350, 440, 485, 430, 300, 180, 110, 104, 108, 115, 125, 115, 108, 98, 90, 98, 108, 112]
    ; ts = indgen(nframe)*1
    ts = indgen(50)*1 ; 0:0.1:49
    aif = interpol(aifs,tt,ts,/spline)

endif else begin
    nframe = 20
    tt = indgen(25)*2  ; 1:2:49
    aifs = [0, 0, 0, 0, 25, 105, 220, 350, 440, 485, 430, 300, 180, 110, 104, 108, 115, 125, 115, 108, 98, 90, 98, 108, 112]
    ; ts = indgen(nframe)*1
    ts = indgen(50)*1 ; 0:0.1:49
    aif = interpol(aifs,tt,ts,/spline)
    ts  = ts[0:nframe-1]
    aif = aif[0:nframe-1]
endelse

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

recon = 0
reconimgs = fltarr(nrcols,nrrows,nrplanes,nframe)
tt=systime(1)    
for iter = 1,nriter do begin
    subsetnum = 0
    for subiter = 0,nrsub - 1 do begin
        ;recon 
        print, 'Recon at iteration ' + nistring(iter) + ' subset ' + nistring(subiter)      
        subset = NIsino_subset(nrsub, totangles, subsetnum, subsetmask, $
        previoussubset)

        for iframe = 0, nframe-1 do begin
            reconimgs[*,*,*,iframe] += baseline
        endfor
        for iframe = 0, nframe-1 do begin
            restore, filename= pad + nistring(iframe)
            ssino = sinor[*, subset, *]  
            recon = reconimgs[*,*,*,iframe]
            if n_elements(where(recon gt 0)) eq 1 then recon += 1
            if tvv eq 1 then priord = NIdef_prior(/tvl1, markovw = 5.0e-1, markovparm = 0, $
                               nriter = tviters, opencl = 0, neighbor = nw)
            recon = NImaposem(ssino, projd=projd_recon, recon=recon, nriter=1, nrsub=1, presubset=subset, /norm_per_sub,priord=priod)      
            if iframe eq 0 then baseline = recon
            reconimgs[*,*,*,iframe] = recon - baseline
        endfor
        ttt = systime(1) - tt
        print, 'Time elapsed ' + nistring(ttt) + 's...'
        ; stop

    ;    if fit eq 1 and tips eq 1 then begin
    ;       for iplane = 0, sz[2]-1 do $
    ;       reconimgs[*,*,iplane,*] = transpose(bilateral(transpose(reform(reconimgs[*,*,iplane,*]),[2,0,1]),3,0.000016),[1,2,0])
    ;     ;   stop
    ;    endif

        ;fit patlak
        ; compute patlak for each voxel
        if fit eq 1 then begin
            print, 'Fitting... '  
            ; stopt  = ts[40]
            slopeimg     = fltarr(sz[0],sz[1],sz[2])
            interceptimg = fltarr(sz[0],sz[1],sz[2])

            for k = 0, sz[2]-1 do begin
                ; print, 'processing plane number ' + nistring(k)
                for j = 0, sz[1]-1 do begin
                    for i = 0, sz[0]-1 do begin
                        ; tac = getTac(reconimgs, sz, i, j, k)
                        tac = reconimgs[i,j,k,*]
                        
                        ;  stop

                        plasma_c = aif
                        tissue_c = reform(tac)                     ; subject to change, e7imgs
                        ; if total(tissue_c) lt 1e3 then continue                ;ignore background
                        if tac[nframe-1] gt 0.0003*1 then begin          ; plot the curve of myocardium
                            modeling, ts, tissue_c, ts, plasma_c, parms, sigma, startt=ts[startf], stopt=stopt, /plotcurve 
                            ; stop 
                        endif else begin
                            modeling, ts, tissue_c, ts, plasma_c, parms, sigma, startt=ts[startf], stopt=stopt
                        endelse
                        slopeimg[i,j,k]     = parms[1]
                        interceptimg[i,j,k] = parms[0]

                    endfor
                endfor
            endfor


            ttt = systime(1) - tt
            print, 'Time elapsed ' + nistring(ttt) + 's...'
            ; stop
        endif



        ; tv
        if fit eq 1 and tv eq 1 then begin
            print, 'Regulerization... '   
            tmp1 = slopeimg ;*(-1)  ; SHOULD REMOVE LATER
            tmp2 = interceptimg
            priord = NIdef_prior(/tvl1, markovw = 3.0e-3, markovparm = 0, $
                        nriter = tviters, opencl = 0, neighbor = nw)
            mapstep1 = tmp1 > (mean(tmp1)*1e-4)           ; from Nipost_prior
            mapstep2 = tmp2 > (mean(tmp2)*1e-4)           ;mean(tmp2)
            niprior_post, tmp1, priord, mapstep=mapstep1, /allow_neg     
            niprior_post, tmp2, priord, mapstep=mapstep2, /allow_neg 
            slopeimg = tmp1 ;*(-1)  ; SHOULD REMOVE LATER
            interceptimg = tmp2
            ttt = systime(1) - tt
            print, 'Time elapsed ' + nistring(ttt) + 's...'   
            ; stop
        endif



        if iter eq 20 and (subiter mod 19 eq 0) and fit eq 1 and tips eq 1 then begin     ; only after every iteration
        ;   for iplane = 0, sz[2]-1 do $
        ;   reconimgs[*,*,iplane,*] = transpose(bilateral(transpose(reform(reconimgs[*,*,iplane,*]),[2,0,1]),3,0.00012),[1,2,0])
          
          reconimgs_filter= fltarr(sz[0],sz[1],sz[2],nframe+2)
          reconimgs_tofilter= fltarr(sz[0],sz[1],sz[2],nframe+2)
          reconimgs_tofilter[*,*,*,0:nframe-1] = reconimgs
          reconimgs_tofilter[*,*,*,nframe] = slopeimg
          reconimgs_tofilter[*,*,*,nframe+1] = interceptimg
          for iplane = 0, sz[2]-1 do $
              reconimgs_filter[*,*,iplane,*] = transpose(bilateral(transpose(reform(reconimgs_tofilter[*,*,iplane,*]),[2,0,1]),3,0.00012),[1,2,0])
        ;   reconimgs =  reconimgs_filter[*,*,*,0:nframe-1]   ; not useful?
          slopeimg = reconimgs_filter[*,*,*,nframe] 
          interceptimg = reconimgs_filter[*,*,*,nframe+1]
          stop
        endif



        ; update current image 
        if fit eq 1 then begin
            ; imgs_integral = imgs*0
            print, 'Updating... '  
            backup = reconimgs
            for iframe=startf,nframe-1 do begin   ; was 0,nframe-1
            ; if iframe eq startf then reconimgs[*,*,*,iframe] = slopeimg*float(aif[iframe]) + interceptimg*aif[iframe]
            ; if iframe gt startf then reconimgs[*,*,*,iframe] = slopeimg*int_tabulated(float(ts[startf:iframe]),float(aif[startf:iframe])) + interceptimg*aif[iframe]
            if iframe eq startf then reconimgs[*,*,*,iframe] = slopeimg*float(aif[iframe]) + interceptimg*aif[iframe]
            if iframe gt startf then reconimgs[*,*,*,iframe] = slopeimg*int_tabulated(float(ts[startf:iframe]),float(aif[startf:iframe])) + interceptimg*aif[iframe]
            ; reconimgs[*,*,*,iframe] *= (1-0.75)/(1-0.55)      ; hematocrit correction factor
            endfor
            ttt = systime(1) - tt
            print, 'Time elapsed ' + nistring(ttt) + 's...'   
            stop
        endif
  endfor
endfor

    ; post fitting...
    if fit eq 0 then begin
    ;    stop
    ;    if tips eq 1 then begin
    ;       for iplane = 0, sz[2]-1 do $
    ;       reconimgs[*,*,iplane,*] = transpose(bilateral(transpose(reform(reconimgs[*,*,iplane,*]),[2,0,1]),3,0.00012),[1,2,0])
    ;     ;   stop
    ;    endif
      
            print, 'Fitting... '  
            ; if test ne 1 then startf = 30 ; 25
            ; if test eq 1 then startf = 10
            ; stopt  = ts[40]
            slopeimg     = fltarr(sz[0],sz[1],sz[2])
            interceptimg = fltarr(sz[0],sz[1],sz[2])

            for k = 0, sz[2]-1 do begin
                ; print, 'processing plane number ' + nistring(k)
                for j = 0, sz[1]-1 do begin
                    for i = 0, sz[0]-1 do begin
                        ; tac = getTac(reconimgs, sz, i, j, k)
                        tac = reconimgs[i,j,k,*]
                        
                        ; stop

                        plasma_c = aif
                        tissue_c = reform(tac)                     ; subject to change, e7imgs
                        ; if total(tissue_c) lt 1e3 then continue                ;ignore background
                        if tac[nframe-1] gt 0.0003*1 then begin          ; plot the curve of myocardium
                            modeling, ts, tissue_c, ts, plasma_c, parms, sigma, startt=ts[startf], stopt=stopt, /plotcurve 
                            ; stop 
                        endif else begin
                            modeling, ts, tissue_c, ts, plasma_c, parms, sigma, startt=ts[startf], stopt=stopt
                        endelse
                        slopeimg[i,j,k]     = parms[1]
                        interceptimg[i,j,k] = parms[0]

                    endfor
                endfor
            endfor


            ; imgs_integral = imgs*0
            print, 'Updating... '  
            backup = reconimgs
            for iframe=startf,nframe-1 do begin   ; was 0,nframe-1
            ; if iframe eq startf then reconimgs[*,*,*,iframe] = slopeimg*float(aif[iframe]) + interceptimg*aif[iframe]
            ; if iframe gt startf then reconimgs[*,*,*,iframe] = slopeimg*int_tabulated(float(ts[startf:iframe]),float(aif[startf:iframe])) + interceptimg*aif[iframe]
            if iframe eq startf then reconimgs[*,*,*,iframe] = slopeimg*float(aif[iframe]) + interceptimg*aif[iframe]
            if iframe gt startf then reconimgs[*,*,*,iframe] = slopeimg*int_tabulated(float(ts[startf:iframe]),float(aif[startf:iframe])) + interceptimg*aif[iframe]
            reconimgs[*,*,*,iframe] *= (1-0.55)/(1-0.75)     ; hematocrit correction factor

            endfor
            ttt = systime(1) - tt
            print, 'Time elapsed ' + nistring(ttt) + 's...'   
            stop


        if tips eq 1 then begin
        ;   for iplane = 0, sz[2]-1 do $
        ;   reconimgs[*,*,iplane,*] = transpose(bilateral(transpose(reform(reconimgs[*,*,iplane,*]),[2,0,1]),3,0.00012),[1,2,0])
          
          reconimgs_filter= fltarr(sz[0],sz[1],sz[2],nframe+2)
          reconimgs_tofilter= fltarr(sz[0],sz[1],sz[2],nframe+2)
          reconimgs_tofilter[*,*,*,0:nframe-1] = reconimgs
          reconimgs_tofilter[*,*,*,nframe] = slopeimg
          reconimgs_tofilter[*,*,*,nframe+1] = interceptimg
          for iplane = 0, sz[2]-1 do $
              reconimgs_filter[*,*,iplane,*] = transpose(bilateral(transpose(reform(reconimgs_tofilter[*,*,iplane,*]),[2,0,1]),3,0.00012),[1,2,0])
        endif

            ttt = systime(1) - tt
            print, 'Time elapsed ' + nistring(ttt) + 's...'
            ; stop
    endif


; alternate updates
  
  ;recon


  ;fit



  ;regulerization



  ;update

stop
End

  