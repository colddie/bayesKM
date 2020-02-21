    ; restore,filename = 'tmp.sav'
    ; restore,filename = 'tmp_delay.sav'
    ; restore,filename = 'tmp_new.sav'
    restore,filename = 'tmp_delay_new.sav'
    nframe = (size(tacall))[4]
    tacall_transform = tacall*0

    simul_motion = 0
    add_noise    = 0
    use_fdk      = 1
    use_mltr     = 0
    raytracing   = 0



    ;;;;;;;
    ; simple motion applied on each frame
    if simul_motion then begin
        posefile = 'rat_pose_data_5s.dat'
        print, 'Reading motion from ' + posefile
        rigmotion = motion_from_posefile(posefile, $
        totangles = nframe, $
        shiftfactor = 3, $
        repeated = repeated)
        rigmotion[1:2,*] = 0
        rigmotion[5,*] = 0
        rigmotion /= 5.
        ;rigmotion[5,200:400] = 0   ; one motion one time

        for iframe = 0,nframe-1 do begin
        if iframe mod 10 eq 0 then  print, 'Process frame' + nistring(iframe)
            curimg = tacall[*,*,*,iframe]
            motion = [rigmotion[*,iframe],1,1,1,0,0,0]
            img_transform = nitransform(curimg,parms=motion,/interp)
            tacall_transform[*,*,*,iframe] = img_transform
        endfor

        ; tacall = tacall_transform
        ; save, filename = 'tmp_delay_transform0.sav', tacall
        
        ; for iframe = 0,nframe-1 do $
        ;    tacall_transform[*,*,*,iframe] -= baselineall

    endif else tacall_transform = tacall


    ;;;;;;;
    ; TAC plot
    cgplot,congrid(ts,50),congrid(aif,nframe),yrange=[-10,500]     
    cgplot,congrid(ts,50),tacall[103,119,151,*],col=1,yrange=[-10,500],/noerase
    cgplot,congrid(ts,50),tacall[105,84,151,*],col=4,yrange=[-10,500],/noerase 
    cgplot,congrid(ts,50),tacall[104,51,151,*],col=5,yrange=[-10,500],/noerase
    cgLegend, color=[0,1,4,5],Symsize=1.5, Location=[0.625, 0.9], $
          Titles=['AIF','infact','grey','white'], Length=0.075, /Box, VSpace=2.75
    ; cgplot,congrid(ts/10,50),aif/10,yrange=[-10,250]     
    ; cgplot,congrid(ts/10,50),tacall_transform[155,114,154,*],col=1,yrange=[-10,250],/noerase
    ; cgplot,congrid(ts/10,50),tacall_transform[104,151,154,*],col=2,yrange=[-10,250],/noerase
    ; cgplot,congrid(ts/10,50),tacall_transform[146,99,154,*],col=4,yrange=[-10,250],/noerase 
    ; cgplot,congrid(ts/10,50),tacall_transform[179,119,154,*],col=5,yrange=[-10,250],/noerase
    ; cgLegend, color=[0,1,2,4,5],Symsize=1.5, Location=[0.625, 0.9], $
    ;       Titles=['AIF','infact','severe-infact','grey','white'], Length=0.075, /Box, VSpace=2.75



stop


    ;;;;;;;
    ;  simulate projections
    ncols   = (size(tacall))[1]
    nrows   = (size(tacall))[2]
    nplanes = (size(tacall))[3]
    nframes = (size(tacall))[4] /1.

    slicenumber = 32                          ;40  ; 16
    flyingfocus         = 0
    angles_per_rotation = 150*4 ;150 ;250 ;500 ;1160?
    pitch      = 1.0 ;0.8 ;1.0; 2.0                        ; 0.5 oversample projection error will not converge

    pixelsizemm = 1. /2
    planesepmm = 1.                                       ; 1.5, 3 mm
    slicesepmm = 1.                                         ;1.  ; 1.5 mm
    tablefeed  = pitch * slicenumber * slicesepmm

    nrrotations = nplanes * planesepmm / (tablefeed) + 1.0

    totangles = fix(angles_per_rotation * nrrotations)
    ;print, 'totangles = ' + NIstring(totangles)

    tubeangle = -findgen(totangles) * (2.0 * !pi) / angles_per_rotation

    tablepos = fltarr(totangles, slicenumber)
    tpos     = -findgen(slicenumber) * slicesepmm
    tpos    -= mean(tpos)
    relpos   = -findgen(totangles) * tablefeed / angles_per_rotation
    relpos  -= mean(relpos)
    for i = 0, slicenumber-1 do $
      tablepos[*,i] = -(relpos + tpos[i])
    rowoffset = 0.0
    coloffset = 0.0
    alignment = tubeangle * 0.0


    ;------
    projd = NIdef_projdistdspiralct(tubeangle, tablepos, alignment, $
      ncols, nrows, nplanes, $
      pixelsizemm, planesepmm, $
      zalignment = zalignment, $
      det_rebin = det_rebin, $
      angle_rebin = angle_rebin, $
      ctmodel = ctmodel, $
      rowoffset = rowoffset, $
      coloffset = coloffset, $
      raytracing = raytracing)
      
    pix = [pixelsizemm, pixelsizemm, planesepmm]

    sinoall = fltarr(projd.ndetcols,totangles,projd.ndetplanes,nframes) 
    blankall = sinoall *0
    transall = sinoall *0
    for iframe = 0, nframes-1 do begin
      print, 'Projecting frame '+nistring(iframe)+'...'
       sinor = 0
       img = tacall_transform[*,*,*,iframe]
       img = NIimgmu2hounsfield(img,/h2m)/10.
       for i = 0, totangles-1 do $    
          NIproj, img, sinor, subset = i, projd = projd

       ;;;;;;;
       ; Simulate noise
      if add_noise eq 1 then begin
          print, 'Adding noise to frame '+nistring(iframe)+'...'
          blank = sinor * 0 + 1e5                                           ; less counts more noise, 1e6,1e5,1e4,1e3,5e2,3e2,1e2
          trans = blank * exp(-sinor)
          trans = nipoisson(seed, trans)
          sinor = -alog(trans/blank)  ; add back noise to sinogram
          if n_elements(where(~finite(sinor))) gt 1 then  sinor[where(~finite(sinor))] = 0.0
      endif else begin
          blank = sinor * 0 + 1e6
          trans = blank * exp(-sinor)
      endelse
      
      sinoall[*,*,*,iframe]  = sinor
      blankall[*,*,*,iframe] = blank 
      transall[*,*,*,iframe] = trans
    endfor
    ; stop




    ;;;;;;;
    ; Reconstructimg dynamic frames
    nriter = 4
    nrsub  = 40
    reconall = fltarr(ncols,nrows,nplanes,nframes)
    for iframe = 0, nframes-1 do begin
      print, 'Reconstructing frame '+nistring(iframe)+'...'
      recon = 0
      sinor = sinoall[*,*,*,iframe]


      if use_fdk eq 1 then begin
          print, 'FDK...'
          recon =  nispiralfdk(projd, sinor, $
                 angles_per_rotation=angles_per_rotation, /postscale) > 0.
      endif else if use_mltr eq 1 then begin
          print, 'MLTR w/o motion correction...'
          blank = float(sinor gt 0) * 1e5
          ;blank = sinor * 0 + 1e6
          trans = blank * exp(-sinor)
          recon = NImapostr2(trans, blank, projd = projd, nriter=nriter, $
            nrsub=nrsub, /shows)
      endif else begin
          print, 'MLEM w/o motion correction...'
          recon = NImaposem(sinor, projd = projd, nriter=nriter, $
            nrsub=nrsub, /norm_per_sub, /shows)
     endelse
      
      recon = NIimgmu2hounsfield(recon*10, /m2h)
      reconall[*,*,*,iframe] = recon
    endfor


    if add_noise eq 0 then  save, filename='simulatedphantom.sav', reconall
       stop 
    if simul_motion eq 0 then save, filename='simulatedphantom_motionfree.sav', reconall
    if simul_motion eq 1 then save, filename='simulatedphantom_moving.sav', reconall 



End
