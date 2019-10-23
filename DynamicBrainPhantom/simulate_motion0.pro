    ; restore,filename = 'tmp.sav'
    restore,filename = 'tmp_delay.sav'
    nframe = (size(tacall))[4]
    tacall_transform = tacall*0


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




    ;;;;;;;
    ; simple motion applied on each frame
    if 1 then begin
        for iframe = 0,nframe-1 do begin
        if iframe mod 10 eq 0 then  print, 'Process frame' + nistring(iframe)
            curimg = tacall[*,*,*,iframe]
            motion = [rigmotion[*,iframe],1,1,1,0,0,0]
            img_transform = nitransform(curimg,parms=motion,/interp)
            tacall_transform[*,*,*,iframe] = img_transform
            ; stop
        endfor

        tacall = tacall_transform
        save, filename = 'tmp_delay_transform0.sav', tacall
        
        for iframe = 0,nframe-1 do $
           tacall_transform[*,*,*,iframe] -= baselineall

        stop
    endif 


    ;;;;;;;
    ;  motion applied on projections
    if 0 then begin
    
        stop
    endif



 cgplot,ts,aif/10,yrange=[-10,250]     
 cgplot,congrid(ts/10,50),tacall_transform[155,114,154,*],col=1,yrange=[-10,250],/noerase
 cgplot,congrid(ts/10,50),tacall_transform[104,151,154,*],col=2,yrange=[-10,250],/noerase
 cgplot,congrid(ts/10,50),tacall_transform[146,99,154,*],col=4,yrange=[-10,250],/noerase 
 cgplot,congrid(ts/10,50),tacall_transform[179,119,154,*],col=5,yrange=[-10,250],/noerase
 cgLegend, color=[0,1,2,4,5],Symsize=1.5, Location=[0.625, 0.9], $
      Titles=['AIF','infact','severe-infact','grey','white'], Length=0.075, /Box, VSpace=2.75


End