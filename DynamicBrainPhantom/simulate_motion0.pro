    restore,filename = 'tmp.sav'
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
  ;rigmotion[5,200:400] = 0   ; one motion one time




    ;;;;;;;
    ; simple motion applied on each frame
    if 1 then begin
        for iframe = 0,nframe-1 do begin
            curimg = tacall[*,*,*,iframe]
            motion = [rigmotion[*,iframe],1,1,1,0,0,0]
            img_transform = nitransform(curimg,parms=motion,/interp)
            tacall_transform[*,*,*,iframe] = img_transform
            ; stop
        endfor

        tacall = tacall_transform
        save, filename = 'tmp_transform0.sav', tacall
        
        for iframe = 0,nframe-1 do $
           tacall_transform[*,*,*,iframe] -= baselineall

        stop
    endif 


    ;;;;;;;
    ;  motion applied on projections
    if 0 then begin
    
        stop
    endif



End