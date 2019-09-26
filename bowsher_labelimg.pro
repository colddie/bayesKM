 
Pro bowsher_labelimg

 ; Computing threshold image for the Bowsher prior
  ;================================================
  
mr = 'phantomK1true_way.nii'
; restore, file= 'labelsimage.sav'            
labelsimage = niread_nii('labelsimg.nii',/no_reorient)
labelsimage = niread_nii(mr,/no_reorient)
bowsher_nrneighbors = 6; 6
bowsher_maxdiff = 0.5
; callexternal = 1;


     ; Check labelsimage
     if n_elements(labelsimage) le 1 then begin
        print, 'NIdef_prior, labelsimage has to be provided if bowsher_nrneighbors is set.'
      ;   return, 0
     endif
     
	   ; Default neighborhood
     if n_elements(neighborweights) le 1 then $
      ;   if (size(labelsimage))[0] eq 2 $
            ; neighborweights = [[0.,1.,0.], [1.,0,1.], [0.,1.,0.]] / 4. 
            neighborweights = [[[0,0,0],[0,1,0],[0,0,0]],$
                                [[0,1,0],[1,0,1],[0,1,0]],$
                                [[0,0,0],[0,1,0],[0,0,0]]] / 6.
  
     labelsthreshold = NIbowsher_th(labelsimage,                               $
                                    neighborweights     = neighborweights,     $ 
                                    bowsher_nrneighbors = bowsher_nrneighbors, $
                                    bowsher_maxdiff     = bowsher_maxdiff,     $
                                    callexternal        = NIrelease() ne 'Win32')
                                    ;callexternal        = !version.os ne 'Win32')

									
stop
   help, niread_nifti(mr, header=hdr)
   ; niwrite_nifti,'labelsthreshold.nii', labelsthreshold
   help, niwrite_nii(labelsthreshold, 'labelsthreshold.nii')   ;, nifti_hdr=hdr)


; combine the images latter into one file?

End
  
  