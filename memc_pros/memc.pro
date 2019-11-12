
; ----------------------
; Generated phantom with two-tissue irreversble compartmetn model
pad = '/home/tsun/bin/fsl/install/src/fabber_core/fabber_pet_c1/fdgC2rebin0.500000/'
imgs = niread_nii(pad+'actimg_fdgC2.nii',orientation='RAS') 
; imgs = niread_nii('actimg_fdgC2_0.00000_0.00500000_noise.nii',orientation='RAS') 
tmp = read_ascii(pad+'plasma_t.txt')
plasma_t = reform(tmp.FIELD1,n_elements(tmp.FIELD1),1) 
plasma_t = plasma_t[0:(size(imgs))[4]-1]
tmp = read_ascii(pad+'plasma_c.txt')
; tmp = read_ascii('ref_tissuec.txt')
plasma_c = reform(tmp.FIELD1,n_elements(tmp.FIELD1),1) 
plasma_c = plasma_c[0:(size(imgs))[4]-1]
plasma_c[0] = 0.
plasma_t[0] += 0.0001     ; prevent integrion over nothing
plasma_t0 = [0,plasma_t[0:n_elements(plasma_t)-2]]

; Patlak for phantom
nrcols = (size(imgs))[1]
nrrows = (size(imgs))[2]
nrplanes = (size(imgs))[3]
nframes= n_elements(plasma_t)
tstart = 20
tstop  = 60
output = double(fltarr(5)) 
debug  = 0
llsq_model = 0
isweight   = 0;
patlog     = '/home/tsun/bin/tpcclib-master/build/bin/libmtga_idl.so'
imgs   = reform(imgs,long(nrcols)*nrrows*nrplanes,nframes)
Kiimg  = fltarr(nrcols,nrrows,nrplanes)
Vbimg  = Kiimg*0
Kiimgsd= Kiimg*0
Vbimgsd= Kiimg*0
SWSSimg= Kiimg*0
for jvoxel = 0, long(nrcols)*nrrows*nrplanes-1 do begin
    tac = reform(imgs[jvoxel, *])
    if total(tac) lt 1e-3 then continue ;;n_elements(where(tac ne 0)) le 1 then continue
    if (jvoxel mod 1000 eq 0) then print,'processed '+nistring(jvoxel)+' voxels...'

    ; do fitting
      success = call_external(patlog, 'patlak_idl', long(nframes), double(plasma_t0), $
                  double(plasma_t), double(tac), $
                  double(plasma_c),double(tstart),double(tstop),double(output),long(debug),long(llsq_model),$
                  long(isweight)) 
      ; success = 
      Kiimg[jvoxel]   = output[0]
      Kiimgsd[jvoxel] = output[1]
      Vbimg[jvoxel]   = output[2]
      Vbimgsd[jvoxel] = output[3]
      SWSSimg[jvoxel] = output[4]

;   if jvoxel gt 20000  then stop
endfor

stop


; ----------------------
; Simulating motion in phantoms
; motion-free


; motion1 fast motion (in-plane only)
parms = [10/180.*!pi*0,0.0,0.0,2.0,5.0,0.0,1.,1.,1.,0,0,0]
imovframe = nframes/2
imgs = reform(imgs,nrcols,nrrows,nrplanes,nframes)
imgs1 = fltarr(nrcols,nrrows,nrplanes,nframes)
for iframe = 0, imovframe-1 do  $
    imgs1[*,*,*,iframe] = imgs[*,*,*,iframe]
for iframe = imovframe, nframes-1 do  begin
   tmp = nitransform(imgs[*,*,*,iframe],parms=parms,/interp)   
   imgs1[*,*,*,iframe] = tmp
endfor

; motion2 slow motion

; motion3 real rat motion (scaled)




;Patlak for moving phantom 
; Patlak for phantom
nframes= n_elements(plasma_t)
tstart = 20
tstop  = 60
output = double(fltarr(5)) 
debug  = 0L
llsq_model = 0
isweight   = 0;
patlog     = '/home/tsun/bin/tpcclib-master/build/bin/libmtga_idl.so'
imgs1   = reform(imgs1,long(nrcols)*nrrows*nrplanes,nframes)
Kiimg1  = fltarr(nrcols,nrrows,nrplanes)
Vbimg1  = Kiimg*0
Kiimgsd1= Kiimg*0
Vbimgsd1= Kiimg*0
SWSSimg1= Kiimg*0
for jvoxel = 0, long(nrcols)*nrrows*nrplanes-1 do begin
    tac = reform(imgs1[jvoxel, *])
    if total(tac) lt 1e-3 then continue ;;n_elements(where(tac ne 0)) le 1 then continue
    if (jvoxel mod 1000 eq 0) then print,'processed '+nistring(jvoxel)+' voxels...'

    ; do fitting
      success = call_external(patlog, 'patlak_idl', long(nframes), double(plasma_t0), $
                  double(plasma_t), double(tac), $
                  double(plasma_c),double(tstart),double(tstop),double(output),long(debug),long(llsq_model),$
                  long(isweight)) 
      ; success = 
      Kiimg1[jvoxel]   = output[0]
      Kiimgsd1[jvoxel] = output[1]
      Vbimg1[jvoxel]   = output[2]
      Vbimgsd1[jvoxel] = output[3]
      SWSSimg1[jvoxel] = output[4]

;   if jvoxel gt 20000  then stop
endfor
imgs1 = reform(imgs1,nrcols,nrrows,nrplanes,nframes)

help, niwrite_nii(imgs1,pad+'memc_movingPhantom.nii',orientation='RAS')
stop

; ----------------------
; Call external optimization program 
debug = 0L
tstart = float(tstart)
tstop = float(tstop)
nframe = (size(imgs1))[4]
imgfilename = pad+'memc_movingPhantom.nii'
lib = 'build/libmeKineticRigid.so'
model = 0   ; 0=patlak, 1=logan
fitmethod = 2L
plasma_tt = [[0], plasma_t[0:n_elements(plasma_t)-2]] 
parms0 = parms *0   ;;
; parms0[0] = parms[1]
; parms0[1] = parms[2]
; parms0[2] = parms[0]
; parms0[3] = parms[5]
; parms0[4] = parms[3]
; parms0[5] = parms[4]
nframeToFit = 1 ;;
rigmotion = fltarr(6,nframeToFit)
fitIndex = lonarr(nframe)
fitIndex[imovframe:nframe-1] = 1    ;;[0.0,...,1,1,...]
success = call_external(lib, 'meKineticRigid', nframe, imgfilename, parms0, tstart, tstop, $
                        double(plasma_tt),double(plasma_t), double(plasma_c), model, fitmethod,rigmotion, $
                         fitIndex, debug)





stop



End







