;;;;;
; read data
pad = '/home/tsun/work/MK6240_datasets_for_Tao/'
patient = 'SMK011'   ; SMK005, SMK011
pad = pad + patient
doregister = 1

; read blood 
inputfile = pad+'/blood_data/plasma1_parent_dc_Bq_tD.txt'
dataStruct = { plasma_t:0.0,  plasma_c:0.0}
nrows = File_Lines(inputfile)
data = Replicate(dataStruct, nrows)
OpenR, lun, inputfile, /GET_LUN
ReadF, lun, data
Free_Lun, lun
plasma_t = (data.plasma_t)   ; truncate equilibrium
plasma_c = (data.plasma_c)

stop

; read dicoms
; dicomfolder_mr = 'Dicom_MPRAGE'
; dicomfolder_pet ='Dicom_PET_dyn_1'
pixpet = [1.1719,1.1719,2.8]
if patient eq 'SMK005' then pet = readraw(pad+'/mk.raw',256,256,4806,'int')
if patient eq 'SMK011' then pet = swap_endian(readraw(pad+'/mk.raw',256,256,4806,'int'))

petsize= size(pet)
nrcols = (size(pet))[1]
nrrows = (size(pet))[2]
nrplanes = 89
nrframe = 54     ;4806/nrplanes
petall = fltarr(petsize[1],petsize[2],nrplanes,nrframe)
for iframe = 0, nrframe-1 do petall[*,*,*,iframe]=pet[*,*,iframe*nrplanes:(iframe+1)*nrplanes-1]

if doregister then begin
    petall0 = petall

    ;register pets, use last frame as reference
    ref = petall[*,*,*,nrframe/2]
    for iframe = 0, nrframe-1 do begin
        if iframe eq nrframe/2 then continue
        if (iframe mod 20) eq 0 then print, 'Registering PET at frame ' + nistring(iframe)
        parms = 0
        tmp = petall[*,*,*,iframe]
        reg = niregisrigid(tmp,ref,method='mi',parms=parms, /norot, res=2)
        reg = niregisrigid(tmp,ref,method='mi',parms=parms, res=1)
        petall[*,*,*,iframe] = reg
    endfor

    ;register mr to pet
    pixmr = [1.0547,1.0547,1.2]
    if patient eq 'SMK005' then    mr = readraw(pad+'/mprage.raw',256,240,176,'int')   
    if patient eq 'SMK011' then    mr = swap_endian(readraw(pad+'/mprage.raw',240,256,176,'int'))
    if patient eq 'SMK005' then    mr = transpose(mr,[1,2,0])
    if patient eq 'SMK011' then    mr = transpose(mr,[2,0,1])   
    ; mr = congrid(mr, nrcols,nrrows,nrplanes)
    parms = 0
    mrreg = niregisrigid(mr,ref,method='mi',parms=parms, /scale, res=2, pix1=pixmr, pix2=pixpet)
    mrreg = niregisrigid(mr,ref,method='mi',parms=parms, /scale, res=1, pix1=pixmr, pix2=pixpet)
endif
stop

; reformat data 




plasma_t[0] = 0.01      ; prevent zero integration!!
petall1 = fltarr(nrcols,nrrows,nrplanes,nrframe/2)
for iframe = 0,nrframe-1,2 do $
    petall1[*,*,*,iframe/2] = petall[*,*,*,iframe]
nrframe /=2
petall = petall1

;;;;;
;  do the Logan analysis or VI or MCMC
plasma_c = congrid(plasma_c, nrframe)
plasma_t = congrid(plasma_t, nrframe)
plasma_t0 = [0,plasma_t[0:n_elements(plasma_t)-2]]

fitmode = 'nlls'
frameNr= n_elements(plasma_t0)
tstart = 40
tstop  = 90
output = double(fltarr(5)) 
debug  = 0
llsq_model = 0
isweight   = 0;
weights = fltarr(frameNr) + 1.0
logan_mode = 0            ; approximation with replacemnt of CR(t) as CP(t)
k2 = -1.0;    no reference tissue 
patlog     = '/home/tsun/bin/tpcclib-master/build/bin/libmtga_idl.so'
petall   = reform(petall,nrcols*nrrows*nrplanes,nrframe)
Kiimg  = fltarr(nrcols,nrrows,nrplanes)
Vbimg  = Kiimg*0
Kiimgsd= Kiimg*0
Vbimgsd= Kiimg*0
SWSSimg= Kiimg*0

for jvoxel = 0, long(nrcols)*nrrows*nrplanes-1 do begin
    tac = reform(petall[jvoxel, *])
    ; if  total(tac) lt 10000. then continue ;;n_elements(where(tac ne 0)) le 1 then continue,    n_elements(where(tac gt 0.)) le frameNr/2 or
    if (jvoxel mod 1000 eq 0) then print,'processed '+nistring(jvoxel)+' voxels...'

    if fitmode eq 'nlls' then begin
    ; do fitting
      success = call_external(patlog, 'logan_idl', long(frameNr), double(plasma_t0), $
                  double(plasma_t), double(tac), $
                  double(plasma_c),double(tstart),double(tstop),double(output),long(debug),long(llsq_model),double(k2),$
                  long(isweight),double(weights),long(logan_mode)) 
      Kiimg[jvoxel]   = output[0]
      Kiimgsd[jvoxel] = output[1]
      Vbimg[jvoxel]   = output[2]
      Vbimgsd[jvoxel] = output[3]
      SWSSimg[jvoxel] = output[4]
    endif
    ; stop
endfor


stop




















End