; % creates 4D perfusion phantom brain volume
; % adapted from
; % author: michael manhart | michael.manhart@cs.fau.de
; %         pattern recognition lab, university of erlangen-nuremberg


;.. below function is to be in c++
; Function GetTissueTac, aif, t, dt, mtt, cbf
;     cbf         = cbf/6000.              ;% ml/100ml/min = 1/(100*60s)
;     tmin        = mtt 
;     r           = cbf *exp( -(t-mtt)) 
;     r[where(t<tmin)]   = cbf                 ;%   ...residue function
;     ; tac         = dt * convol(aif,r)         ; % indicator-dilution theory
;     tac = aif *0
;     lib = '/home/tsun/bin/tpcclib-master/build/bin/libmtga_idl.so'
;     success = call_external(lib,'convolut_idl',double(aif),long(n_elements(aif)),double(r),long(n_elements(r)),double(tac))
;     tac         = tac[0:n_elements(aif)-1]      ;   % trf should have same size as aif
;     stop
;     return, tac
; End



Pro create_phantom_pct

dt       = 1    ; %s
dtt      = 0.1
addskull = 1
ncol     = 256
nrow     = 256 
nplane   = 256
tstart   = 0 ; %s
tstop    = 49 ; %s
delay    = 0

; % define perfusion parameters for different matters
; % gm -> gray matter wm -> white matter
; % gm_cbv -> healthy gray matter blood flow
; % gmr_cbv -> reduced gray matter blood flow
; % gmsr_cbv -> severely reduced gray matter blood flow etc...
; % cbv in ml/100ml  cbf in ml/100 ml/min  mtt in seconds
gm_cbv = 3.3
gm_cbf = 53; 
gm_mtt = (gm_cbv/gm_cbf)*60;
gm_delay = 1.0;   
gmr_cbv = 3;
gmr_cbf = 16;
gmr_mtt = (gmr_cbv/gmr_cbf)*60; 
gmr_delay = 1.0*4;    
gmsr_cbv = 0.71;
gmsr_cbf = 5.3;
gmsr_mtt = (gmsr_cbv/gmsr_cbf)*60;   
gmsr_delay = 3.0*3;   

wm_cbv = 1.9;
wm_cbf = 25;
wm_mtt = (wm_cbv/wm_cbf)*60;  
wm_delay = 1.0;    
wmr_cbv = 1.7;   
wmr_cbf = 7.5; 
wmr_mtt = (wmr_cbv/wmr_cbf)*60;    
wmr_delay = 2.0*4;   
wmsr_cbv = 0.42;
wmsr_cbf = 2.5;
wmsr_mtt = (wmsr_cbv/wmsr_cbf)*60;  
wmsr_delay = 4.0*3;   

; % define perfusion parameter variation w.r.t. mr data
cbf_var = 14;
mtt_var = 0.7;  
cbfr_var = 4.25;
mttr_var = 0.75;
cbfsr_var = 1.4;
mttsr_var = 1;   
delay_var = 0.0
delayr_var = 0.5
delaysr_var = 1.0     

; % tissue attenuation offsets
HU_GM = 35
HU_WM = 29
HU_V = 50
HU_CSF = 12;
HU_GM_VAR = 3
HU_WM_VAR = 3
HU_CSF_VAR = 3;
; % classes of tissue
GM = 1
WM = 2
GMR = 3          ; reduced
GMSR = 4        ; severely reduced
WMR = 5
WMSR = 6
AI = 7
VO = 8
CSF = 9
HEGM = 11     ;;;??? matlab generate phantom with wrong index!
HEWM = 10

; % load brain maps
; brain.mat not found! create stroke annotation with strokecreator before running create_phantom!
; brain = readraw('brain.raw',ncol,nrow,nplane,'byte')
load_mat,'brain.mat', store_level=-1
; restore,filename='brain.sav'
brain = byte(brain)
brain[where(brain eq HEWM)] = WM;
brain[where(brain eq HEGM)] = GM;
; % load brain mr data for perfusion parameter variation
mrbrain = readraw('mrbrain.raw',ncol,nrow,nplane,'uint')
; % load skull data
skull = readraw('skull.raw',ncol,nrow,nplane,'int')

; % check if phantom files to create already exist
; % create aif 0..49 s in 0.1 sec sampling interval
nframe = (tstop-tstart)/dt+1
t = indgen(nframe)*dt;     ;0:1:49
to = indgen(25)*2 +1  ;1:2:49;
aifs = [0, 0, 0, 0, 25, 105, 220, 350, 440, 485, 430, 300, 180, 110, 104, 108, 115, 125, 115, 108, 98, 90, 98, 108, 112];
ts = indgen(tstop*10+1)*dtt   ;0:0.1:49;
aif = interpol(float(aifs),to,ts, /spline);

stop

; % outputs
baselineall = brain *0
tacall = fltarr(ncol,nrow,nplane,nframe)
cbfall = float(brain) *0
cbvall = float(brain) *0
mttall = float(brain) *0
delayall = float(brain) *0
ttpall = float(brain) *0


; % Iterative over each planes to derive TACs
; size = size(brain)
for z=0,nplane-1 do begin
  print, 'Processing slice '+nistring(z)
  cbf =0
  mtt =0
  cbv =0
  ttp =0
  delay =0
  baselineslice = fltarr(ncol,nrow)
  tacslice = fltarr(ncol,nrow,nframe)
  cbfslice = baselineslice *0
  cbvslice = baselineslice *0
  mttslice = baselineslice *0
  delayslice = baselineslice *0
  ttpslice = baselineslice *0
  mrslice = mrbrain[*,*,z]
  brainslice = brain[*,*,z]

  ; % Process white matter
;    wm_idx = where(brainslice eq WM || brainslice eq WMR || brainslice eq WMSR)
;    wm_idx = SetUnion(SetUnion(where(brainslice eq WM),where(brainslice eq WMR)),where(brainslice eq WMSR))
    wm_idx = [where(brainslice eq WM),where(brainslice eq WMR),where(brainslice eq WMSR)]
    wm_idx = [wm_idx[where(wm_idx ne -1)]]
   if n_elements(wm_idx) ne 1 then begin
     mrvalues = mrslice[wm_idx]
     mrvalues = mrvalues - mean(mrvalues)
     mrvalues_bound = 2*sqrt(variance(mrvalues));
     mrvalues[where(mrvalues lt -mrvalues_bound)] = -mrvalues_bound;
     mrvalues[where(mrvalues gt  mrvalues_bound)] =  mrvalues_bound;
     mrvalues = mrvalues/mrvalues_bound;
     baselineslice[wm_idx] = HU_WM - (mrvalues*HU_WM_VAR);

    for i=0,n_elements(wm_idx)-1 do begin
      idx = wm_idx[i]
      case brainslice[idx] of
        WM: begin
            cbf = wm_cbf - mrvalues[i]*cbf_var;
            mtt = wm_mtt + mrvalues[i]*mtt_var;
            delay = wm_delay + mrvalues[i]*delay_var
        end
        WMR: begin
            cbf = wmr_cbf - mrvalues[i]*cbfr_var;
            mtt = wmr_mtt + mrvalues[i]*mttr_var;
            delay = wmr_delay + mrvalues[i]*delayr_var           
        end
        WMSR: begin
            cbf = wmsr_cbf - mrvalues[i]*cbfsr_var;
            mtt = wmsr_mtt + mrvalues[i]*mttsr_var; 
            delay = wmsr_delay + mrvalues[i]*delaysr_var
        end
      endcase

        ; tac = GetTissueTac(aif,ts,0.1,mtt,cbf);
        tac = double(fltarr(n_elements(ts)))
        lib = '/home/tsun/bin/tpcclib-master/build/bin/libmtga_idl.so'
        success = call_external(lib,'simpct_idl',double(ts),double(aif),long(n_elements(ts)),double(cbf),double(mtt),double(delay),tac)
        for ti = 0,nframe-1 do begin
            tac_index = round(t[ti]/0.1) ;+ 1;
            sliceposition = array_indices(baselineslice,idx);
            ; [slice_x slice_y] = ind2sub(size(baselineslice),idx);
            tacslice[sliceposition[0],sliceposition[1],ti] = tac[tac_index]+baselineslice[idx];
            ; stop
        endfor   
        cbfslice[idx] = cbf;
        cbvslice[idx] = cbf*mtt/60;
        mttslice[idx] = mtt;
        delayslice[idx] = delay;
        ttp = where(tac eq max(tac));
        ; [~, ttp] = max(tac); 
        ttpslice[idx] = ttp;
    endfor
   endif

  ; % Process grey matter
;   gm_idx = where(brainslice eq GM || brainslice eq GMR || brainslice eq GMSR)
;   gm_idx = SetUnion(SetUnion(where(brainslice eq GM),where(brainslice eq GMR)),where(brainslice eq GMSR))
    gm_idx = [where(brainslice eq GM),where(brainslice eq GMR),where(brainslice eq GMSR)]
    gm_idx = [gm_idx[where(gm_idx ne -1)]]
  if n_elements(gm_idx) ne 1 then begin
     mrvalues = mrslice[gm_idx]
     mrvalues = mrvalues - mean(mrvalues)
     mrvalues_bound = 2*sqrt(variance(mrvalues));
     mrvalues[where(mrvalues lt -mrvalues_bound)] = -mrvalues_bound;
     mrvalues[where(mrvalues gt  mrvalues_bound)] =  mrvalues_bound;
     mrvalues = mrvalues/mrvalues_bound;
     baselineslice[gm_idx] = HU_GM - (mrvalues*HU_GM_VAR);

    for i=0,n_elements(gm_idx)-1 do begin
      idx = gm_idx[i]
      case brainslice[idx] of
        GM: begin
            cbf = gm_cbf - mrvalues[i]*cbf_var;
            mtt = gm_mtt + mrvalues[i]*mtt_var;
            delay = gm_delay + mrvalues[i]*delay_var
        end
        GMR: begin
            cbf = gmr_cbf - mrvalues[i]*cbfr_var;
            mtt = gmr_mtt + mrvalues[i]*mttr_var;
            delay = gmr_delay + mrvalues[i]*delayr_var
        end
        GMSR: begin
            cbf = gmsr_cbf - mrvalues[i]*cbfsr_var;
            mtt = gmsr_mtt + mrvalues[i]*mttsr_var; 
            delay = wmsr_delay + mrvalues[i]*delaysr_var
        end
      endcase

        ; tac = GetTissueTac(aif,ts,0.1,mtt,cbf);
        tac = double(fltarr(n_elements(ts)))
        lib = '/home/tsun/bin/tpcclib-master/build/bin/libmtga_idl.so'
        success = call_external(lib,'simpct_idl',double(ts),double(aif),long(n_elements(ts)),double(cbf),double(mtt),double(delay),tac)
        for ti = 0,nframe-1 do begin
            tac_index = round(t[ti]/0.1) ;+ 1;
            sliceposition = array_indices(baselineslice,idx);
            ; [slice_x slice_y] = ind2sub(size(baselineslice),idx);
            tacslice[sliceposition[0],sliceposition[1],ti] = tac[tac_index]+baselineslice[idx];
        endfor   
        cbfslice[idx] = cbf;
        cbvslice[idx] = cbf*mtt/60;
        mttslice[idx] = mtt;
        delayslice[idx] = delay;
        ttp = where(tac eq max(tac));
        ; [~, ttp] = max(tac); 
        ttpslice[idx] = ttp*0.1;

        ; if cbf gt 0 then stop
    endfor
   endif


  ; % Process arteries
  ai_idx = where(brainslice eq AI)
  if  n_elements(ai_idx) gt 1 then begin
    baselineslice[ai_idx] = HU_V;
    for i=0,n_elements(ai_idx)-1 do begin
        idx = ai_idx[i]
        for ti=0,nframe-1  do begin
        tac_index = round(t[ti]/0.1) ;+ 1;
        sliceposition = array_indices(baselineslice,idx);
        ; [slice_x slice_y] = ind2sub(size(baselineslice),idx);
        tacslice[sliceposition[0],sliceposition[1],ti] = aif[tac_index]+baselineslice[idx];
        endfor
    endfor
  endif
  
  

  ; % Process csf
  csf_idx = where(brainslice eq CSF)
  if n_elements(csf_idx) gt 1 then begin
        mrvalues = mrslice[csf_idx];
        mrvalues = mrvalues - mean(mrvalues);
        mrvalues_bound = 2*sqrt(variance(mrvalues));
        mrvalues[where(mrvalues lt -mrvalues_bound)] = -mrvalues_bound;
        mrvalues[where(mrvalues gt  mrvalues_bound)] =  mrvalues_bound;
        mrvalues = mrvalues/mrvalues_bound;
        baselineslice[csf_idx] = HU_CSF - mrvalues*HU_CSF_VAR;
        for i=0,n_elements(csf_idx)-1 do begin
          idx = csf_idx[i]
          for ti=0,nframe-1  do begin
            sliceposition = array_indices(baselineslice,idx);
            ; [slice_x slice_y] = ind2sub(size(baselineslice),idx);
            tacslice[sliceposition[0],sliceposition[1],ti] = baselineslice[idx];
          endfor
        endfor
  endif


    ; % save baseline slice
   skullslice = skull[*,*,z];
   skullslice[baselineslice ne 0.] = 0.;
   baselineslice = baselineslice + skullslice;   
   baselineall[*,*,z]=baselineslice


    ; % save temporal slices
    for i=0, nframe-1 do begin
      tacslice[*,*,i] = tacslice[*,*,i] + skullslice;
      tacall[*,*,z,i]=tacslice[*,*,i]
    endfor

    ; % save reference perfusion maps
    cbfall[*,*,z]=cbfslice
    cbvall[*,*,z]=cbvslice
    mttall[*,*,z]=mttslice
    delayall[*,*,z]=delayslice
    ttpall[*,*,z]=ttpslice

    ; if n_elements(where(cbfslice gt 0.0)) gt 1 then stop

endfor



stop
; save,filename='tmp_delay.sav', cbfall,tacall,mttall,baselineall,delayall,aif,ts
save,filename='tmp_delay_new.sav', cbfall,tacall,mttall,baselineall,delayall,aif,ts


; plasma_t, plasma_c
  fname= 'plasma_t.txt'
  OPENW,1,fname 
  PRINTF,1,plasma_t   ; ts?
  CLOSE,1
  fname= 'plasma_c.txt'
  OPENW,1,fname 
  PRINTF,1,plasma_c  ; aif?
  CLOSE,1

; create mask image
mask[where(cbfall gt 0.0)] = 1.0
mask1=morph_close(mask,replicate(1,5,5,5))
masksub=mask[*,*,153:157]
delaysub=delayall[*,*,153:157]
baselinesub=baselineall[*,*,153:157]
img=tacall[*,*,153:157,*]

save,filename='tmpsub_delay.sav',masksub,baselinesub,delaysub,img

stop



; comapre generated TACs with ones from matlab
 tacall_ref = tacall
 for i=1, nframe do $
      tacall_ref[*,*,*,i-1]  = readraw('ref/'+nistring(i),ncol,nrow,nplane,'float')
 

 voxelj = [ncol/2,nrow/2,nplane/2]
 plot, tacall[voxelj[0],voxelj[1],voxelj[2],*]
 oplot, tacall_ref[voxelj[0],voxelj[1],voxelj[2],*],col=1



 cbf_ref  = readraw('ref/cbf',ncol,nrow,nplane,'float')
 mtt_ref  = readraw('ref/mtt',ncol,nrow,nplane,'float')
 cbv_ref  = readraw('ref/cbv',ncol,nrow,nplane,'float')
 ttp_ref  = readraw('ref/ttp',ncol,nrow,nplane,'float')





End