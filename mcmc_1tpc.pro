
test = 0
debug = 0

if 1 then begin

  pad = '/home/tsun/bin/fsl/install/src/fabber_core/fabber_pet_c1/fluropirdizeC1rebin0.5/'
  ; restore,filename='img.sav'
  ; restore,filename='img_noise.sav'
  ; restore, filename= 'actimg_1comp_noise_turku.sav'
  ; img = niread_nii('actimg_1comp.nii', orientation='RAS')
  ; img = niread_nii('actimg_fdg.nii', orientation='RAS')
  ; img = niread_nii('actimg_wayformcmc.nii', orientation='RAS')
  ; img = niread_nii('actimg_waybound.nii', orientation='RAS')
  ; img = niread_nii('actimg_way_0.3_0.01_noise.nii', orientation='RAS')
  ; img = niread_nii('actimg_way_0_0.005_noise.nii', orientation='RAS')
  ; img = niread_nii('actimg_way.nii', orientation='RAS')
  ; img = niread_nii('actimg_waysrtm.nii', orientation='RAS')
  ; img = niread_nii('actimg_wayC2.nii', orientation='RAS')
  ; img = niread_nii(pad+'actimg_fluropirdizeC1.nii', orientation='RAS')
  img = niread_nii(pad+'actimg_fluropirdizeC1_0.00000_0.00500000_noise.nii', orientation='RAS')
  ; img = niread_nii('actimg_waysrtm.nii', orientation='RAS')
  ; img = niread_nii('actimg_wayC2_0_0.005_noise.nii', orientation='RAS')
  ; restore,filename='actimg_way.sav'
  ; img=tmpimg
  restore,filename='/home/tsun/bin/fsl/install/src/fabber_core/fabber_pet_c1/mask.sav'
  
  ; plasma_t = [0.166667,0.283333,0.366667,0.466667,0.550000,0.616667,0.716667,0.833333,0.950000,1.06667,1.18333,1.31667,1.46667,1.60000,1.78333,1.95000,3.70000,6.76667,11.8833,20.1167,32.0333,45.5833,57.8333]
  ; plasma_c = [0.00000,0.00000,0.00000,0.0375032,0.331662,0.474564,0.592856,0.416336,0.0932330,0.0784635,0.0513207,0.0438369,0.0386417,0.0339184,0.0288506,0.0261968,0.0188040,0.0222774,0.0185163,0.0176231,0.0123044,0.00832840,0.00534190]
  tmp = read_ascii(pad+'plasma_t.txt')
  plasma_t = reform(tmp.FIELD1,n_elements(tmp.FIELD1),1) 
  plasma_t = plasma_t[0:(size(img))[4]-1]
  tmp = read_ascii(pad+'plasma_c.txt')
  ; tmp = read_ascii('ref_tissuec.txt')
  plasma_c = reform(tmp.FIELD1,n_elements(tmp.FIELD1),1) 
  plasma_c = plasma_c[0:(size(img))[4]-1]
  plasma_c[0] = 0.
endif

if 0 then begin
  
  pad = '/home/tsun/work/bayesKM/DynamicBrainPhantom/'
  ; restore,filename='tmpsub.sav'
  restore,filename='tmpsub_delay.sav'
  ; restore,filename=pad+'tmpsub.sav'
  ; restore,filename=pad+'masksub.sav'
  ; restore,filename=pad+'baselinesub.sav'

  for iframe = 0,(size(img))[4]-1 do  img[*,*,*,iframe] -= baselinesub    ; subtract the baseline image
  tstart = 0 
  tstop = 49
  t = indgen(50)*1;     ;0:1:49
  to = indgen(25)*2 +1  ;1:2:49;
  ts = indgen(tstop*10+1)*0.1 
    
  ; plasma_t = [0.166667,0.283333,0.366667,0.466667,0.550000,0.616667,0.716667,0.833333,0.950000,1.06667,1.18333,1.31667,1.46667,1.60000,1.78333,1.95000,3.70000,6.76667,11.8833,20.1167,32.0333,45.5833,57.8333]
  ; plasma_c = [0.00000,0.00000,0.00000,0.0375032,0.331662,0.474564,0.592856,0.416336,0.0932330,0.0784635,0.0513207,0.0438369,0.0386417,0.0339184,0.0288506,0.0261968,0.0188040,0.0222774,0.0185163,0.0176231,0.0123044,0.00832840,0.00534190]
  tmp = read_ascii(pad+'plasma_t.txt')
  plasma_t = reform(tmp.FIELD1,n_elements(tmp.FIELD1),1) 
  plasma_t = plasma_t[0:n_elements(ts)-1]
  tmp = read_ascii(pad+'plasma_c.txt')
  ; tmp = read_ascii('ref_tissuec.txt')
  plasma_c = reform(tmp.FIELD1,n_elements(tmp.FIELD1),1) 
  plasma_c = plasma_c[0:n_elements(ts)-1]
  plasma_c[0] = 0.

endif




nsample = n_elements(plasma_t)
; weight = total(total(total(img,1),1),1) / total(img)*nsample   ; ought to be duration but we infer that from the image in simulation
; weight   = (dta->weight)[i]	=  (dta->tissue_c)[i]/arma::accu(dta->tissue_c)*nsample;
; weight   = double(nsample) + 1.0
weight = (plasma_t -shift(plasma_t,1)) / plasma_t[n_elements(plasma_t)-1] * nsample
weight[0] = 0.
; weight = fltarr(nsample) + 1.
model = 1L   ; 1,2,3
tracer = 'fluropirdize'   ; 'way','fdg', 'fluropirdize', 'fdopa'
mcmc = 0L;  rwmh, hmc
useprior = 0LL

case model of 
  1: begin     ;K1,k2
    if tracer eq 'fluropirdize' then begin  
      initialK = [0.3, 0.4]            
      lb = [0., 0.]
      ub = [2., 4.]    ;[2., 4.]
    endif
  end

  2: begin    ;K1,k2,k3,k4
    if tracer eq 'fdg' then begin
      initialK = [0.5, 1.8, 0.06, 0.002]
      lb = [0., 0., 0., 0.]
      ub = [2., 4., 0.12, 0.01]
      ; prior = initialK
    endif
    if tracer eq 'way' then begin
      initialK = [0.14,0.4,0.2,0.05]   ;[0.14, 0.29, 0.18, 0.02]
      lb = [0., 0., 0., 0.]
      ub = [1., 1., 0.5, 0.1]
      ; prior = initialK
    endif
    if tracer eq 'fdopa' then stop

    if useprior then begin
      phantomK1 =  niread_nii('phantomK1true_way.nii', orientation='RAS')
      phantomk2 =  niread_nii('phantomk2true_way.nii', orientation='RAS')
      phantomk3 =  niread_nii('phantomk3true_way.nii', orientation='RAS')
      phantomk4 =  niread_nii('phantomk4true_way.nii', orientation='RAS')
    endif

  end

  3: begin ;R1,k2,BP
    if tracer eq 'way' then begin
      initialK = [1.0, 0.1, 3.,0]           
      lb = [0., 0., 0.,0]
      ub = [2., 1., 15.,0]
    endif

      if usepriorimg then begin
      phantomK1 =  niread_nii('srtmK1true_way.nii', orientation='RAS')
      phantomk2 =  niread_nii('srtmk2true_way.nii', orientation='RAS')
      phantomk3 =  niread_nii('srtmk3true_way.nii', orientation='RAS')
      phantomk4 =  niread_nii('srtmk4true_way.nii', orientation='RAS')
    endif
  end

  4: begin ;R1,k2,k3,k4
    if tracer eq 'way' then begin
      initialK = [1.0, 0.29, 0.18, 0.02]           
      lb = [0., 0., 0., 0.]
      ub = [2., 1., 0.5, 0.1]
    endif
  end


  ; 5,6: begin
  ; end

  7: begin ;CBF, MTT, delay
      initialK = [30., 10., 0.]            
      lb = [0., 0., 0.]
      ub = [100., 100., 10.]    ;[2., 4.]     
  end
    
endcase



; to reduce the variance, increase samples, decrease scale, extend bounds, TRAINING?
initialK = double(initialK)
lb = double(lb)
ub = double(ub)
rwmh_par_scale = double(0.006)    ; was 0.06 for 1tpc 
hmc_step_size  = double(0.001)
rwmh_n_burnin  = 10000L *2          ; ?
rwmh_n_draws   = 10000L 
output = fltarr(rwmh_n_draws, n_elements(initialK))
prior = dblarr(n_elements(initialK))
sens = double(1e0)
mcmc_tac = '/home/tsun/bin/mcmc-master/tests/example/librwmh_tac.so'

; 0.06 for 2tpc
; rwmh_n_burnin  = 10000L *2           
; rwmh_n_draws   = 10000L *4
; 0.001 for srtm
; rwmh_par_scale = double(0.001)    ; was 1 larger the less acceptance rate
; hmc_step_size  = double(0.001)  ;?
; rwmh_n_burnin  = 10000L *8           
; rwmh_n_draws   = 10000L *8

ncol   = (size(img))[1]
nrow   = (size(img))[2]
nplane = (size(img))[3]
meanimg = fltarr(ncol, nrow, nplane, n_elements(initialK))
varimg  = fltarr(ncol, nrow, nplane, n_elements(initialK))
meanimg1 = fltarr(ncol, nrow, nplane, n_elements(initialK))
varimg1  = fltarr(ncol, nrow, nplane, n_elements(initialK))

tt = systime(1)
for iplane = 0, nplane-1 do begin
    print, 'plane' + nistring(iplane)

  ; mpi_rank, rank
  ; start = 0
  ; stop  = nrow
  ; ;; divide task among MPI ranks
  ; mpi_helper, start, stop, i_start, i_stop

  for irow = 0, nrow-1 do begin
    for icol = 0, ncol-1 do begin

    if test then begin
      icol = 128 ;79
      irow = 173 ;126
      iplane = 1
    endif

    if mask[icol, irow, iplane] eq 0.0 then continue      ;masksub

        tissue_c = reform(img[icol,irow,iplane,*])            ;reform(img[79,126,1,*])

        case model of 
          1: begin ;1 tissue compartment model
            if useprior then $
              prior = [ phantomK1[icol,irow,iplane], phantomk2[icol,irow,iplane], 0.0, 0.0 ]
              
            success = call_external(mcmc_tac, 'rwmh_tac_2tpc', long(nsample), double(tissue_c), double(plasma_t), double(plasma_c), $   ; was rwmh_tac_1tpc
                            double(weight), double(prior), output, $
                            double(initialK), double(lb), double(ub), $
                            rwmh_par_scale,hmc_step_size,rwmh_n_burnin,rwmh_n_draws, model,debug, mcmc,useprior,sens) 

            print, mean(output(*,0)),stddev(output(*,0)),mean(output(*,1)),stddev(output(*,1))
           
            meanimg[icol,irow,iplane,0] = mean(output(*,0))
            meanimg[icol,irow,iplane,1] = mean(output(*,1))
            varimg[icol,irow,iplane,0]  = stddev(output(*,0))
            varimg[icol,irow,iplane,1]  = stddev(output(*,1))

            ; pdf0 = histogram(output(*,0), LOCATIONS=xbin0, nbins=100)
            ; pdf1 = histogram(output(*,1), LOCATIONS=xbin1, nbins=100)
            ; yfit0 = GaussFit(xbin0, pdf0, coeff, NTERMS=3)
            ; yfit1 = GaussFit(xbin1, pdf1, coeff, NTERMS=3)
            ; p0 = plot(xbin0, yfit0, THICK=2)
            ; p1 = PLOT(xbin1, yfit1, THICK=2)

            ; output1 = double(fltarr(12))   ; 3+3+6
            ; debug  = 0
            ; llsq_model = 0
            ; t1 = plasma_t
            ; t0 = [[0], t1[0:n_elements(t1)-2]] 
            ; isweight = 1
            ; ; weights = fltarr(frameNr) + 1.0
            ; weights = (t1-t0)/total(t1-t0)  ; * nsample
            ; ; directbp = 1
            ; def_pmin = [[lb],0.0]
            ; def_pmax = [[ub],0.0]  
            ; fVb = 0.0     ; fix Vb to 0 as plasma presented
            ; doSD = 1
            ; doCL = 0
            ; bootstrapIter = 200  ; has to be larger than 100!
            ; matrix = double(fltarr(bootstrapIter*4)) ; change with num_param

            ; patlog = '/home/tsun/bin/tpcclib-master/build/bin/libmtga_idl.so'
            ; success = call_external(patlog, 'tcm1_idl', long(nsample), double(t0), double(t1), double(tissue_c), $
            ;                 double(plasma_c),output1,long(debug), $
            ;                 long(isweight),double(weights),double(def_pmin),double(def_pmax), $
            ;                 double(fVb),long(doSD),long(doCL), $
            ;                 long(bootstrapIter),matrix) 

            ; print, output1[0], output1[5], output1[1], output1[6]
            ; meanimg1[icol,irow,iplane,0] = output1[0]
            ; meanimg1[icol,irow,iplane,1] = output1[1]
            ; varimg1[icol,irow,iplane,0]  = output1[5]
            ; varimg1[icol,irow,iplane,1]  = output1[6]

          end
          2: begin  ;2 tissue compartment model
           if useprior then $
              prior = [ phantomK1[icol,irow,iplane], phantomk2[icol,irow,iplane], phantomk3[icol,irow,iplane], phantomk4[icol,irow,iplane] ]

            success = call_external(mcmc_tac, 'rwmh_tac_2tpc', nsample, double(tissue_c), double(plasma_t), double(plasma_c), $
                            double(weight), double(prior), output, $
                            double(initialK), double(lb), double(ub), $
                            rwmh_par_scale,hmc_step_size,rwmh_n_burnin,rwmh_n_draws, model,debug, mcmc,useprior,sens) 

            print, mean(output(*,0)),stddev(output(*,0)),mean(output(*,1)),stddev(output(*,1)),$
                   mean(output(*,2)),stddev(output(*,2)),mean(output(*,3)),stddev(output(*,3))
            meanimg[icol,irow,iplane,0] = mean(output(*,0))
            meanimg[icol,irow,iplane,1] = mean(output(*,1))
            meanimg[icol,irow,iplane,2] = mean(output(*,2))
            meanimg[icol,irow,iplane,3] = mean(output(*,3))
            varimg[icol,irow,iplane,0]  = stddev(output(*,0))
            varimg[icol,irow,iplane,1]  = stddev(output(*,1))
            varimg[icol,irow,iplane,2]  = stddev(output(*,2))
            varimg[icol,irow,iplane,3]  = stddev(output(*,3))
          end
          3: begin  ;simplified reference tissue compartment model
            if usepriorimg then begin
              initialK = [ phantomK1[icol,irow,iplane], phantomk2[icol,irow,iplane], phantomk3[icol,irow,iplane], phantomk4[icol,irow,iplane] ]
              ; lb = initialK-[0.2,0.005,1,0]
              ; ub = initialK+[0.2,0.005,1,0]
            endif
            
            if useprior then $
              prior = [ phantomK1[icol,irow,iplane], phantomk2[icol,irow,iplane], phantomk3[icol,irow,iplane], phantomk4[icol,irow,iplane] ]
            
            success = call_external(mcmc_tac, 'rwmh_tac_2tpc', nsample, double(tissue_c), double(plasma_t), double(plasma_c), $
                            double(weight),double(prior), output, $
                            double(initialK), double(lb), double(ub), $
                            rwmh_par_scale,hmc_step_size,rwmh_n_burnin,rwmh_n_draws, model,debug, mcmc, useprior,sens) 


            print, mean(output(*,0)),stddev(output(*,0)),mean(output(*,1)),stddev(output(*,1)),$
                   mean(output(*,2)),sqrt(stddev(output(*,2)))
            meanimg[icol,irow,iplane,0] = mean(output(*,0))
            meanimg[icol,irow,iplane,1] = mean(output(*,1))
            meanimg[icol,irow,iplane,2] = mean(output(*,2))
            varimg[icol,irow,iplane,0]  = stddev(output(*,0))
            varimg[icol,irow,iplane,1]  = stddev(output(*,1))
            varimg[icol,irow,iplane,2]  = stddev(output(*,2))
          end
          4: begin  ;full reference tissue compartment model
            if useprior then $
              prior = [ phantomK1[icol,irow,iplane], phantomk2[icol,irow,iplane], phantomk3[icol,irow,iplane], phantomk4[icol,irow,iplane] ]

            success = call_external(mcmc_tac, 'rwmh_tac_2tpc', nsample, double(tissue_c), double(plasma_t), double(plasma_c), $
                            double(weight),double(prior), output, $
                            double(initialK), double(lb), double(ub), $
                            rwmh_par_scale,hmc_step_size,rwmh_n_burnin,rwmh_n_draws, model,debug, mcmc, useprior,sens) 


            print, mean(output(*,0)),stddev(output(*,0)),mean(output(*,1)),stddev(output(*,1)),$
                   mean(output(*,2)),stddev(output(*,2)),mean(output(*,3)),stddev(output(*,3))
            meanimg[icol,irow,iplane,0] = mean(output(*,0))
            meanimg[icol,irow,iplane,1] = mean(output(*,1))
            meanimg[icol,irow,iplane,2] = mean(output(*,2))
            meanimg[icol,irow,iplane,3] = mean(output(*,3))
            varimg[icol,irow,iplane,0]  = stddev(output(*,0))
            varimg[icol,irow,iplane,1]  = stddev(output(*,1))
            varimg[icol,irow,iplane,2]  = stddev(output(*,2))
            varimg[icol,irow,iplane,3]  = stddev(output(*,3))
          end

                ;  5,6: begin ; patlak, logan
         
      ;  end

        7: begin  ; perfusion CT
            if test eq 1 then begin
              icol=103
              irow=151
              iplane = 1
            endif
            tissue_c = reform(img[icol,irow,iplane,*])
            tac1 =  interpol(tissue_c,t,ts, /spline);
            success = call_external(mcmc_tac, 'rwmh_tac_2tpc', long(nsample), double(tac1), double(plasma_t), double(plasma_c), $   ; was rwmh_tac_1tpc
                            double(weight), double(prior), output, $
                            double(initialK), double(lb), double(ub), $
                            rwmh_par_scale,hmc_step_size,rwmh_n_burnin,rwmh_n_draws, model,debug, mcmc,useprior,sens) 

            print, mean(output(*,0)),stddev(output(*,0)),mean(output(*,1)),stddev(output(*,1))
           
            meanimg[icol,irow,iplane,0] = mean(output(*,0))
            meanimg[icol,irow,iplane,1] = mean(output(*,1))
            meanimg[icol,irow,iplane,2] = mean(output(*,2))
            varimg[icol,irow,iplane,0]  = stddev(output(*,0))
            varimg[icol,irow,iplane,1]  = stddev(output(*,1))
            varimg[icol,irow,iplane,2]  = stddev(output(*,2))
        end

      endcase

        if test then     stop
    end
  end
  save, filename='mcmc_noiseresults'+nistring(iplane)+'_'+nistring(model)+tracer+nistring(mcmc)+nistring(useprior)+'.sav', meanimg, varimg
end
print, 'elapsed time... ' + nistring(systime(1)-tt) + 'seconds'









stop

cgHistoplot,output(*,0), POLYCOLOR='olive', /FILL, YRANGE=yrange, $
    XRANGE=xrange, BINSIZE=binsize, HISTDATA=h, POSITION=position, LOCATIONS=loc
binCenters = loc + (binsize / 2.0)
yfit = GaussFit(binCenters, h, coeff, NTERMS=3, sigma=s)
thick = (!D.Name EQ 'PS') ? 5 : 3
cgPlot, binCenters, yfit, COLOR='dodger blue', THICK=thick, /OVERPLOT


window
cgHistoplot,output(*,1), POLYCOLOR='red', /FILL, YRANGE=yrange, $
    XRANGE=xrange, BINSIZE=binsize, HISTDATA=h, POSITION=position, LOCATIONS=loc
binCenters = loc + (binsize / 2.0)
yfit = GaussFit(binCenters, h, coeff, NTERMS=3)
thick = (!D.Name EQ 'PS') ? 5 : 3
cgPlot, binCenters, yfit, COLOR='dodger blue', THICK=thick, /OVERPLOT

window
cgHistoplot,output(*,2), POLYCOLOR='red', /FILL, YRANGE=yrange, $
    XRANGE=xrange, BINSIZE=binsize, HISTDATA=h, POSITION=position, LOCATIONS=loc
binCenters = loc + (binsize / 2.0)
yfit = GaussFit(binCenters, h, coeff, NTERMS=3)
thick = (!D.Name EQ 'PS') ? 5 : 3
cgPlot, binCenters, yfit, COLOR='dodger blue', THICK=thick, /OVERPLOT

window
cgHistoplot,output(*,3), POLYCOLOR='red', /FILL, YRANGE=yrange, $
    XRANGE=xrange, BINSIZE=binsize, HISTDATA=h, POSITION=position, LOCATIONS=loc
binCenters = loc + (binsize / 2.0)
yfit = GaussFit(binCenters, h, coeff, NTERMS=3)
thick = (!D.Name EQ 'PS') ? 5 : 3
cgPlot, binCenters, yfit, COLOR='dodger blue', THICK=thick, /OVERPLOT


End