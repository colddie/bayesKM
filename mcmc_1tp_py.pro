
test = 1 ;1

; restore,filename='img.sav'
; restore,filename='img_noise.sav'
; restore, filename= 'actimg_1comp_noise_turku.sav'
; img = niread_nii('actimg_1comp.nii', orientation='RAS')
; img = niread_nii('actimg_fdg.nii', orientation='RAS')
img = niread_nii('actimg_way.nii', orientation='RAS')
restore,filename='mask.sav'

; plasma_t = [0.166667,0.283333,0.366667,0.466667,0.550000,0.616667,0.716667,0.833333,0.950000,1.06667,1.18333,1.31667,1.46667,1.60000,1.78333,1.95000,3.70000,6.76667,11.8833,20.1167,32.0333,45.5833,57.8333]
; plasma_c = [0.00000,0.00000,0.00000,0.0375032,0.331662,0.474564,0.592856,0.416336,0.0932330,0.0784635,0.0513207,0.0438369,0.0386417,0.0339184,0.0288506,0.0261968,0.0188040,0.0222774,0.0185163,0.0176231,0.0123044,0.00832840,0.00534190]
tmp = read_ascii('plasma_t.txt')
plasma_t = tmp.FIELD01 ;/20
tmp = read_ascii('plasma_c.txt')
plasma_c = tmp.FIELD01
nsample = double(n_elements(plasma_t))
weight = total(total(total(img,1),1),1) / total(img)*nsample   ; ought to be duration but we infer that from the image in simulation
; weight   = (dta->weight)[i]	=  (dta->tissue_c)[i]/arma::accu(dta->tissue_c)*nsample;
; weight   = double(nsample) + 1.0
model = 2   ; 1,2,3
tracer = 'way'   ; 'way','fdg', 'fluropirdize', 'fdopa'
mcmc = 0;  rwmh, hmc


case model of 
  1: begin     ;K1,k2
    if tracer eq 'fluropirdize' then begin  
      initialK = [0.3, 1.6]            
      lb = [0., 0.]
      ub = [2., 4.]    ;[2., 4.]
    endif
  end

  2: begin    ;K1,k2,k3,k4
    if tracer eq 'fdg' then begin
      initialK = [0.5, 1.8, 0.06, 0.002]
      lb = [0., 0., 0., 0.]
      ub = [2., 4., 0.12, 0.01]
    endif
    if tracer eq 'way' then begin
      initialK = [0.14, 0.25, 0.18, 0.02]
      lb = [0., 0., 0., 0.]
      ub = [1., 1., 0.5, 0.1]
    endif
    if tracer eq 'fdopa' then stop
  end

  3: begin ;R1,k2,BP
    if tracer eq 'way' then begin
      initialK = [0.5, 1.8, 0.06]           
      lb = [0., 0., 0.]
      ub = [2., 4., 0.12]
    endif
  end
endcase

initialK = double(initialK)
lb = double(lb)
ub = double(ub)
rwmh_par_scale = double(0.1)    ; was 1
hmc_step_size  = double(0.01)  ;?
rwmh_n_burnin  = 10000L              ; ?
rwmh_n_draws   = 10000L 
output = fltarr(rwmh_n_draws, n_elements(initialK))
debug = 1L ;1L

mcmc_tac = 'librwmh_tac.so'




ncol   = (size(img))[1]
nrow   = (size(img))[2]
nplane = (size(img))[3]
meanimg = fltarr(ncol, nrow, nplane, n_elements(initialK))
varimg  = fltarr(ncol, nrow, nplane, n_elements(initialK))


tt = systime(1)
for iplane = 0, nplane-1 do begin
    print, 'plane' + nistring(iplane)
  for irow = 0, nrow-1 do begin
    for icol = 0, ncol-1 do begin

    if test then begin
      icol = 79
      irow = 126
      iplane = 1
    endif

    if mask[icol, irow, iplane] eq 0.0 then continue
        tissue_c = reform(img[icol,irow,iplane,*])            ;reform(img[79,126,1,*])

        case model of 
          1: begin ;1 tissue compartment model
            success = call_external(mcmc_tac, 'rwmh_tac_1tpc', nsample, double(tissue_c), double(plasma_t), double(plasma_c), $
                                weight, output, $
                                initialK[0],initialK[1], $
                                rwmh_par_scale,hmc_step_size,rwmh_n_burnin,rwmh_n_draws, lb[0], lb[1], ub[0], ub[1], debug, mcmc) 


            print, mean(output(*,0)),sqrt(stddev(output(*,0))),mean(output(*,1)),sqrt(stddev(output(*,1)))
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
          end
          2: begin  ;2 tissue compartment model
            model=Python.import('stanidl')              
            output=model.rwmh_tac_2tpc(tissue_c,plasma_c,plasma_t,weight,initialK, $
                          lb,ub,rwmh_par_scale,hmc_step_size,rwmh_n_burnin,rwmh_n_draws) 








            ; success = call_external(mcmc_tac, 'rwmh_tac_2tpc', nsample, double(tissue_c), double(plasma_t), double(plasma_c), weight, output, $
            ;                 initialK[0],initialK[1],initialK[2],initialK[3], $
            ;                 rwmh_par_scale,hmc_step_size,rwmh_n_burnin,rwmh_n_draws, lb[0],lb[1],lb[2],lb[3], ub[0],ub[1],ub[2],ub[3], model,debug, mcmc) 


            ; print, mean(output(*,0)),sqrt(stddev(output(*,0))),mean(output(*,1)),sqrt(stddev(output(*,1))),$
            ;        mean(output(*,2)),sqrt(stddev(output(*,2))),mean(output(*,3)),sqrt(stddev(output(*,3)))
            ; meanimg[icol,irow,iplane,0] = mean(output(*,0))
            ; meanimg[icol,irow,iplane,1] = mean(output(*,1))
            ; meanimg[icol,irow,iplane,2] = mean(output(*,2))
            ; meanimg[icol,irow,iplane,3] = mean(output(*,3))
            ; varimg[icol,irow,iplane,0]  = stddev(output(*,0))
            ; varimg[icol,irow,iplane,1]  = stddev(output(*,1))
            ; varimg[icol,irow,iplane,2]  = stddev(output(*,2))
            ; varimg[icol,irow,iplane,3]  = stddev(output(*,3))
          end
          3: begin  ;simplified reference tissue compartment model
            success = call_external(mcmc_tac, 'rwmh_tac_2tpc', nsample, double(tissue_c), double(plasma_t), double(plasma_c), weight, output, $
                            initialK[0],initialK[1],initialK[2], $
                            rwmh_par_scale,hmc_step_size,rwmh_n_burnin,rwmh_n_draws, lb[0],lb[1],lb[2], ub[0],ub[1],ub[2], model,debug, mcmc) 


            print, mean(output(*,0)),sqrt(stddev(output(*,0))),mean(output(*,1)),sqrt(stddev(output(*,1))),$
                   mean(output(*,2)),sqrt(stddev(output(*,2)))
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
  save, filename='mcmc_noiseresults'+nistring(iplane)+'.sav', meanimg, varimg
end
print, 'elapsed time... ' + nistring(systime(1)-tt) + 'seconds'









stop

cgHistoplot,output(*,0), POLYCOLOR='olive', /FILL, YRANGE=yrange, $
    XRANGE=xrange, BINSIZE=binsize, HISTDATA=h, POSITION=position, LOCATIONS=loc
binCenters = loc + (binsize / 2.0)
yfit = GaussFit(binCenters, h, coeff, NTERMS=3)
thick = (!D.Name EQ 'PS') ? 5 : 3
cgPlot, binCenters, yfit, COLOR='dodger blue', THICK=thick, /OVERPLOT


window
cgHistoplot,output(*,1), POLYCOLOR='red', /FILL, YRANGE=yrange, $
    XRANGE=xrange, BINSIZE=binsize, HISTDATA=h, POSITION=position, LOCATIONS=loc
binCenters = loc + (binsize / 2.0)
yfit = GaussFit(binCenters, h, coeff, NTERMS=3)
thick = (!D.Name EQ 'PS') ? 5 : 3
cgPlot, binCenters, yfit, COLOR='dodger blue', THICK=thick, /OVERPLOT




End