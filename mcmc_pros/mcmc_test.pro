

restore,filename='act_img_tomatlab.sav'
tissue_c = reform(img[79,126,1,*])
plasma_t = [0.166667,0.283333,0.366667,0.466667,0.550000,0.616667,0.716667,0.833333,0.950000,1.06667,1.18333,1.31667,1.46667,1.60000,1.78333,1.95000,3.70000,6.76667,11.8833,20.1167,32.0333,45.5833,57.8333]
plasma_c = [0.00000,0.00000,0.00000,0.0375032,0.331662,0.474564,0.592856,0.416336,0.0932330,0.0784635,0.0513207,0.0438369,0.0386417,0.0339184,0.0288506,0.0261968,0.0188040,0.0222774,0.0185163,0.0176231,0.0123044,0.00832840,0.00534190]

help, tissue_c, plasma_c, plasma_t
nsample = double(n_elements(plasma_t))


initialK = [0.5, 1.8, 0.06, 0.002]
initialK = double(initialK)
rwmh_par_scale = double(1)    ; was 1
hmc_step_size  = double(0.5)  ;?
rwmh_n_burnin  = 50000L *2
rwmh_n_draws   = 50000L *2
output = fltarr(rwmh_n_draws, n_elements(initialK))
debug = 1L

; stop
mcmc_tac = 'librwmh_tac.so'

success = call_external(mcmc_tac, 'rwmh_tac', nsample, double(tissue_c), double(plasma_t), double(plasma_c), output, $
                     initialK[0],initialK[1],initialK[2],initialK[3], $
                     rwmh_par_scale,hmc_step_size,rwmh_n_burnin,rwmh_n_draws, debug) 


print, mean(output(*,0)),mean(output(*,1)),mean(output(*,2)),mean(output(*,3))

stop






End