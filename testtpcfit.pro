; test tpc pixel-wised patlak, logan with given TAC, input function

; first we test the original turku program to see what happens:
; patlak --verbose=10 dynamic_1tac.tac plasma.bld 30 90 ut1234.res       RETURN KI
; patlak --verbose=10 dynamic_1tac.tac ambf_dx_Pl001 30 90 ut1234.res     RETURN KIREF
; logan logan_1tac.tac plasma.bld 30 90 ut2345.res  RETURN DV
; logan -k2=0.163 logan_1tac.tac ref_0_Pl01 30 90 ut2345.res   RETURN DVR/BP
; regfur --verbose=10 plasma.bld  dynamic_1tac.tac30 90 ut1234.res     ; plasma first!
; regfur plasma.bld  dynamic_1tac.tac30 90 ut1234.res -curve=ut.res


; USER is responsible to make the tac and input mathc in frames and time!
; t0 - frame start time
; t1 - frame end time
; tac - time-activity-curve 
; ctt - input function

; in case metabolic rate is of interest- MRf=100.*Ca/(density*LC);
;   " -Ca=<value>",
;   "     Concentration of native substrate in arterial plasma (mM),",
;   "     for example plasma glucose in [18F]FDG studies.",
;   "     With this option the metabolic rate (umol/(min*100 g)) is calculated.",
;   " -LC=<value>",
;   "     Lumped Constant in MR calculation; default is 1.0.",
;   " -density=<value>",
;   "     Tissue density in MR calculation; default is 1.0 g/ml.",

  time0 = systime(1)
  useref = 0;   use reference region as input?
  testmode = 'pCT';   '1TCM','2TCM','patlak','logan','fur','mrtm','sim_patlak','sim_logan','SRTM','2TCMrev','MBF', 'pCT'

case testmode of
    'patlak': begin
    ;*****************************************************************************/
    ; PATLAK
    ;*****************************************************************************/
        t0 = [0.,0.25,0.5,0.75,1.0,1.5,2.0,3.0,4.0,5.0,10.0,15.0,20.0,25.0,30.,35,40.,45.,50.,55,60,70,80,90,100,110];
        t1 = [0.25,0.5,0.75,1.0,1.5,2.0,3.0,4.0,5.0,10.0,15.0,20.0,25.0,30.,35,40.,45.,50.,55,60,70,80,90,100,110,120]; 
        ctt = [0.06415,1.879,84.8,66.15,41.86,26.86,21.36,19.87,17.15,13.54,10.89,9.305,8.147,7.294,6.655,6.167,5.786,5.479,$
        5.225,5.008,4.729,4.406,4.114,3.836,3.564,3.295];
        tac = [2.99e-4,5.463e-3,4.679e-1,1.523e0,2.443563e+00,3.236325e+00,4.015189e+00,4.897500e+00,5.656930e+00, $
        7.241067e+00,9.052315e+00,1.013719e+01,1.076269e+01,1.107761e+01,1.118487e+01,1.115529e+01,1.103736e+01,1.086400e+01,1.065737e+01,1.043218e+01, $
        1.008137e+01,9.610046e+00,9.152024e+00,8.711683e+00,8.287338e+00,7.875306e+00];

        if useref then begin
        
        endif

        frameNr= n_elements(t0)
        tstart = 30
        tstop  = 90 
        output = double(fltarr(5)) 
        debug  = 10
        llsq_model = 0
        isweight = 0;

        patlog = '/home/tsun/bin/tpcclib-master/build/bin/libmtga_idl.so'
        success = call_external(patlog, 'patlak_idl', long(frameNr), double(t0), double(t1), double(tac), $
                        double(ctt),double(tstart),double(tstop),double(output),long(debug),long(llsq_model),$
                        long(isweight)) 

        if useref then begin
            
        endif

  end

    'logan': begin
    ;*****************************************************************************/
    ; LOGAN
    ;*****************************************************************************/
        ; model := C3S K1 := 0.3 k2 := 0.15 k3 := 0 k4 := 0 k5 := 0 k6 := 0
        t0 = [0.,0.25,0.5,0.75,1.0,1.5,2.0,3.0,4.0,5.0,10.0,15.0,20.0,25.0,30.,35,40.,45.,50.,55,60,70,80,90,100,110];
        t1 = [0.25,0.5,0.75,1.0,1.5,2.0,3.0,4.0,5.0,10.0,15.0,20.0,25.0,30.,35,40.,45.,50.,55,60,70,80,90,100,110,120]; 
        ctt = [0.06415,1.879,84.8,66.15,41.86,26.86,21.36,19.87,17.15,13.54,10.89,9.305,8.147,7.294,6.655,6.167,5.786,5.479,$
        5.225,5.008,4.729,4.406,4.114,3.836,3.564,3.295];
        tac= [1.775532e-03,3.239127e-02,2.774877e+00,8.931059e+00,1.391659e+01,1.773189e+01,2.070444e+01,2.346612e+01,$
        2.531194e+01,2.660934e+01,2.522953e+01,2.245561e+01,1.972853e+01,1.740301e+01,1.553417e+01,1.406894e+01,1.292658e+01,$
        1.202974e+01,1.131428e+01,1.073048e+01,1.004793e+01,9.288321e+00,8.647856e+00,8.065366e+00,7.510165e+00,6.966928e+00]

        if useref then begin
        ctt= [1.775532e-03,3.239127e-02,2.774877e+00,8.931059e+00,1.391659e+01,1.773189e+01,2.070444e+01,2.346612e+01,$
        2.531194e+01,2.660934e+01,2.522953e+01,2.245561e+01,1.972853e+01,1.740301e+01,1.553417e+01,1.406894e+01,1.292658e+01,$
        1.202974e+01,1.131428e+01,1.073048e+01,1.004793e+01,9.288321e+00,8.647856e+00,8.065366e+00,7.510165e+00,6.966928e+00]
        tac = [2.358366e-03,4.303226e-02,3.684540e+00,1.179610e+01,1.816910e+01,2.279551e+01,2.601503e+01,2.873279e+01, $
        3.032440e+01,3.030020e+01,2.726783e+01,2.387788e+01,2.098287e+01,1.863618e+01,1.677277e+01,1.530138e+01,1.413542e+01,$
        1.320196e+01,1.244281e+01,1.181292e+01,1.106495e+01,1.022359e+01,9.511543e+00,8.865074e+00,8.250660e+00,7.650896e+00]; 
        endif

        frameNr= n_elements(t0)
        tstart = 30
        tstop  = 90 
        output = double(fltarr(5)) 
        debug  = 10
        llsq_model = 0
        k2 = -0.163          ; <0 if not using reference region
        isweight = 1
        weights = (t1-t0)/total(t1-t0)
        ; weights = fltarr(frameNr) + 1.0
        logan_mode = 0            ; approximation with replacemnt of CR(t) as CP(t)

        patlog = '/home/tsun/bin/tpcclib-master/build/bin/libmtga_idl.so'
        success = call_external(patlog, 'logan_idl', long(frameNr), double(t0), double(t1), double(tac), $
                        double(ctt),double(tstart),double(tstop),double(output),long(debug),long(llsq_model),double(k2),$
                        long(isweight),double(weights),long(logan_mode)) 


        if k2 le 0 then DV = output[0]
        if k2 gt 0 then DVR = output[0]
        ;    if k2 le 0 then BPp = output[0]   ;?

        print, output
  end

    'fur': begin
    ;*****************************************************************************/
    ; FUR
    ;*****************************************************************************/
        t0 = [0.,0.25,0.5,0.75,1.0,1.5,2.0,3.0,4.0,5.0,10.0,15.0,20.0,25.0,30.,35,40.,45.,50.,55,60,70,80,90,100,110];
        t1 = [0.25,0.5,0.75,1.0,1.5,2.0,3.0,4.0,5.0,10.0,15.0,20.0,25.0,30.,35,40.,45.,50.,55,60,70,80,90,100,110,120]; 
        ctt = [0.06415,1.879,84.8,66.15,41.86,26.86,21.36,19.87,17.15,13.54,10.89,9.305,8.147,7.294,6.655,6.167,5.786,5.479,$
        5.225,5.008,4.729,4.406,4.114,3.836,3.564,3.295];
        tac = [2.99e-4,5.463e-3,4.679e-1,1.523e0,2.443563e+00,3.236325e+00,4.015189e+00,4.897500e+00,5.656930e+00, $
        7.241067e+00,9.052315e+00,1.013719e+01,1.076269e+01,1.107761e+01,1.118487e+01,1.115529e+01,1.103736e+01,1.086400e+01,1.065737e+01,1.043218e+01, $
        1.008137e+01,9.610046e+00,9.152024e+00,8.711683e+00,8.287338e+00,7.875306e+00];


        frameNr= n_elements(t0)
        tstart = 30
        tstop  = 90 
        frame_used = n_elements(where( (t0+t1)/2 gt tstart and (t0+t1)/2 lt tstop ))
        output = double(fltarr(10+1)) 
        debug  = 10
        fur_mode = 0

        patlog = '/home/tsun/bin/tpcclib-master/build/bin/libmtga_idl.so'
        success = call_external(patlog, 'regfur_idl', long(frameNr), double(t0), double(tac), $   ;, double(t1)
                        double(ctt),double(tstart),double(tstop),double(output),long(debug),long(fur_mode) ) 
    

        print, output
  end

    'MRTM': begin
    
    ;*****************************************************************************/
    ; MRTM
    ;*****************************************************************************/
    ; model := C3S K1 := 0.3 k2 := 0.15 k3 := 0 k4 := 0 k5 := 0 k6 := 0
        t0 = [0.,0.25,0.5,0.75,1.0,1.5,2.0,3.0,4.0,5.0,10.0,15.0,20.0,25.0,30.,35,40.,45.,50.,55,60,70,80,90,100,110];
        t1 = [0.25,0.5,0.75,1.0,1.5,2.0,3.0,4.0,5.0,10.0,15.0,20.0,25.0,30.,35,40.,45.,50.,55,60,70,80,90,100,110,120]; 
        ctt = [0.06415,1.879,84.8,66.15,41.86,26.86,21.36,19.87,17.15,13.54,10.89,9.305,8.147,7.294,6.655,6.167,5.786,5.479,$
        5.225,5.008,4.729,4.406,4.114,3.836,3.564,3.295];
        tac= [1.775532e-03,3.239127e-02,2.774877e+00,8.931059e+00,1.391659e+01,1.773189e+01,2.070444e+01,2.346612e+01,$
        2.531194e+01,2.660934e+01,2.522953e+01,2.245561e+01,1.972853e+01,1.740301e+01,1.553417e+01,1.406894e+01,1.292658e+01,$
        1.202974e+01,1.131428e+01,1.073048e+01,1.004793e+01,9.288321e+00,8.647856e+00,8.065366e+00,7.510165e+00,6.966928e+00]

        if useref then begin
        ctt= [1.775532e-03,3.239127e-02,2.774877e+00,8.931059e+00,1.391659e+01,1.773189e+01,2.070444e+01,2.346612e+01,$
        2.531194e+01,2.660934e+01,2.522953e+01,2.245561e+01,1.972853e+01,1.740301e+01,1.553417e+01,1.406894e+01,1.292658e+01,$
        1.202974e+01,1.131428e+01,1.073048e+01,1.004793e+01,9.288321e+00,8.647856e+00,8.065366e+00,7.510165e+00,6.966928e+00]
        tac = [2.358366e-03,4.303226e-02,3.684540e+00,1.179610e+01,1.816910e+01,2.279551e+01,2.601503e+01,2.873279e+01, $
        3.032440e+01,3.030020e+01,2.726783e+01,2.387788e+01,2.098287e+01,1.863618e+01,1.677277e+01,1.530138e+01,1.413542e+01,$
        1.320196e+01,1.244281e+01,1.181292e+01,1.106495e+01,1.022359e+01,9.511543e+00,8.865074e+00,8.250660e+00,7.650896e+00]; 
        endif

        frameNr= n_elements(t0)
        tstart = 30
        tstop  = 120 
        output = double(fltarr(3)) 
        debug  = 10
        llsq_model = 0
        isweight = 1
        weights = (t1-t0)/total(t1-t0)
        ; weights = fltarr(frameNr) + 1.0
        directbp = 1 ; no division

        patlog = '/home/tsun/bin/tpcclib-master/build/bin/libmtga_idl.so'
        success = call_external(patlog, 'mrtm_idl', long(frameNr), double(t0), double(t1), double(tac), $
                        double(ctt),double(tstart),double(tstop),double(output),long(debug), $
                        long(isweight),double(weights),long(directbp)) 


        if directbp then print, 'BP '+nistring(output[1])
        if not directbp then print, 'R1 k2 k2A '+nistring(output) ; BP=R1*k2/k2A-1
  end


    'MA1': begin
    
  end

    'MA2': begin
    
  end

    'srtm' : begin     ; nonlinear model
    ;*****************************************************************************/
    ; SRTM
    ;*****************************************************************************/
    ; model := C3S K1 := 0.3 k2 := 0.15 k3 := 0 k4 := 0 k5 := 0 k6 := 0
    t0 = [0.,0.25,0.5,0.75,1.0,1.5,2.0,3.0,4.0,5.0,10.0,15.0,20.0,25.0,30.,35,40.,45.,50.,55,60,70,80,90,100,110];
    t1 = [0.25,0.5,0.75,1.0,1.5,2.0,3.0,4.0,5.0,10.0,15.0,20.0,25.0,30.,35,40.,45.,50.,55,60,70,80,90,100,110,120]; 
    ctt = [0.06415,1.879,84.8,66.15,41.86,26.86,21.36,19.87,17.15,13.54,10.89,9.305,8.147,7.294,6.655,6.167,5.786,5.479,$
    5.225,5.008,4.729,4.406,4.114,3.836,3.564,3.295];
    tac= [1.775532e-03,3.239127e-02,2.774877e+00,8.931059e+00,1.391659e+01,1.773189e+01,2.070444e+01,2.346612e+01,$
    2.531194e+01,2.660934e+01,2.522953e+01,2.245561e+01,1.972853e+01,1.740301e+01,1.553417e+01,1.406894e+01,1.292658e+01,$
    1.202974e+01,1.131428e+01,1.073048e+01,1.004793e+01,9.288321e+00,8.647856e+00,8.065366e+00,7.510165e+00,6.966928e+00]

    if useref then begin
    ; maybe a bad model?
    ctt= [1.775532e-03,3.239127e-02,2.774877e+00,8.931059e+00,1.391659e+01,1.773189e+01,2.070444e+01,2.346612e+01,$
    2.531194e+01,2.660934e+01,2.522953e+01,2.245561e+01,1.972853e+01,1.740301e+01,1.553417e+01,1.406894e+01,1.292658e+01,$
    1.202974e+01,1.131428e+01,1.073048e+01,1.004793e+01,9.288321e+00,8.647856e+00,8.065366e+00,7.510165e+00,6.966928e+00]
    tac = [2.358366e-03,4.303226e-02,3.684540e+00,1.179610e+01,1.816910e+01,2.279551e+01,2.601503e+01,2.873279e+01, $
    3.032440e+01,3.030020e+01,2.726783e+01,2.387788e+01,2.098287e+01,1.863618e+01,1.677277e+01,1.530138e+01,1.413542e+01,$
    1.320196e+01,1.244281e+01,1.181292e+01,1.106495e+01,1.022359e+01,9.511543e+00,8.865074e+00,8.250660e+00,7.650896e+00]; 
    endif

    frameNr= n_elements(t0)
    ; tstart = 30
    ; tstop  = 90 
    output = double(fltarr(12))   ; 3+3+6
    debug  = 10
    llsq_model = 0
    isweight = 1
    weights = (t1-t0)/total(t1-t0)
    ; weights = fltarr(frameNr) + 1.0
    ; directbp = 1
    def_pmin = [0.001,0.000001,0.0]
    def_pmax = [10.0,10.0,60.0]
    doSD = 1
    doCL = 0
    bootstrapIter = 200  ; has to be larger than 100!
    matrix = double(fltarr(bootstrapIter*3))    ; change with num_param

    patlog = '/home/tsun/bin/tpcclib-master/build/bin/libmtga_idl.so'
    success = call_external(patlog, 'srtm_idl', long(frameNr), double(t0), double(t1), double(tac), $
                    double(ctt),output,long(debug), $
                    long(isweight),double(weights),double(def_pmin),double(def_pmax),long(doSD),long(doCL), $
                    long(bootstrapIter),matrix) 

   print, 'R1 k2 BP '+nistring(output)                  
   stop
  end

    'FRTM' : begin
  
  end


    '1TCM' : begin
    t1 = [ 0.0200000,     0.165000,     0.335000,     0.500000,     0.665000,     0.835000, $
      1.00000,      1.16500,      1.33500,      1.50000,      1.66500,      1.83500, $
      2.00000,      2.16500,      2.33500,      2.50000,      2.66500,      2.83500, $
      3.00000,      3.16500,      3.33500,      3.99000,      7.74000,      9.50000, $
      11.5000,      13.5000,      15.5000,      17.5000,      19.5000,      21.5000, $
      23.5000,      25.5000,      27.5000,      29.5000,      31.5000,      33.5000, $
      35.5000,      37.5000,      39.5000,      41.5000,      43.5000,      45.5000, $
      47.5000,      49.5000,      51.5000,      53.5000,      55.5000,      57.5000, $
      59.5000,      61.5000,      63.5000,      65.5000,      67.5000,      69.5000, $
      71.5000,      73.5000,      75.5000,      77.5000,      79.5000,      81.5000, $
      83.5000,      85.5000,      87.5000,      89.5000,      91.5000,      93.5000, $
      95.5000,      97.5000,      99.5000,      101.500,      103.500,      105.500, $
      107.500,      109.500,      111.500,      113.500,      115.500,      117.500, $
      119.500 ]
    t0 = [[0], t1[0:n_elements(t1)-2]] 
    ctt = [  0.0269768,    0.0784673,     0.154387,      20.2208,      110.244,      69.2936, $
      54.0065,      44.0441,      38.6799,      32.5098,      27.7044,      25.0393, $
      23.9469,      21.8217,      20.9813,      22.1366,      20.3976,      20.8287, $
      20.6408,      20.2345,      19.3993,      18.7402,      12.9408,      12.0988, $
      11.2533,      10.5118,      9.86051,      9.28776,      8.78324,      8.33802, $
      7.94436,      7.59550,      7.28559,      7.00951,      6.76285,      6.54175, $
      6.34286,      6.16328,      6.00046,      5.85223,      5.71665,      5.59206, $
      5.47703,      5.37030,      5.27076,      5.17748,      5.08963,      5.00650, $
      4.92745,      4.85195,      4.77952,      4.70976,      4.64230,      4.57685, $
      4.51312,      4.45090,      4.38995,      4.33013,      4.27127,      4.21323, $
      4.15591,      4.09919,      4.04300,      3.98725,      3.93188,      3.87683, $
      3.82205,      3.76751,      3.71315,      3.65895,      3.60487,      3.55089, $
      3.49700,      3.44317,      3.38938,      3.33563,      3.28189,      3.22816, $
      3.17443 ]
    tac = [  3.23076e-05,  0.000935636,   0.00323977,     0.201574,      1.46567,      3.21735, $
      4.31376,      5.12866,      5.78688,      6.29235,      6.67452,      6.98037, $
      7.23085,      7.44186,      7.62236,      7.79484,      7.95604,      8.10354, $
      8.24434,      8.37479,      8.49232,      8.85495,      9.20914,      8.70112, $
      8.13596,      7.60048,      7.10422,      6.65097,      6.24108,      5.87285, $
      5.54347,      5.24963,      4.98786,      4.75475,      4.54707,      4.36184, $
      4.19635,      4.04818,      3.91516,      3.79538,      3.68714,      3.58896, $
      3.49955,      3.41777,      3.34262,      3.27323,      3.20887,      3.14886, $
      3.09263,      3.03970,      2.98961,      2.94200,      2.89654,      2.85294, $
      2.81096,      2.77037,      2.73100,      2.69268,      2.65526,      2.61862, $
      2.58266,      2.54728,      2.51241,      2.47796,      2.44389,      2.41013, $
      2.37664,      2.34338,      2.31032,      2.27742,      2.24466,      2.21202, $
      2.17947,      2.14700,      2.11459,      2.08223,      2.04990,      2.01761, $
      1.98533 ]

    frameNr= n_elements(t0)
    ; tstart = 30
    ; tstop  = 90 
    output = double(fltarr(12))   ; 3+3+6
    debug  = 10
    llsq_model = 0
    isweight = 1
    ; weights = fltarr(frameNr) + 1.0
    weights = (t1-t0)/total(t1-t0)
	t = (t0+t1)/2
    ; directbp = 1
    def_pmin = [0.0,0.00001,0.0]
    def_pmax = [5.0,10.0,0.0]  
    fVb = 0.0     ; fix Vb to 0 as plasma presented
    doSD = 1
    doCL = 0
    bootstrapIter = 200  ; has to be larger than 100!
    matrix = double(fltarr(bootstrapIter*4)) ; change with num_param

    patlog = '/home/tsun/bin/tpcclib-master/build/bin/libmtga_idl.so'
    success = call_external(patlog, 'tcm1_idl', long(frameNr), double(t), double(tac), $  ; double(t0), double(t1)
                    double(ctt),output,long(debug), $
                    long(isweight),double(weights),double(def_pmin),double(def_pmax), $
                    double(fVb),long(doSD),long(doCL), $
                    long(bootstrapIter),matrix) 

    ; print, 'True k1 k2' + nistring(0.12) + nistring(0.2)
  end


    '2TCM' : begin

    t0 = [0.,0.25,0.5,0.75,1.0,1.5,2.0,3.0,4.0,5.0,10.0,15.0,20.0,25.0,30.,35,40.,45.,50.,55,60,70,80,90,100,110];
    t1 = [0.25,0.5,0.75,1.0,1.5,2.0,3.0,4.0,5.0,10.0,15.0,20.0,25.0,30.,35,40.,45.,50.,55,60,70,80,90,100,110,120]; 
    ctt = [0.06415,1.879,84.8,66.15,41.86,26.86,21.36,19.87,17.15,13.54,10.89,9.305,8.147,7.294,6.655,6.167,5.786,5.479,$
    5.225,5.008,4.729,4.406,4.114,3.836,3.564,3.295];
    tac = [2.99e-4,5.463e-3,4.679e-1,1.523e0,2.443563e+00,3.236325e+00,4.015189e+00,4.897500e+00,5.656930e+00, $
    7.241067e+00,9.052315e+00,1.013719e+01,1.076269e+01,1.107761e+01,1.118487e+01,1.115529e+01,1.103736e+01,1.086400e+01,1.065737e+01,1.043218e+01, $
    1.008137e+01,9.610046e+00,9.152024e+00,8.711683e+00,8.287338e+00,7.875306e+00];

    ; if useref then begin
    ; ; maybe a bad model?
    ; ctt= [1.775532e-03,3.239127e-02,2.774877e+00,8.931059e+00,1.391659e+01,1.773189e+01,2.070444e+01,2.346612e+01,$
    ; 2.531194e+01,2.660934e+01,2.522953e+01,2.245561e+01,1.972853e+01,1.740301e+01,1.553417e+01,1.406894e+01,1.292658e+01,$
    ; 1.202974e+01,1.131428e+01,1.073048e+01,1.004793e+01,9.288321e+00,8.647856e+00,8.065366e+00,7.510165e+00,6.966928e+00]
    ; tac = [2.358366e-03,4.303226e-02,3.684540e+00,1.179610e+01,1.816910e+01,2.279551e+01,2.601503e+01,2.873279e+01, $
    ; 3.032440e+01,3.030020e+01,2.726783e+01,2.387788e+01,2.098287e+01,1.863618e+01,1.677277e+01,1.530138e+01,1.413542e+01,$
    ; 1.320196e+01,1.244281e+01,1.181292e+01,1.106495e+01,1.022359e+01,9.511543e+00,8.865074e+00,8.250660e+00,7.650896e+00]; 
    ; endif

    frameNr= n_elements(t0)
    ; tstart = 30
    ; tstop  = 90 
    output = double(fltarr(12))   ; 3+3+6
    debug  = 10
    llsq_model = 0
    isweight = 1
    ; weights = fltarr(frameNr) + 1.0
    weights = (t1-t0)/total(t1-t0)
    t = (t0+t1)/2
    ; directbp = 1
    def_pmin = [0.0,0.00001,0.0,0.0]
    def_pmax = [5.0,10.0,2.0,0.08]
    fVb = 0.0     ; fix Vb to 0 as plasma presented
    doSD = 1
    doCL = 0
    bootstrapIter = 200  ; has to be larger than 100!
    matrix = double(fltarr(bootstrapIter*4)) ; change with num_param

    patlog = '/home/tsun/bin/tpcclib-master/build/bin/libmtga_idl.so'
    success = call_external(patlog, 'tcm2_idl', long(frameNr), double(t), double(tac), $  ; double(t0), double(t1)
                    double(ctt),output,long(debug), $
                    long(isweight),double(weights),double(def_pmin),double(def_pmax), $
                    double(fVb),long(doSD),long(doCL), $
                    long(bootstrapIter),matrix) 

;  /* K1   
;  /* K1/k2 
;   /* k3    
;   /* Vb    
  stop
  end

    '2TCMrev' : begin
    t0 = [0.,0.25,0.5,0.75,1.0,1.5,2.0,3.0,4.0,5.0,10.0,15.0,20.0,25.0,30.,35,40.,45.,50.,55,60,70,80,90,100,110];
    t1 = [0.25,0.5,0.75,1.0,1.5,2.0,3.0,4.0,5.0,10.0,15.0,20.0,25.0,30.,35,40.,45.,50.,55,60,70,80,90,100,110,120]; 
    ctt = [0.06415,1.879,84.8,66.15,41.86,26.86,21.36,19.87,17.15,13.54,10.89,9.305,8.147,7.294,6.655,6.167,5.786,5.479,$
    5.225,5.008,4.729,4.406,4.114,3.836,3.564,3.295];
    tac= [1.775532e-03,3.239127e-02,2.774877e+00,8.931059e+00,1.391659e+01,1.773189e+01,2.070444e+01,2.346612e+01,$
    2.531194e+01,2.660934e+01,2.522953e+01,2.245561e+01,1.972853e+01,1.740301e+01,1.553417e+01,1.406894e+01,1.292658e+01,$
    1.202974e+01,1.131428e+01,1.073048e+01,1.004793e+01,9.288321e+00,8.647856e+00,8.065366e+00,7.510165e+00,6.966928e+00]
    frameNr= n_elements(t0)
    ; tstart = 30
    ; tstop  = 90 
    output = double(fltarr(15))   ; 3+3+6
    debug  = 10
    llsq_model = 0
    isweight = 1
    ; weights = fltarr(frameNr) + 1.0
    weights = (t1-t0)/total(t1-t0)
    t = (t0+t1)/2
    ; directbp = 1
    def_pmin = [0.0,0.00001,0.0,0.0,0.0]
    def_pmax = [5.0,10.0,2.0,10.0,0.08]
    fVb = 0.0     ; fix Vb to 0 as plasma presented
    doSD = 1
    doCL = 0
    bootstrapIter = 200  ; has to be larger than 100!
    matrix = double(fltarr(bootstrapIter*5)) ; change with num_param

    patlog = '/home/tsun/bin/tpcclib-master/build/bin/libmtga_idl.so'
    success = call_external(patlog, 'tcm2re_idl', long(frameNr), double(t), double(tac), $
                    double(ctt),output,long(debug), $
                    long(isweight),double(weights),double(def_pmin),double(def_pmax), $
                    double(fVb),long(doSD),long(doCL), $
                    long(bootstrapIter),matrix) 

;  /* K1   
;  /* K1/k2 
;   /* k3    
;   /* Vb   
  stop 
  end
    '3TCM' : begin
  
  end


  'MBF' : begin
    t1 = [ 0.0200000,     0.165000,     0.335000,     0.500000,     0.665000,     0.835000, $
      1.00000,      1.16500,      1.33500,      1.50000,      1.66500,      1.83500, $
      2.00000,      2.16500,      2.33500,      2.50000,      2.66500,      2.83500, $
      3.00000,      3.16500,      3.33500,      3.99000,      7.74000,      9.50000, $
      11.5000,      13.5000,      15.5000,      17.5000,      19.5000,      21.5000, $
      23.5000,      25.5000,      27.5000,      29.5000,      31.5000,      33.5000, $
      35.5000,      37.5000,      39.5000,      41.5000,      43.5000,      45.5000, $
      47.5000,      49.5000,      51.5000,      53.5000,      55.5000,      57.5000, $
      59.5000,      61.5000,      63.5000,      65.5000,      67.5000,      69.5000, $
      71.5000,      73.5000,      75.5000,      77.5000,      79.5000,      81.5000, $
      83.5000,      85.5000,      87.5000,      89.5000,      91.5000,      93.5000, $
      95.5000,      97.5000,      99.5000,      101.500,      103.500,      105.500, $
      107.500,      109.500,      111.500,      113.500,      115.500,      117.500, $
      119.500 ]
    t0 = [[0], t1[0:n_elements(t1)-2]] 
    ctt = [  0.0269768,    0.0784673,     0.154387,      20.2208,      110.244,      69.2936, $
      54.0065,      44.0441,      38.6799,      32.5098,      27.7044,      25.0393, $
      23.9469,      21.8217,      20.9813,      22.1366,      20.3976,      20.8287, $
      20.6408,      20.2345,      19.3993,      18.7402,      12.9408,      12.0988, $
      11.2533,      10.5118,      9.86051,      9.28776,      8.78324,      8.33802, $
      7.94436,      7.59550,      7.28559,      7.00951,      6.76285,      6.54175, $
      6.34286,      6.16328,      6.00046,      5.85223,      5.71665,      5.59206, $
      5.47703,      5.37030,      5.27076,      5.17748,      5.08963,      5.00650, $
      4.92745,      4.85195,      4.77952,      4.70976,      4.64230,      4.57685, $
      4.51312,      4.45090,      4.38995,      4.33013,      4.27127,      4.21323, $
      4.15591,      4.09919,      4.04300,      3.98725,      3.93188,      3.87683, $
      3.82205,      3.76751,      3.71315,      3.65895,      3.60487,      3.55089, $
      3.49700,      3.44317,      3.38938,      3.33563,      3.28189,      3.22816, $
      3.17443 ]
    tac = [  3.23076e-05,  0.000935636,   0.00323977,     0.201574,      1.46567,      3.21735, $
      4.31376,      5.12866,      5.78688,      6.29235,      6.67452,      6.98037, $
      7.23085,      7.44186,      7.62236,      7.79484,      7.95604,      8.10354, $
      8.24434,      8.37479,      8.49232,      8.85495,      9.20914,      8.70112, $
      8.13596,      7.60048,      7.10422,      6.65097,      6.24108,      5.87285, $
      5.54347,      5.24963,      4.98786,      4.75475,      4.54707,      4.36184, $
      4.19635,      4.04818,      3.91516,      3.79538,      3.68714,      3.58896, $
      3.49955,      3.41777,      3.34262,      3.27323,      3.20887,      3.14886, $
      3.09263,      3.03970,      2.98961,      2.94200,      2.89654,      2.85294, $
      2.81096,      2.77037,      2.73100,      2.69268,      2.65526,      2.61862, $
      2.58266,      2.54728,      2.51241,      2.47796,      2.44389,      2.41013, $
      2.37664,      2.34338,      2.31032,      2.27742,      2.24466,      2.21202, $
      2.17947,      2.14700,      2.11459,      2.08223,      2.04990,      2.01761, $
      1.98533 ]

    frameNr= n_elements(t0)
    ; tstart = 30
    ; tstop  = 90 
    output = double(fltarr(12))   ; 3+3+6
    debug  = 10
    llsq_model = 0
    isweight = 1
    ; weights = fltarr(frameNr) + 1.0
    weights = (t1-t0)/total(t1-t0)
	  t = (t0+t1)/2
    ; directbp = 1
    def_pmin = [0.0,0.00001,0.0]
    def_pmax = [5.0,10.0,1.0]  
    fVb = 0.0     ; fix Vb to 0 as plasma presented
    doSD = 1
    doCL = 0
    bootstrapIter = 200  ; has to be larger than 100!
    matrix = double(fltarr(bootstrapIter*4)) ; change with num_param

    patlog = '/home/tsun/bin/tpcclib-master/build/bin/libmtga_idl.so'
    success = call_external(patlog, 'mbf_idl', long(frameNr), double(t), double(tac), $  ; double(t0), double(t1),
                    double(ctt),output,long(debug), $
                    long(isweight),double(weights),double(def_pmin),double(def_pmax), $
                    double(fVb),long(doSD),long(doCL), $
                    long(bootstrapIter),matrix) 

    print, output
    ; K1,k2,Vb,wss,aic,sd[k1],sd[k2],sd[vb]
    ; print, 'True k1 k2' + nistring(0.12) + nistring(0.2)
  end




  'pCT' : begin
    dt       = 1    ; %s
    addskull = 1
    tstart   = 0 ; %s
    tstop    = 49 ; %s
    frameNr = (tstop-tstart)/dt+1
    t = indgen(frameNr)*dt;     ;0:1:49
    to = indgen(25)*2 +1  ;1:2:49;
    aifs = [0, 0, 0, 0, 25, 105, 220, 350, 440, 485, 430, 300, 180, 110, 104, 108, 115, 125, 115, 108, 98, 90, 98, 108, 112];
    ts = indgen(tstop*2+1)*0.5   ;0:0.1:49;
    aif = interpol(float(aifs),to,ts, /spline);
    
    cbf = 20         ; too small mtt svd make underestimated cbf
    mtt = 10       
    delay = 5.0       ; too small svd not working
    tac = double(fltarr(n_elements(ts)))
    lib = '/home/tsun/bin/tpcclib-master/build/bin/libmtga_idl.so'
    success = call_external(lib,'simpct_idl',double(ts),double(aif),long(n_elements(ts)),double(cbf),double(mtt),double(delay),tac)


    debug = 10
    isweight = 1
    weights = fltarr(n_elements(ts))+1.0
    def_pmin = [0.0,0.0,0.0]     ; cbf, mtt; cbv and ttp are calculated 
    def_pmax = [100.0,50.0,20.0]  
    doSD = 1
    doCL = 0
    bootstrapIter = 200  ; has to be larger than 100!
    matrix = double(fltarr(bootstrapIter*n_elements(def_pmin))) ; change with num_param
    output = double(fltarr(8))
    success = call_external(lib,'pCT_idl',long(n_elements(ts)),double(ts),double(tac), $
                    double(aif),output,long(debug),$
                    long(isweight),double(weights),double(def_pmin),double(def_pmax), $
                    long(doSD),long(doCL), $
                    long(bootstrapIter),matrix) 
  
    
    cbf = output[0]
    mtt = output[1]
    delay = output[2]
    cbv = cbf*mtt/60;
    ttp = where(tac eq max(tac));
    print, output,' ', cbv,' ', ttp/10.


    ; stop

    ; Apply conventional block-circulant SVD approach
    if 0 then begin
      lambda = 0.2   ; truncation
      mpad   = 2
      mask   = 0
      dt    /= 1.    ; /=10.
      first  = 0
      last   = 200     ; depends ont mtt
      ; tac = congrid(tac,10)
      ; aif = congrid(aif,10)
      pct_bsvd, tac, aif, dt, lambda, mpad, mask, $
                  cbf=cbfmap, cbv=cbvmap, mtt=mttmap,delay=delaymap,k=k,$
                  first=first,last=last

      print, cbfmap, mttmap/10., cbvmap, delaymap, ttp/10.

      stop 
    endif

  end

















    'sim_patlak': begin
    ;*****************************************************************************/
    ; SIM_PATLAK
    ;*****************************************************************************/
        t0 = [0.,0.25,0.5,0.75,1.0,1.5,2.0,3.0,4.0,5.0,10.0,15.0,20.0,25.0,30.,35,40.,45.,50.,55,60,70,80,90,100,110];
        t1 = [0.25,0.5,0.75,1.0,1.5,2.0,3.0,4.0,5.0,10.0,15.0,20.0,25.0,30.,35,40.,45.,50.,55,60,70,80,90,100,110,120]; 
        ctt = [0.06415,1.879,84.8,66.15,41.86,26.86,21.36,19.87,17.15,13.54,10.89,9.305,8.147,7.294,6.655,6.167,5.786,5.479,$
        5.225,5.008,4.729,4.406,4.114,3.836,3.564,3.295];
        tacref = [2.99e-4,5.463e-3,4.679e-1,1.523e0,2.443563e+00,3.236325e+00,4.015189e+00,4.897500e+00,5.656930e+00, $
        7.241067e+00,9.052315e+00,1.013719e+01,1.076269e+01,1.107761e+01,1.118487e+01,1.115529e+01,1.103736e+01,1.086400e+01,1.065737e+01,1.043218e+01, $
        1.008137e+01,9.610046e+00,9.152024e+00,8.711683e+00,8.287338e+00,7.875306e+00];

        frameNr= n_elements(t0)
        tstart = 30
        tstop  = 90 
        frame_used = where( (t0+t1)/2 gt tstart and (t0+t1)/2 lt tstop ) ; linear part
        output = double(fltarr(frameNr)) 
        Ki = 0.00504756 
        Vb = 1.48738
        debug = 1

        success = call_external(patlog, 'simPatlak', long(frameNr),double(Ki), double(Vb), double(t0), double(t1), $
                        double(ctt),double(tstart),double(tstop),double(output),long(debug) ) 


        print, output[frame_used]
  end

    'sim_logan': begin
    ;*****************************************************************************/
    ; SIM_LOGAN
    ;*****************************************************************************/
        t0 = [0.,0.25,0.5,0.75,1.0,1.5,2.0,3.0,4.0,5.0,10.0,15.0,20.0,25.0,30.,35,40.,45.,50.,55,60,70,80,90,100,110];
        t1 = [0.25,0.5,0.75,1.0,1.5,2.0,3.0,4.0,5.0,10.0,15.0,20.0,25.0,30.,35,40.,45.,50.,55,60,70,80,90,100,110,120]; 
        ctt = [0.06415,1.879,84.8,66.15,41.86,26.86,21.36,19.87,17.15,13.54,10.89,9.305,8.147,7.294,6.655,6.167,5.786,5.479,$
        5.225,5.008,4.729,4.406,4.114,3.836,3.564,3.295];
        tac= [1.775532e-03,3.239127e-02,2.774877e+00,8.931059e+00,1.391659e+01,1.773189e+01,2.070444e+01,2.346612e+01,$
        2.531194e+01,2.660934e+01,2.522953e+01,2.245561e+01,1.972853e+01,1.740301e+01,1.553417e+01,1.406894e+01,1.292658e+01,$
        1.202974e+01,1.131428e+01,1.073048e+01,1.004793e+01,9.288321e+00,8.647856e+00,8.065366e+00,7.510165e+00,6.966928e+00]

        if useref then begin
        ctt= [1.775532e-03,3.239127e-02,2.774877e+00,8.931059e+00,1.391659e+01,1.773189e+01,2.070444e+01,2.346612e+01,$
        2.531194e+01,2.660934e+01,2.522953e+01,2.245561e+01,1.972853e+01,1.740301e+01,1.553417e+01,1.406894e+01,1.292658e+01,$
        1.202974e+01,1.131428e+01,1.073048e+01,1.004793e+01,9.288321e+00,8.647856e+00,8.065366e+00,7.510165e+00,6.966928e+00]
        tac = [2.358366e-03,4.303226e-02,3.684540e+00,1.179610e+01,1.816910e+01,2.279551e+01,2.601503e+01,2.873279e+01, $
        3.032440e+01,3.030020e+01,2.726783e+01,2.387788e+01,2.098287e+01,1.863618e+01,1.677277e+01,1.530138e+01,1.413542e+01,$
        1.320196e+01,1.244281e+01,1.181292e+01,1.106495e+01,1.022359e+01,9.511543e+00,8.865074e+00,8.250660e+00,7.650896e+00]; 
        endif
        
        frameNr= n_elements(t0)
        tstart = 30
        tstop  = 90 
        frame_used = where( (t0+t1)/2 gt tstart and (t0+t1)/2 lt tstop ) ; linear part
        output = double(fltarr(frameNr)) 
        Dv =  2.0179289
        Ic = -28.5  ;-11 no approximate
        debug = 1
        k2 = -0.163     ; no reference

        success = call_external(patlog, 'simLogan', long(frameNr),double(Dv), double(Ic), double(t0), double(t1), $
                        double(ctt),double(tstart),double(tstop),double(output),long(debug),double(k2) ) 


        print, output[frame_used]
  end

  else: begin 
    print, 'Unsupportted model!'
  end
endcase



  print, 'elapsed time ' + nistring(systime(1)-time0) + ' seconds...'


stop


End

