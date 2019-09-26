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

    useref = 0;   use reference region as input?
    testmode = 'patlak';    'logan','fur','mrtm','sim_patlak','sim_logan','SRTM','2TCMrev'


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
        output = double(fltarr(frame_used+1)) 
        debug  = 10
        fur_mode = 0

        patlog = '/home/tsun/bin/tpcclib-master/build/bin/libmtga_idl.so'
        success = call_external(patlog, 'regfur_idl', long(frameNr), double(t0), double(t1), double(tac), $
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
    ; directbp = 1
    def_pmin = [0.0,0.00001,0.0,0.0]
    def_pmax = [5.0,10.0,2.0,0.08]
    fVb = 0.0     ; fix Vb to 0 as plasma presented
    doSD = 1
    doCL = 0
    bootstrapIter = 200  ; has to be larger than 100!
    matrix = double(fltarr(bootstrapIter*4)) ; change with num_param

    patlog = '/home/tsun/bin/tpcclib-master/build/bin/libmtga_idl.so'
    success = call_external(patlog, 'tcm2_idl', long(frameNr), double(t0), double(t1), double(tac), $
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
    ; directbp = 1
    def_pmin = [0.0,0.00001,0.0,0.0,0.0]
    def_pmax = [5.0,10.0,2.0,10.0,0.08]
    fVb = 0.0     ; fix Vb to 0 as plasma presented
    doSD = 1
    doCL = 0
    bootstrapIter = 200  ; has to be larger than 100!
    matrix = double(fltarr(bootstrapIter*5)) ; change with num_param

    patlog = '/home/tsun/bin/tpcclib-master/build/bin/libmtga_idl.so'
    success = call_external(patlog, 'tcm2re_idl', long(frameNr), double(t0), double(t1), double(tac), $
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







stop


End