
Function rsq, y,yfit
    ; %Calculates the coefficient of determination (r-squared)
    yresid = y - yfit;
    SSresid = total(yresid^2);
    SStotal = (n_elements(y)-1) * variance(y);
    r2 = 1 - SSresid/SStotal;
End

Function pct_timedelay,aif,delay
    ; delay of the AIF and its sum
    aif = shift(aif,delay)     
    aif[0:delay-1] = 0
End

Pro  pct_bbbp, ctmap, cbvmap,delaymap=delaymap, AIF, dt, rho, first, last, mask, bbbmap, x, ymap, R
; %PCT_BBBP Computes a permeability map (BBBP).
; %
; %   This algorithm calculates blood-brain barrier (BBB) permeability
; %   using the Patlak model.
; %
; %   USAGE: [BBBP X YMAP R] = PCT_BBBP(CTMAP,CBVMAP,AIF,DT,RHO,FIRST,LAST,MASK);
; %
; %   PRE:
; %       CTMAP   - Untruncated, preprocessed CT sequence in HU [T x Y x X]
; %       CBVMAP  - A CBV (Cerebral Blood Volume) map in mL/100g [Y x X]
; %       AIF     - A preprocessed Arterial Input Function [T x 1]
; %       DT      - Time interval between samples in seconds [Scalar]
; %       RHO     - Average brain tissue density in g/mL [Scalar]
; %       FIRST   - The first frame to include in the calculations [Scalar]
; %       LAST    - The last frame to include in the calculations [Scalar]
; %       MASK    - Optional parameter. A logical mask [Y x X] that indicates
; %                 which pixels are to be processed (processes those that are
; %                 TRUE).
; %
; %   POST:
; %       BBBP    - A map of brain permeability in mL/100g/min [Y x X]
; %       X       - The independent variable of the patlak plot [T x 1]
; %       YMAP    - A map of dependent variables of the patlak plot [T x Y x X]
; %       R       - A map of coefficients of determination (R^2)
; %
; %   Kolbeinn Karlsson 06/05/12
; %   Advanced Multimedia Processing (AMP) Lab, Cornell University

; %Get mapsize
sizet = (size(ctmap))[1:4];
height = sizet[0]
width  = sizet[1]
depth  = sizet[2]
len    = sizet[3]

if n_elements(mask) lt 2 then $
    mask = boolarr(height, width, depth)+byte(1);

if last eq -1 then $
    last = len;

; %Pre-allocate the output variable
bbbmap = fltarr(height,width,depth);      
R      = fltarr(height,width,depth); 
ymap   = fltarr(height,width,depth,len); 


; %Calculate the AIF integral for each time step
aif_sum = fltarr(n_elements(AIF));
tmp = indgen(n_elements(AIF))
for n = 0,n_elements(AIF)-1 do begin
    aif_sum[n] = dt*int_tabulated(tmp,AIF[0:n_elements(AIF)-1]);
    ; %aif_sum(n) = sum(AIF(1:n));
endfor

; %Preconditioning the AIF (We don't want to divide by zero)
AIF[where(AIF lt 1)] = 1;

aif_sum_orig = aif_sum
AIF_orig = AIF
; ; %Calculate x, the independent variable of the Patlak plot
; x = aif_sum / AIF;

; %Calculate the y map, the dependent variables of the Patlak plot
for p = 0,depth-1 do begin
    print, 'Processing Plane '+ nistring(p)+' ..'
    for i = 0,height-1 do begin
        for j = 0,width-1 do begin
            ; %Pass through mask
            if mask[i,j,p] eq 1 then begin

                if keyword_set(delaymap) then begin
                    delay = ttpmap(i,j)*1/dt-ttp_aif;
                    if delay lt 0 then delay = 0;
                    AIF =  pct_timedelay(AIF_orig,delay);
                    aif_sum = pct_timedelay(AIF_sum_orig,delay);
                endif else begin
                    delay = 0 
                    AIF = AIF_orig
                    aif_sum = aif_sum_orig
                endelse
            
                ; %Calculate x, the independent variable of the Patlak plot
                x = aif_sum / AIF;

                ; %Calculate the y-map
                RIF = reform(ctmap[i,j,p,*]); 
                if i eq 128 and j eq 128 and p eq 100 then stop
                ymap[i,j,p,*] = RIF / AIF * 1/rho - cbvmap[i,j,p]/100;
                ; %Calculate the permeability
                xx = x[(first+delay):last];
                yy = reform(ymap[i,j,p,(first+delay):last]);
                bbbmap[i,j,p] = mean(yy/xx)    ;la_least_squares((xx,yy)    ;yy \ xx;
                ; %Calculate the coefficient of determination
                R[i,j,p] = rsq(yy,bbbmap[i,j,p]*xx)
            endif else begin
                ; %Set entries to zero
                ymap[i,j,p,*] = 0;
                R[i,j,p] = 0;
                bbbmap[i,j,p] = 0;
            endelse
        endfor
    endfor
endfor

stop

; %Convert the permeability map from mL/g/s to mL/100g/min
bbbmap = bbbmap * 60. * 100;

; %Filter negative values
bbbmap[where(bbbmap lt 0)]=0;

; ; % Remove nonbrain tissue
; bbbmap[where(mask ne 0)] = 0;

End



