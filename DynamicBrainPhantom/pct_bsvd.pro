Pro pct_bsvd, inmap,aif,dt,lambda,m,mask, $
              cbf=cbfmap,cbv=cbvmap,mtt=mttmap,delay=delaymap, k=k, $
              first=first,last=last
; %PCT_BSVD Deconvolution of the Indicator-Dilution equation using bSVD
; % by circular deconvolution and circulating back the last values
; %
; %   USAGE:  K = PCT_BSVD(INMAP,AIF,DT,LAMBDA,M,MASK)
; %
; %   INPUT:
; %       INMAP   - A [T x Y x X] matrix, so 2D?
; %       AIF     - The Arterial Input Function [T x 1]
; %       DT      - The sampling interval in seconds [Scalar]
; %       LAMBDA  - Truncation parameter for the SVD. It denotes the fraction of
; %                 the lowest singular values that are set to zero. [Scalar]
; %       M       - How often to extend the input. If input is of length T, then
; %                 it will be zeropadded to length MT. [Scalar]
; %       MASK    - A logical [Y x X] mask. The computation will only be performed
; %                 for the voxels that are logical 1.
; %
; %   OUTPUT:
; %       K       - BF * R. The impulse residue function scaled by the cerebral
; %                 blood flow. [T x Y x X].
; %
; %   This function solves the Indicator-Dilution equation
; %
; %       C = F * conv( C_a,R )
; %
; %   using block-circulant SVD. See "Tracer Arrival Timing-Insensitive Technique
; %   for Estimating Flow in MR Perfusion-Weighted Imaging Using Singular Value
; %   Decomposition With a Block-Circulant Deconvolution Matrix" by Ona Wu,
; %   Leif stergaard, Robert M. Weisskoff, Thomas Benner, Bruce R. Rosen, and
; %   A. Gregory Sorensen for more detail.
; %
; %   Kolbeinn Karlsson, 08/06/12
; %   Ruogu Fang, 06/19/2014 add circulating back the RIF
; %   Advanced Multimedia Processing (AMP) Lab, Cornell University

; Delay in AIF
; For the delay correction of R, I found that by zero padding the
; signals, Ca and C, we could get the right R even AIF if delay by t_d.
; However the computed R will be circularly delayed. For instance, if
; t_d = 1 sec, R should be

; R = [1 0.5 0 0]

; when there is no delay.

; But when there is delay, R_d will be

; R_d = [0.5 0 0 1]

; This is the result output by the MATLAB code bcsvd_long.m

; In that case, we should circularly shift R_d right by 1 sec to make it
; R_shifted = [1 0.5 0 0].



; Temporal interpolation

if (size(inmap))[0] eq 1 then begin
    print, "Input: 1D voxel supported"  
    nframe = (size(inmap))[1]
    height = 1
    width  = 1
endif else if (size(inmap))[0] eq 2 then begin
   print, "Input: 2D image supported"  
    nframe = (size(inmap))[1]
    height = (size(inmap))[2]
    width  = (size(inmap))[3]
endif else begin
   print, "Input: 3D image not supported"  
endelse

; Initialization
delaymap = fltarr(height,width)
cbfmap   = delaymap*0
cbvmap   = delaymap*0
mttmap   = delaymap*0

if n_elements(mask) lt 2 then $
    mask = fltarr(height,width)+1.0 

; Extend AIF by zero-padding
N  = m*nframe
Ca = fltarr(N)
Ca[0:nframe-1] = aif

; Create block-circulant matrix
D = fltarr(N,N)
for i=0,N-1 do $
    D[i,*] = shift(Ca,i)
D = transpose(D)

; Get inverse D
; stop
LA_SVD, D, S, U, V
; SVDC, D, S, U, V
maxS = max(diag_matrix(S))
S[where(S lt lambda*maxS)] = 0
invD = transpose(V) # pinv(diag_matrix(S)) # U ;invD = V * pinv(diag_matrix(S)) * transpose(U)

k = fltarr(N,height,width)

; Deconvolution
for i=0,height-1 do begin
  for j=0,width-1 do begin
    if mask[i,j] gt 0 then begin
      R=invD ## [inmap[*,i,j], fltarr(N-nframe)]   ;; #
      if R[n_elements(R)-1] gt R[0] then begin
         idx = where(R eq max(R))
         delay = N - idx + 1                    ;; see delay.txt
         R = shift(R,delay)
         delaymap[i,j] = delay          ; store AIF delay 
      endif else delay = 0.0
    endif
    k[*,i,j] = R
  endfor
endfor

k = k[0:nframe-1,*,*]

; Correct sampling rate
Rmap = k/dt
if not keyword_set(first) then first = 0
if not keyword_set(last) then last = nframe-1
Rmap1 = Rmap[first:last,*,*];

; stop

; Get CBF from Residuals
rho = 1.0   ; brain tissue density in g/ml

Rcol    = height   ;(size(R))[2]
Rrow    = width    ;(size(R))[3]
mask = fltarr(Rcol,Rrow) + 1.0

scaling_factor = 60*100/rho;
cbfmap = scaling_factor * reform(max(Rmap1))
cbfmap = cbfmap > 0.0
; cbfmap[where(cbfmap lt 0)] = 0.0;
; cbfmap[where(mask eq 0)] = 0.0;


; Get CBV
scaling_factor = 100/(rho);
cbvmap = scaling_factor * reform(total(Rmap1,1));
cbvmap = cbvmap > 0.0
; cbvmap[where(cbvmap lt 0)] = 0.0;
; cbvmap[where(mask eq 0)] = 0.0;



; Get MTT
; first = 0
; last  = 200 ;;nframe-1
eps = 1e-7
; mttmap = cbvmap / (cbfmap+eps) * 60
mttmap = reform(total(Rmap1,1)) / (reform(max(Rmap1))+eps);
mttmap = mttmap > 0.0
; mttmap[where(mttmap lt 0)] = 0.0;
; mttmap[where(mask eq 0)] = 0.0;



stop



End