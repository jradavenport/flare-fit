function ffit,t,f,e,ncomp,tstart = tstart, tstop = tstop,bic=bic, num=num

;--- generalize the iterative flare fit stuff.
; - dont include the flare model itself, just the subtract and mpfit
; - output should be some grid of fits

if not keyword_set(tstart) then tstart = min(t)
if not keyword_set(tstop) then tstop = max(t)

rng = where(t ge tstart and t le tstop)

fmax = double(max(f[rng],b))
tmax = double(t[rng[b]])
tdur = double(max(t[rng]) - min(t[rng]))

yy = f[rng]
xx = t[rng]
ee = e[rng]



BIC_out = fltarr(ncomp)         ;-- BIC values for all 7 fits
param_out = fltarr(3.*ncomp)    ;-- final params, all 7 fits
nflare = indgen(ncomp)+1        ;-- # of flare components
bic = fltarr(n_elements(nflare))

ptmp = [tmax, tdur/15d0, fmax] ; first big flare
param_tmp = dblarr(ncomp, 3.*ncomp)*0.-1.
ostatus = fltarr(n_elements(nflare))

nsmo = 2 ; floor(float(n_elements(yy))/20. + 1)

for i=0,ncomp-1 do begin
   if i gt 0 then begin
      ; subtract previous model
      ydif = smooth(yy - real(aflare(xx,pout_i),0),nsmo,/edge)

      ; find new peak (residual)
      fmax_i = double(max(ydif,bi))
      tmax_i = double(xx[bi])
      ; add another flare
      maxtmp = ydif[bi]
      if maxtmp lt mean(ee) then maxtmp = abs(median(yy))

      durtmp = pout_i[1] * abs(maxtmp/pout_i[2]) ; scale width
      if durtmp lt 1.39d-3 then durtmp = 3d-3

      ptmp = [ptmp, $
              tmax_i, durtmp, maxtmp]

      ; shrink est. of previous flares
      ;; for k = 0,i-1 do $
      ;;    ptmp[[k,k+1,k+2]] = ptmp[[k,k+1,k+2]]*[1., 0.6, 0.8]
   endif


   ;-- set ranges for fitting
   parinfo = replicate({limits:[0d0,0d0], limited:[1,1]}, n_elements(ptmp))
   for j=0,nflare[i]-1 do begin
      parinfo[0+j*3].limits=[min(xx)-0.01d0, max(xx)]
      parinfo[1+j*3].limits=[1.39d-3, tdur*0.1d0]           ; width lim
      parinfo[2+j*3].limits=[abs(mean(ee)*2d0), max(yy)*1.1] ; peak lim
   endfor

   ;-- minimize!
   pout_i = mpfitfun('aflare', xx, yy, ee, ptmp, $ ;-- new analytic fxn
                   /nan,perror=outerr,bestnorm=bestnorm,status=status,$
                   /double,parinfo=parinfo,maxiter=1000,/quiet)
   
   ;; plot,t,f
   ;; oplot,t,aflare(t,ptmp),color=90
   ;; oplot,t,aflare(t,pout_i),color=160


   ; update input params
   ptmp = pout_i
   param_tmp[ i, 0:(n_elements(pout_i)-1)] = pout_i

   ostatus[i] = status ;-- mpfit output status


   ; BIC == chisq + k*ln(n), k=#DOF, n=#data
   ; K free parameters (per component)
   kfp = 3. + 8. ;-- 3 for fit params, 8 for model
   bic[i] = total(((yy-aflare(xx,pout_i))/ee)^2.0) + $         ; chisq
            float((kfp*nflare[i])*alog(float(n_elements(yy)))) ; k*ln(n)


   ;-- if MPFIT didn't converge well, then increase BIC
   if ostatus[i] le 0 then bic[i] = bic[i]*2.
   if ostatus[i] gt 3 then bic[i] = bic[i]*2.

endfor



if ncomp gt 1 then begin
;-- pick the preferred # of components for the event
   THENUM = 1
; find where BIC no longer decreases by 10% of the previous value
;; dBIC = BIC[1:*]/BIC lt 0.9
;; dBIC = BIC[1:*]-BIC lt -1d-1*BIC[0]

;; turn = where(dBIC eq 0)
;; if turn[0] eq 0 then THENUM = 1
;; if turn[0] eq -1 then THENUM = ncomp
;; if turn[0] ne -1 and turn[0] ne 0 then THENUM = turn[0]+1


;-- try different way...
   dBIC = BIC/BIC[1:*]
   bett = where(dBIC gt 1)
   if bett[0] ne -1 then THENUM = max(bett)+1
   
   num = THENUM
endif

if ncomp eq 1 then begin
   num = 1
endif
   
return,param_tmp ;(param_tmp[THENUM-1,*])[*]
end
