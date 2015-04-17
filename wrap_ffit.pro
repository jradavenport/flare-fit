pro wrap_ffit,file
  
  loadct,39,/silent
  set_plot,'X'
  plotstuff,/set,/silent
  !p.font=0



  readcol, file, t, f, e, /silent

  tstart = 1.7d4                ; by eye start time
  tstop = 1.78d4                ; by eye stop time
  NFL = 3                       ; the # of flare components to try



;-----------------------------
dur = tstop - tstart
p = FFIT(t,f,e, NFL, tstart=tstart,tstop=tstop,bic=bic,num=num)

print,num,' best fit'


!p.charsize=1.2
set_plot,'ps'
device,filename='flare_fit.eps',/encap,/color

xrng = [tstart-dur*0.05, tstop + dur*0.3]

clr = (findgen(NFL)+1)/float(NFL)*250
;clr = [20,80,150,212]

plot,t,f,/xsty,xtitle='!7Time - t!d0!n',ytitle='!7Flux',psym=8,xrange=xrng;,$
for n=0l,(size(p))[1] - 1 do begin
   ptmp = p[n,0:(3+n*3-1)]
   
   oplot,t,aflare(t,(ptmp)),color=clr[n],thick=5

endfor

device,/close

;-- component plots
device,filename='flare_3comp.eps',/encap,/color
n = 2
plot,t-tstart,f,/xsty,xtitle='!7t - t!dstart!n (seconds)',ytitle='!7Flux',$
     psym=8,xrange=xrng-tstart

for i=0,n do $
   oplot, t-tstart, aflare(t, p[n, [0,1,2]+3*i]), color = clr[i],thick=5,psym=10

loadct,0,/silent
oplot, t-tstart, aflare(t, p[ n, 0:(3*(n+1)-1) ] ), color=150,thick=3,psym=10
device,/close
loadct,39,/silent

;--
device,filename='flare_2comp.eps',/encap,/color
n = 1
plot,t-tstart,f,/xsty,xtitle='!7t - t!dstart!n (seconds)',ytitle='!7Flux',$
     psym=8,xrange=xrng-tstart

for i=0,n do $
   oplot, t-tstart, aflare(t, p[n, [0,1,2]+3*i]), color = clr[i],thick=5,psym=10

loadct,0,/silent
oplot, t-tstart, aflare(t, p[ n, 0:(3*(n+1)-1) ] ), color=150,thick=3,psym=10
device,/close
loadct,39,/silent


;--
device,filename='flare_1comp.eps',/encap,/color
n = 0
plot,t-tstart,f,/xsty,xtitle='!7t - t!dstart!n (seconds)',ytitle='!7Flux',$
     psym=8,xrange=xrng-tstart

for i=0,n do $
   oplot, t-tstart, aflare(t, p[n, [0,1,2]+3*i]), color = clr[i],thick=5,psym=10

loadct,0,/silent
oplot, t-tstart, aflare(t, p[ n, 0:(3*(n+1)-1) ] ), color=150,thick=3,psym=10
device,/close
loadct,39,/silent


set_plot,'X'

stop
end
