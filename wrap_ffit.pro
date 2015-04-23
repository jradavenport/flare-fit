pro wrap_ffit,file
  
  loadct,39,/silent
  set_plot,'X'
  plotstuff,/set,/silent
  !p.font=0



  ;readcol, file, t, f, e, /silent

  f=mrdfits('phot3500.fit',9)
  t=mrdfits('time3500.fit',9)
  e = f*0.+0.02
  f = f-1.2
  
  tstart = -10                  ; by eye start time
  tstop = 750                   ; by eye stop time
  NFL = 6                       ; the # of flare components to try
  


;-----------------------------
dur = tstop - tstart
p = FFIT(t,f,e, NFL, tstart=tstart,tstop=tstop,bic=bic,num=num)

print,num,' best fit'


!p.charsize=1.2
set_plot,'ps'
device,filename='flare_models.eps',/encap,/color

xrng = [tstart-dur*0.15, tstop + dur*0.3]

clr = (findgen(NFL)+1)/float(NFL)*250
;clr = [20,80,150,212]

plot,t,f,/xsty,xtitle='!7Time - t!d0!n',ytitle='!7Flux',psym=8,xrange=xrng;,$
for n=0l,(size(p))[1] - 1 do begin
   ptmp = p[n,0:(3+n*3-1)]
   
   oplot,t,aflare(t,(ptmp)),color=clr[n],thick=5

endfor

device,/close

;-- component plots
device,filename='flare_comps.eps',/encap,/color
n = num
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
