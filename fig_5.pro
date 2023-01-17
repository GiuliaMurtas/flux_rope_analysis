; ------------------------------
; ## CALCULATION OF VARIABLES ##
; ------------------------------

;rdmpi,ds,datapath='../../../../../../media/giulia/Plasma/kink_instability_PIP_0',time_step=[18],var=['pr_p','pr_n','ro_p','ro_n']

;Tp = (5./3.)*ds.pr_p[*,*,441,0]/(2.*ds.ro_p[*,*,441,0])/((5./3.)*ds.pr_p[0,0,441,0]/(2.*ds.ro_p[0,0,441,0]))
;Tn = (5./3.)*ds.pr_n[*,*,441,0]/(ds.ro_n[*,*,441,0])/((5./3.)*ds.pr_n[0,0,441,0]/ds.ro_n[0,0,441,0])
;T_drift = Tn-Tp

;rdmpi,ds,datapath='../../../../../../media/giulia/Plasma/kink_instability_PIP_0',time_step=[18],var=['vx_p','vx_n','vy_p','vy_n','vz_p','vz_n']
;vd=(ds.vx_n[*,*,441,0]-ds.vx_p[*,*,441,0])^2. + (ds.vy_n[*,*,441,0]-ds.vy_p[*,*,441,0])^2. + (ds.vz_n[*,*,441,0]-ds.vz_p[*,*,441,0])^2.

;rdmpi,ds,datapath='../../../../../../media/giulia/Plasma/kink_instability_PIP_0',time_step=[18],var=['bx','by','bz']
;current_3,ds.x,ds.y,ds.z,ds.t,ds.bx,ds.by,ds.bz,jx,jy,jz
;J = sqrt((jx[*,*,441,0])^2.0 + (jy[*,*,441,0])^2.0 + (jz[*,*,441,0])^2.0)

;x=ds.x
;y=ds.y

;save,T_drift,vd,J,x,y,filename='PIP_1_mosaic.sav'

; ## Change directory below ##

;rdmpi,ds,datapath='../../../../../../media/giulia/Plasma/kink_instability_PIP_0',time_step=[18],var=['ro_p','ro_n']

popul=0
if popul eq 1 then begin
	rdmpi,ds,datapath='../kink_instability_PIP_0',time_step=[18],var=['ro_n']

	ion = ds.ion[*,*,441,0]
	rec = ds.rec[*,*,441,0]

	save,ion,rec,filename='fig_5_ion_rec.sav'
endif
if (popul eq 0) then begin
	restore,'fig_5_ion_rec.sav'
endif

; -----------------------
; ## CONTOUR MAP PIP 1 ##
; -----------------------

; ## Load the file 'PIP_1_mosaic.sav' ##
restore,'PIP_1_mosaic.sav'

levels=findgen(201)/200*0.3-0.15	; Temperature difference
levels1=findgen(201)/200*0.19+0.0	; velocity drift
levels2=findgen(201)/200*5.0+0.0	; current density magnitude
levels3=findgen(201)/200*1+0.0		; ionisation-recombination rates

w=window(dimension=[1920,1080],/buffer)

ca=contour(reform(T_drift)>levels[0]+0.0001,x,y,c_value=levels,/fill,xstyle=1,ystyle=1,xr=[-1.5 ,1.5],yr=[-1.5,1.5],rgb_table=69,aspect_ratio=1, AXIS_STYLE=2,col=0,xtitle='x',ytitle='y',xtickfont_size=20,ytickfont_size=20,font_size=20,title='$T_n$ - $T_p$',C_LABEL_SHOW=0,position=[0.05,0.55,0.23,0.9],/current)
p=plot([-0.6,-0.6],[0.2,0.9],thick=2,color='black',/overplot)
p=plot([0.3,0.3],[0.2,0.9],thick=2,color='black',/overplot)
p=plot([-0.6,0.3],[0.2,0.2],thick=2,color='black',/overplot)
p=plot([-0.6,0.3],[0.9,0.9],thick=2,color='black',/overplot)

cb=contour(reform(sqrt(vd))>levels1[0]+0.0001,x,y,c_value=levels1,/fill,xstyle=1,ystyle=1,xr=[-1.5,1.5],yr=[-1.5,1.5],rgb_table=59,aspect_ratio=1, AXIS_STYLE=2,col=0,xtitle='x',xtickfont_size=20,ytickfont_size=20,font_size=20,title='|$v_D$|',C_LABEL_SHOW=0,position=[0.23,0.55,0.41,0.9],/current)
p=plot([-0.6,-0.6],[0.2,0.9],thick=2,color='black',/overplot)
p=plot([0.3,0.3],[0.2,0.9],thick=2,color='black',/overplot)
p=plot([-0.6,0.3],[0.2,0.2],thick=2,color='black',/overplot)
p=plot([-0.6,0.3],[0.9,0.9],thick=2,color='black',/overplot)

cc=contour(reform(J)>levels2[0]+0.0001,x,y,c_value=levels2,/fill,xstyle=1,ystyle=1,xr=[-1.5,1.5],yr=[-1.5,1.5],rgb_table=56,aspect_ratio=1, AXIS_STYLE=2,col=0,xtitle='x',xtickfont_size=20,ytickfont_size=20,font_size=20,title='J',C_LABEL_SHOW=0,position=[0.41,0.55,0.59,0.9],/current)
p=plot([-0.6,-0.6],[0.2,0.9],thick=2,color='black',/overplot)
p=plot([0.3,0.3],[0.2,0.9],thick=2,color='black',/overplot)
p=plot([-0.6,0.3],[0.2,0.2],thick=2,color='black',/overplot)
p=plot([-0.6,0.3],[0.9,0.9],thick=2,color='black',/overplot)

c1=contour(reform(alog10(ion)),x,y,/fill,xstyle=1,ystyle=1,xr=[-1.5,1.5],yr=[-1.5,1.5],rgb_table=60,aspect_ratio=1, AXIS_STYLE=2,col=0,xtitle='x',xtickfont_size=20,ytickfont_size=20,font_size=20,title='$\Gamma_{ion}$',C_LABEL_SHOW=0,position=[0.59,0.55,0.77,0.9],/current)
p=plot([-0.6,-0.6],[0.2,0.9],thick=2,color='black',/overplot)
p=plot([0.3,0.3],[0.2,0.9],thick=2,color='black',/overplot)
p=plot([-0.6,0.3],[0.2,0.2],thick=2,color='black',/overplot)
p=plot([-0.6,0.3],[0.9,0.9],thick=2,color='black',/overplot)

c2=contour(reform(alog10(rec)),x,y,/fill,xstyle=1,ystyle=1,xr=[-1.5,1.5],yr=[-1.5,1.5],rgb_table=57,aspect_ratio=1, AXIS_STYLE=2,col=0,xtitle='x',xtickfont_size=20,ytickfont_size=20,font_size=20,title='$\Gamma_{rec}$',C_LABEL_SHOW=0,position=[0.77,0.55,0.95,0.9],/current)
p=plot([-0.6,-0.6],[0.2,0.9],thick=2,color='black',/overplot)
p=plot([0.3,0.3],[0.2,0.9],thick=2,color='black',/overplot)
p=plot([-0.6,0.3],[0.2,0.2],thick=2,color='black',/overplot)
p=plot([-0.6,0.3],[0.9,0.9],thick=2,color='black',/overplot)

aa=ca.axes
aa[0].tickvalue=[-1.0,-0.5,0.0,0.5,1.0]
aa[1].tickvalue=[-1.0,-0.5,0.0,0.5,1.0]

ab=cb.axes
ab[0].tickvalue=[-1.0,-0.5,0.0,0.5,1.0]
ab[1].tickvalue=[-1.0,-0.5,0.0,0.5,1.0]
ab[1].tickname=['','','','','']

ac=cc.axes
ac[0].tickvalue=[-1.0,-0.5,0.0,0.5,1.0]
ac[1].tickvalue=[-1.0,-0.5,0.0,0.5,1.0]
ac[1].tickname=['','','','','']

a1=c1.axes
a1[0].tickvalue=[-1.0,-0.5,0.0,0.5,1.0]
a1[1].tickvalue=[-1.0,-0.5,0.0,0.5,1.0]
a1[1].tickname=['','','','','']

a2=c2.axes
a2[0].tickvalue=[-1.0,-0.5,0.0,0.5,1.0]
a2[1].tickvalue=[-1.0,-0.5,0.0,0.5,1.0]
a2[1].tickname=['','','','','']

cd=contour(reform(T_drift)>levels[0]+0.0001,x,y,c_value=levels,/fill,xstyle=1,ystyle=1,xr=[-0.6,0.3],yr=[0.2,0.9],rgb_table=69,aspect_ratio=1, AXIS_STYLE=2,col=0,xtitle='x',ytitle='y',xtickfont_size=20,ytickfont_size=20,font_size=20,C_LABEL_SHOW=0,position=[0.05,0.2,0.23,0.55],/current)

ce=contour(reform(sqrt(vd))>levels1[0]+0.0001,x,y,c_value=levels1,/fill,xstyle=1,ystyle=1,xr=[-0.6,0.3],yr=[0.2,0.9],rgb_table=59,aspect_ratio=1, AXIS_STYLE=2,col=0,xtitle='x',xtickfont_size=20,ytickfont_size=20,font_size=20,C_LABEL_SHOW=0,position=[0.23,0.2,0.41,0.55],/current)

cf=contour(reform(J)>levels2[0]+0.0001,x,y,c_value=levels2,/fill,xstyle=1,ystyle=1,xr=[-0.6,0.3],yr=[0.2,0.9],rgb_table=56,aspect_ratio=1, AXIS_STYLE=2,col=0,xtitle='x',xtickfont_size=20,ytickfont_size=20,font_size=20,C_LABEL_SHOW=0,position=[0.41,0.2,0.59,0.55],/current)

c3=contour(reform(alog10(ion)),x,y,/fill,xstyle=1,ystyle=1,xr=[-0.6,0.3],yr=[0.2,0.9],rgb_table=60,aspect_ratio=1, AXIS_STYLE=2,col=0,xtitle='x',xtickfont_size=20,ytickfont_size=20,font_size=20,C_LABEL_SHOW=0,position=[0.59,0.2,0.77,0.55],/current)

c4=contour(reform(alog10(rec)),x,y,/fill,xstyle=1,ystyle=1,xr=[-0.6,0.3],yr=[0.2,0.9],rgb_table=57,aspect_ratio=1, AXIS_STYLE=2,col=0,xtitle='x',xtickfont_size=20,ytickfont_size=20,font_size=20,C_LABEL_SHOW=0,position=[0.77,0.2,0.95,0.55],/current)

ad=cd.axes
ad[0].tickvalue=[-0.4,-0.2,0.0,0.2]
ad[1].tickvalue=[0.3,0.5,0.7,0.9]

ae=ce.axes
ae[0].tickvalue=[-0.4,-0.2,0.0,0.2]
ae[1].tickvalue=[0.3,0.5,0.7,0.9]
ae[1].tickname=['','','','']

af=cf.axes
af[0].tickvalue=[-0.4,-0.2,0.0,0.2]
af[1].tickvalue=[0.3,0.5,0.7,0.9]
af[1].tickname=['','','','']

a3=c3.axes
a3[0].tickvalue=[-0.4,-0.2,0.0,0.2]
a3[1].tickvalue=[0.3,0.5,0.7,0.9]
a3[1].tickname=['','','','']

a4=c4.axes
a4[0].tickvalue=[-0.4,-0.2,0.0,0.2]
a4[1].tickvalue=[0.3,0.5,0.7,0.9]
a4[1].tickname=['','','','']

bar1=COLORBAR(TARGET=cd,rgb_table=69,ORIENTATION=0,position=[0.07,0.1,0.21,0.15],TEXTPOS = 0,major=3,taper=0,FONT_SIZE=18, range=[-0.15,0.15],TICKVALUES=[-0.1,0.0,0.1],TICKNAME=['-0.1','0.0','0.1'])
bar2=COLORBAR(TARGET=ce,rgb_table=59,ORIENTATION=0,position=[0.25,0.1,0.39,0.15],TEXTPOS = 0,major=3,taper=0,FONT_SIZE=18, range=[0,0.19],TICKVALUES=[0,0.09,0.18],TICKNAME=['0.0','0.09','0.18'])
bar3=COLORBAR(TARGET=cf,rgb_table=56,ORIENTATION=0,position=[0.43,0.1,0.57,0.15],TEXTPOS = 0,major=3,taper=0,FONT_SIZE=18, range=[0,5],TICKVALUES=[0,2.5,5],TICKNAME=['0.0','2.5','5'])
bar4=COLORBAR(TARGET=c3,rgb_table=60,ORIENTATION=0,position=[0.61,0.1,0.75,0.15],TEXTPOS = 0,major=3,taper=0,FONT_SIZE=18, range=[0,1],TICKVALUES=[0,0.5,1],TICKNAME=['0.0','0.5','1.0'])
bar5=COLORBAR(TARGET=c4,rgb_table=57,ORIENTATION=0,position=[0.79,0.1,0.93,0.15],TEXTPOS = 0,major=3,taper=0,FONT_SIZE=18, range=[0,1],TICKVALUES=[0,0.5,1],TICKNAME=['0.0','0.5','1.0'])

ca.save,'fig_5.jpg'

END

