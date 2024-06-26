; --------------------------------
; ## ION/REC FRICTIONAL HEATING ##
; --------------------------------

;Does the data need populating? yes=1, no=0
popul=0

if popul eq 1 then begin

;rdmpi,ds,datapath='../kink_instability_PIP_0',time_step=[9,11,13],var=['ro_p','ro_n','vx_p','vy_p','vz_p','vx_n','vy_n','vz_n']
rdmpi,ds,datapath='../kink_instability_PIP_0',time_step=[9],var=['ro_p','ro_n','vx_p','vy_p','vz_p','vx_n','vy_n','vz_n']
vn2_0 = ds.vx_n[*,*,441]^2.0 + ds.vy_n[*,*,441]^2.0 + ds.vz_n[*,*,441]^2.0
vp2_0 = ds.vx_p[*,*,441]^2.0 + ds.vy_p[*,*,441]^2.0 + ds.vz_p[*,*,441]^2.0
vnvp_0 = (ds.vx_n[*,*,441]*ds.vx_p[*,*,441]) + (ds.vy_n[*,*,441]*ds.vy_p[*,*,441]) + (ds.vz_n[*,*,441]*ds.vz_p[*,*,441])
fh2_0 = ds.rec[*,*,441]*ds.ro_p[*,*,441]*vp2_0 - (ds.rec[*,*,441]*ds.ro_p[*,*,441] + ds.ion[*,*,441]*ds.ro_n[*,*,441])*vnvp_0 + ds.ion[*,*,441]*ds.ro_n[*,*,441]*vn2_0

print,ds.t[0]

;delvar,ds

rdmpi,ds,datapath='../kink_instability_PIP_0',time_step=[11],var=['ro_p','ro_n','vx_p','vy_p','vz_p','vx_n','vy_n','vz_n']

vn2_1 = ds.vx_n[*,*,441]^2.0 + ds.vy_n[*,*,441]^2.0 + ds.vz_n[*,*,441]^2.0
vp2_1 = ds.vx_p[*,*,441]^2.0 + ds.vy_p[*,*,441]^2.0 + ds.vz_p[*,*,441]^2.0
vnvp_1 = (ds.vx_n[*,*,441]*ds.vx_p[*,*,441]) + (ds.vy_n[*,*,441]*ds.vy_p[*,*,441]) + (ds.vz_n[*,*,441]*ds.vz_p[*,*,441])
fh2_1 = ds.rec[*,*,441]*ds.ro_p[*,*,441]*vp2_1 - (ds.rec[*,*,441]*ds.ro_p[*,*,441] + ds.ion[*,*,441]*ds.ro_n[*,*,441])*vnvp_1 + ds.ion[*,*,441]*ds.ro_n[*,*,441]*vn2_1

print,ds.t[0]

;delvar,ds

rdmpi,ds,datapath='../kink_instability_PIP_0',time_step=[13],var=['ro_p','ro_n','vx_p','vy_p','vz_p','vx_n','vy_n','vz_n']

vn2_2 = ds.vx_n[*,*,441]^2.0 + ds.vy_n[*,*,441]^2.0 + ds.vz_n[*,*,441]^2.0
vp2_2 = ds.vx_p[*,*,441]^2.0 + ds.vy_p[*,*,441]^2.0 + ds.vz_p[*,*,441]^2.0
vnvp_2 = (ds.vx_n[*,*,441]*ds.vx_p[*,*,441]) + (ds.vy_n[*,*,441]*ds.vy_p[*,*,441]) + (ds.vz_n[*,*,441]*ds.vz_p[*,*,441])
fh2_2 = ds.rec[*,*,441]*ds.ro_p[*,*,441]*vp2_2 - (ds.rec[*,*,441]*ds.ro_p[*,*,441] + ds.ion[*,*,441]*ds.ro_n[*,*,441])*vnvp_2 + ds.ion[*,*,441]*ds.ro_n[*,*,441]*vn2_2

print,ds.t[0]

;delvar,ds

save, vn2_0,vp2_0,vnvp_0,fh2_0, vn2_1,vp2_1,vnvp_1,fh2_1, vn2_2,vp2_2,vnvp_2,fh2_2,filename='fig_6_heatterms.sav'

endif

if popul eq 0 then begin
restore,'fig_6_heatterms.sav'
endif

;print,ds.t[0], ds.t[1], ds.t[2]

;vn2_0 = ds.vx_n[*,*,441,0]^2.0 + ds.vy_n[*,*,441,0]^2.0 + ds.vz_n[*,*,441,0]^2.0
;vp2_0 = ds.vx_p[*,*,441,0]^2.0 + ds.vy_p[*,*,441,0]^2.0 + ds.vz_p[*,*,441,0]^2.0
;vnvp_0 = (ds.vx_n[*,*,441,0]*ds.vx_p[*,*,441,0]) + (ds.vy_n[*,*,441,0]*ds.vy_p[*,*,441,0]) + (ds.vz_n[*,*,441,0]*ds.vz_p[*,*,441,0])
;fh2_0 = ds.rec[*,*,441,0]*ds.ro_p[*,*,441,0]*vp2_0 - (ds.rec[*,*,441,0]*ds.ro_p[*,*,441,0] + ds.ion[*,*,441,0]*ds.ro_n[*,*,441,0])*vnvp_0 + ds.ion[*,*,441,0]*ds.ro_n[*,*,441,0]*vn2_0

;vn2_1 = ds.vx_n[*,*,441,1]^2.0 + ds.vy_n[*,*,441,1]^2.0 + ds.vz_n[*,*,441,1]^2.0
;vp2_1 = ds.vx_p[*,*,441,1]^2.0 + ds.vy_p[*,*,441,1]^2.0 + ds.vz_p[*,*,441,1]^2.0
;vnvp_1 = (ds.vx_n[*,*,441,1]*ds.vx_p[*,*,441,1]) + (ds.vy_n[*,*,441,1]*ds.vy_p[*,*,441,1]) + (ds.vz_n[*,*,441,1]*ds.vz_p[*,*,441,1])
;fh2_1 = ds.rec[*,*,441,1]*ds.ro_p[*,*,441,1]*vp2_1 - (ds.rec[*,*,441,1]*ds.ro_p[*,*,441,1] + ds.ion[*,*,441,1]*ds.ro_n[*,*,441,1])*vnvp_1 + ds.ion[*,*,441,1]*ds.ro_n[*,*,441,1]*vn2_1

;vn2_2 = ds.vx_n[*,*,441,2]^2.0 + ds.vy_n[*,*,441,2]^2.0 + ds.vz_n[*,*,441,2]^2.0
;vp2_2 = ds.vx_p[*,*,441,2]^2.0 + ds.vy_p[*,*,441,2]^2.0 + ds.vz_p[*,*,441,2]^2.0
;vnvp_2 = (ds.vx_n[*,*,441,2]*ds.vx_p[*,*,441,2]) + (ds.vy_n[*,*,441,2]*ds.vy_p[*,*,441,2]) + (ds.vz_n[*,*,441,2]*ds.vz_p[*,*,441,2])
;fh2_2 = ds.rec[*,*,441,2]*ds.ro_p[*,*,441,2]*vp2_2 - (ds.rec[*,*,441,2]*ds.ro_p[*,*,441,2] + ds.ion[*,*,441,2]*ds.ro_n[*,*,441,2])*vnvp_2 + ds.ion[*,*,441,2]*ds.ro_n[*,*,441,2]*vn2_2

; ----------------------------------------
; ## HEATING TERMS AND TEMPERATURE MAPS ##
; ----------------------------------------

; Files to be loaded: T_PIP_1.sav, fh_PIP_1.sav, oh_PIP_1.sav, vd_PIP_1.sav
restore,'T_PIP_1.sav'
restore,'fh_PIP_1.sav'
restore,'oh_PIP_1.sav'
restore,'vd_PIP_1.sav'

levels1=findgen(201)/200*0.001+0.0	; color scale for ion/rec frictional heating
levels2=findgen(201)/200*0.01+0.0	; color scale for collisional frictional heating
levels3=findgen(201)/200*0.2+0.0	; color scale for Ohmic heating
levels4=findgen(201)/200*0.08+0.07	; color scale for temperatures

print,'Max Tn = ',max(Tn_1_2)
print,'Max Tp = ',max(Tp_1_2)
print,'Max FH2 = ',max(fh2_2)

; ## Neutral temperature panels ##

w=window(dimension=[1980,1080],/buffer)

ca=contour(reform(Tn_1_0)>levels4[0]+0.0001,x,y,c_value=levels4,/fill,xstyle=1,ystyle=1,xr=[-1.,1.],yr=[-1.,1.],rgb_table=52,aspect_ratio=1, AXIS_STYLE=2,col=0,xtitle='x',ytitle='y',xtickfont_size=18,ytickfont_size=18,C_LABEL_SHOW=0,position=[0.05,0.65,0.275,0.95],/current)

cb=contour(reform(Tn_1_1)>levels4[0]+0.0001,x,y,c_value=levels4,/fill,xstyle=1,ystyle=1,xr=[-1.,1.],yr=[-1.,1.],rgb_table=52,aspect_ratio=1, AXIS_STYLE=2,col=0,xtitle='x',ytitle='y',xtickfont_size=18,ytickfont_size=18,C_LABEL_SHOW=0,position=[0.05,0.35,0.275,0.65],/current)

cc=contour(reform(Tn_1_2)>levels4[0]+0.0001,x,y,c_value=levels4,/fill,xstyle=1,ystyle=1,xr=[-1.,1.],yr=[-1.,1.],rgb_table=52,aspect_ratio=1, AXIS_STYLE=2,col=0,xtitle='x',ytitle='y',xtickfont_size=18,ytickfont_size=18,C_LABEL_SHOW=0,position=[0.05,0.05,0.275,0.35],/current)

bar=COLORBAR(TARGET=ca,rgb_table=52,ORIENTATION=0,position=[0.07,0.96,0.255,0.98],TEXTPOS = 1,major=2,taper=0,FONT_SIZE=18,TITLE='$T_n$',range=[0.07,0.15],TICKVALUES=[0.07,0.1],TICKNAME=['0.07','0.1'])

aa=ca.axes
aa[0].tickvalue=[-0.5,-0,0.5]
aa[1].tickvalue=[-0.5,-0,0.5]

ab=cb.axes
ab[0].tickvalue=[-0.5,-0,0.5]
ab[1].tickvalue=[-0.5,-0,0.5]

ac=cc.axes
ac[0].tickvalue=[-0.5,-0,0.5]
ac[1].tickvalue=[-0.5,-0,0.5]

; ## Plasma temperature panels ##

cd=contour(reform(Tp_1_0)>levels4[0]+0.0001,x,y,c_value=levels4,/fill,xstyle=1,ystyle=1,xr=[-1.,1.],yr=[-1.,1.],rgb_table=55,aspect_ratio=1, AXIS_STYLE=2,col=0,xtitle='x',xtickfont_size=18,ytickfont_size=18,C_LABEL_SHOW=0,position=[0.275,0.65,0.5,0.95],/current)

ce=contour(reform(Tp_1_1)>levels4[0]+0.0001,x,y,c_value=levels4,/fill,xstyle=1,ystyle=1,xr=[-1.,1.],yr=[-1.,1.],rgb_table=55,aspect_ratio=1, AXIS_STYLE=2,col=0,xtitle='x',xtickfont_size=18,ytickfont_size=18,C_LABEL_SHOW=0,position=[0.275,0.35,0.5,0.65],/current)

cf=contour(reform(Tp_1_2)>levels4[0]+0.0001,x,y,c_value=levels4,/fill,xstyle=1,ystyle=1,xr=[-1.,1.],yr=[-1.,1.],rgb_table=55,aspect_ratio=1, AXIS_STYLE=2,col=0,xtitle='x',xtickfont_size=18,ytickfont_size=18,C_LABEL_SHOW=0,position=[0.275,0.05,0.5,0.35],/current)

bar2=COLORBAR(TARGET=cd,rgb_table=55,ORIENTATION=0,position=[0.295,0.96,0.48,0.98],TEXTPOS = 1,major=2,taper=0,FONT_SIZE=18,TITLE='$T_p$',range=[0.007,0.15],TICKVALUES=[0.07,0.1],TICKNAME=['0.07','0.1'])

ad=cd.axes
ad[0].tickvalue=[-0.5,-0,0.5]
ad[1].tickvalue=[-0.5,-0,0.5]
ad[0].tickname=['','','']
ad[1].tickname=['','','']

ae=ce.axes
ae[0].tickvalue=[-0.5,-0,0.5]
ae[1].tickvalue=[-0.5,-0,0.5]
ae[0].tickname=['','','']
ae[1].tickname=['','','']

af=cf.axes
af[0].tickvalue=[-0.5,-0,0.5]
af[1].tickvalue=[-0.5,-0,0.5]
af[1].tickname=['','','']

; ## Collisional frictional heating panels (FH1) ##

cg=contour(reform(fh_1_0)>levels2[0]+0.0001,x,y,c_value=levels2,/fill,xstyle=1,ystyle=1,xr=[-1.,1.],yr=[-1.,1.],rgb_table=57,aspect_ratio=1, AXIS_STYLE=2,col=0,xtitle='x',xtickfont_size=18,ytickfont_size=18,C_LABEL_SHOW=0,position=[0.5,0.65,0.725,0.95],/current)

ch=contour(reform(fh_1_1)>levels2[0]+0.0001,x,y,c_value=levels2,/fill,xstyle=1,ystyle=1,xr=[-1.,1.],yr=[-1.,1.],rgb_table=57,aspect_ratio=1, AXIS_STYLE=2,col=0,xtitle='x',xtickfont_size=18,ytickfont_size=18,C_LABEL_SHOW=0,position=[0.5,0.35,0.725,0.65],/current)

ci=contour(reform(fh_1_2)>levels2[0]+0.0001,x,y,c_value=levels2,/fill,xstyle=1,ystyle=1,xr=[-1.,1.],yr=[-1.,1.],rgb_table=57,aspect_ratio=1, AXIS_STYLE=2,col=0,xtitle='x',xtickfont_size=18,ytickfont_size=18,C_LABEL_SHOW=0,position=[0.5,0.05,0.725,0.35],/current)

bar3=COLORBAR(TARGET=cg,rgb_table=57,ORIENTATION=0,position=[0.52,0.96,0.705,0.98],TEXTPOS = 1,major=3,taper=0,FONT_SIZE=18,TITLE='$FH_1$',range=[0,0.01],TICKVALUES=[0.0,0.005,0.01],TICKNAME=['0.0','0.005','0.01'])

ag=cg.axes
ag[0].tickvalue=[-0.5,-0,0.5]
ag[1].tickvalue=[-0.5,-0,0.5]
ag[0].tickname=['','','']
ag[1].tickname=['','','']

ah=ch.axes
ah[0].tickvalue=[-0.5,-0,0.5]
ah[1].tickvalue=[-0.5,-0,0.5]
ah[0].tickname=['','','']
ah[1].tickname=['','','']

ai=ci.axes
ai[0].tickvalue=[-0.5,-0,0.5]
ai[1].tickvalue=[-0.5,-0,0.5]
ai[1].tickname=['','','']

; ## Ion/rec frictional heating panels (FH2) ##

;c1=contour(reform(fh2_0)>levels2[0]+0.0001,x,y,c_value=levels2,/fill,xstyle=1,ystyle=1,xr=[-1.,1.],yr=[-1.,1.],rgb_table=51,aspect_ratio=1, AXIS_STYLE=2,col=0,xtitle='x',xtickfont_size=18,ytickfont_size=18,C_LABEL_SHOW=0,position=[0.59,0.65,0.77,0.95],/current)

;c2=contour(reform(fh2_1)>levels2[0]+0.0001,x,y,c_value=levels2,/fill,xstyle=1,ystyle=1,xr=[-1.,1.],yr=[-1.,1.],rgb_table=51,aspect_ratio=1, AXIS_STYLE=2,col=0,xtitle='x',xtickfont_size=18,ytickfont_size=18,C_LABEL_SHOW=0,position=[0.59,0.35,0.77,0.65],/current)

;c3=contour(reform(fh2_2)>levels2[0]+0.0001,x,y,c_value=levels2,/fill,xstyle=1,ystyle=1,xr=[-1.,1.],yr=[-1.,1.],rgb_table=51,aspect_ratio=1, AXIS_STYLE=2,col=0,xtitle='x',xtickfont_size=18,ytickfont_size=18,C_LABEL_SHOW=0,position=[0.59,0.05,0.77,0.35],/current)

;bara=COLORBAR(TARGET=c3,rgb_table=51,ORIENTATION=0,position=[0.61,0.96,0.75,0.98],TEXTPOS = 1,major=3,taper=0,FONT_SIZE=18,TITLE='$FH_2$',range=[0,0.001],TICKVALUES=[0.0,0.0005,0.001],TICKNAME=['0.0','0.0005','0.001'])

;a1=c1.axes
;a1[0].tickvalue=[-0.5,-0,0.5]
;a1[1].tickvalue=[-0.5,-0,0.5]
;a1[0].tickname=['','','']
;a1[1].tickname=['','','']

;a2=c2.axes
;a2[0].tickvalue=[-0.5,-0,0.5]
;a2[1].tickvalue=[-0.5,-0,0.5]
;a2[0].tickname=['','','']
;a2[1].tickname=['','','']

;a3=c3.axes
;a3[0].tickvalue=[-0.5,-0,0.5]
;a3[1].tickvalue=[-0.5,-0,0.5]
;a3[1].tickname=['','','']

; ## Ohmic heating panels ##

cl=contour(reform(ohm_1_0)>levels3[0]+0.0001,x,y,c_value=levels3,/fill,xstyle=1,ystyle=1,xr=[-1.,1.],yr=[-1.,1.],rgb_table=62,aspect_ratio=1, AXIS_STYLE=2,col=0,xtitle='x',xtickfont_size=18,ytickfont_size=18,C_LABEL_SHOW=0,position=[0.725,0.65,0.95,0.95],/current)

cm=contour(reform(ohm_1_1)>levels3[0]+0.0001,x,y,c_value=levels3,/fill,xstyle=1,ystyle=1,xr=[-1.,1.],yr=[-1.,1.],rgb_table=62,aspect_ratio=1, AXIS_STYLE=2,col=0,xtitle='x',xtickfont_size=18,ytickfont_size=18,C_LABEL_SHOW=0,position=[0.725,0.35,0.95,0.65],/current)

cn=contour(reform(ohm_1_2)>levels3[0]+0.0001,x,y,c_value=levels3,/fill,xstyle=1,ystyle=1,xr=[-1.,1.],yr=[-1.,1.],rgb_table=62,aspect_ratio=1, AXIS_STYLE=2,col=0,xtitle='x',xtickfont_size=18,ytickfont_size=18,C_LABEL_SHOW=0,position=[0.725,0.05,0.95,0.35],/current)

bar4=COLORBAR(TARGET=cl,rgb_table=62,ORIENTATION=0,position=[0.745,0.96,0.93,0.98],TEXTPOS = 1,major=3,taper=0,FONT_SIZE=18,TITLE='$Ohmic heating$',range=[0,0.2],TICKVALUES=[0.0,0.1,0.2],TICKNAME=['0.0','0.1','0.2'])

al=cl.axes
al[0].tickvalue=[-0.5,-0,0.5]
al[1].tickvalue=[-0.5,-0,0.5]
al[0].tickname=['','','']
al[1].tickname=['','','']

am=cm.axes
am[0].tickvalue=[-0.5,-0,0.5]
am[1].tickvalue=[-0.5,-0,0.5]
am[0].tickname=['','','']
am[1].tickname=['','','']

an=cn.axes
an[0].tickvalue=[-0.5,-0,0.5]
an[1].tickvalue=[-0.5,-0,0.5]
an[1].tickname=['','','']

ca.save,'fig_6.jpg'

end
