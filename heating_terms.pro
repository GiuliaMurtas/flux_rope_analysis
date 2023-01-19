;God Damn it, Giulia!
popul=1

filename0 = '../kink_instability_MHD_0/' ;MHD
filename1 = '../kink_instability_PIP_0/' ; P1
filename2 = '../kink_instability_PIP_1/' ; P2

ohm0 = []
ohm1 = []
ohm2 = []

fric1 = []
fric2 = []

ir1 = []
ir2 = []

meantp0=[]
meantp1 = []
meantn1 = []
meantp2 = []
meantn2 = []

time0 = []
time1 = []
time2 = []

;# Ohmic heating##

;# MHD case - M1 ##
if popul eq 1 then begin
for t =0,20 do begin ;20
    print,t
    rdmpi,ds,datapath=filename0,time_step=t,var=['bx','by','bz'],/current
    
    dxm = ds.x[1] - ds.x[0]
    dym = ds.y[1] - ds.y[0]
    dzm = ds.z[1] - ds.z[0]
    jtot2 = ds.j_z^2 + ds.j_y^2 +ds.j_x^2

    ohm_element=total(ds.eta[7:502,7:502,7:874]*jtot2[7:502,7:502,7:874])*dxm*dym*dzm
    ohm0=[ohm0,ohm_element]
    time0=[time0,ds.t[0]]
    
    rdmpi,ds,datapath=filename0,time_step=t,var=['pr_p','ro_p']
    
    Tp = (5.0/6.0)*(ds.pr_p[7:502,7:502,7:874]/ds.ro_p[7:502,7:502,7:874])
    
    meanTp=mean(Tp)
    meantp0=[meantp0,meanTp]
endfor
    
;## PIP case - P1 ##
for t =0,18 do begin ;18
    print,t
    rdmpi,ds,datapath=filename1,time_step=t,var=['bx','by','bz'],/current
    
    dxm = ds.x[1] - ds.x[0]
    dym = ds.y[1] - ds.y[0]
    dzm = ds.z[1] - ds.z[0]
    jtot2 = ds.j_z^2 + ds.j_y^2 +ds.j_x^2

    ohm_element=total(ds.eta[7:502,7:502,7:874]*jtot2[7:502,7:502,7:874])*dxm*dym*dzm
    ohm1=[ohm1,ohm_element]
    time1=[time1,ds.t[0]]

endfor

;## PIP case - P2 ##
for t =0,19 do begin ;19
    print,t
    rdmpi,ds,datapath=filename2,time_step=t,var=['bx','by','bz'],/current
    
    dxm = ds.x[1] - ds.x[0]
    dym = ds.y[1] - ds.y[0]
    dzm = ds.z[1] - ds.z[0]
    jtot2 = ds.j_z^2 + ds.j_y^2 +ds.j_x^2

    ohm_element=total(ds.eta[7:502,7:502,7:874]*jtot2[7:502,7:502,7:874])*dxm*dym*dzm
    ohm2=[ohm2,ohm_element]
    time2=[time2,ds.t[0]]

endfor   
save,ohm0,ohm1,ohm2,meanTp0,time0,time1,time2,filename='fig_8_oh.sav'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;    
;## Collisional frictional heating ##

;## PIP case - P1 ##
for t =0,18 do begin ;18
    print,t
    rdmpi,ds,datapath=filename1,time_step=t,var=['pr_n','pr_p','ro_n','ro_p','vx_p','vy_p','vz_p','vx_n','vy_n','vz_n']

    dxm = ds.x[1] - ds.x[0]
    dym = ds.y[1] - ds.y[0]
    dzm = ds.z[1] - ds.z[0]
        
    Tn = (5.0/3.0)*(ds.pr_n[7:502,7:502,7:874]/ds.ro_n[7:502,7:502,7:874])
    Tp = (5.0/6.0)*(ds.pr_p[7:502,7:502,7:874]/ds.ro_p[7:502,7:502,7:874])
    
    meanTp=mean(Tp)
    meanTn=mean(Tn)
    meantp1=[meantp1,meanTp]
    meantn1=[meantn1,meanTn]
    
    vD = (ds.vx_n[7:502,7:502,7:874]-ds.vx_p[7:502,7:502,7:874])^2 $
    + (ds.vy_n[7:502,7:502,7:874]-ds.vy_p[7:502,7:502,7:874])^2 $
    + (ds.vz_n[7:502,7:502,7:874]-ds.vz_p[7:502,7:502,7:874])^2
    
    alpha = 1.0*sqrt(0.5*(Tn+Tp))*sqrt(1 + (9*!pi/64.0)*(5.0/6.0)*(vD/(Tn+Tp)))

    fric_element = total(0.5*alpha*ds.ro_n[7:502,7:502,7:874]*ds.ro_p[7:502,7:502,7:874]*vD)*dxm*dym*dzm
    fric1=[fric1,fric_element]
    
    vp2 = ds.vx_p[7:502,7:502,7:874]^2 + ds.vy_p[7:502,7:502,7:874]^2 $
    + ds.vz_p[7:502,7:502,7:874]^2
    vn2 = ds.vx_n[7:502,7:502,7:874]^2 + ds.vy_n[7:502,7:502,7:874]^2 $
    + ds.vz_n[7:502,7:502,7:874]^2
    vnvp = ds.vx_n[7:502,7:502,7:874]*ds.vx_p[7:502,7:502,7:874] $
    + ds.vy_n[7:502,7:502,7:874]*ds.vy_p[7:502,7:502,7:874] $
    + ds.vz_n[7:502,7:502,7:874]*ds.vz_p[7:502,7:502,7:874]
    
    ir_element = ds.rec[7:502,7:502,7:874]*ds.ro_p[7:502,7:502,7:874]*vp2 $
    - (ds.rec[7:502,7:502,7:874]*ds.ro_p[7:502,7:502,7:874] $
    + ds.ion[7:502,7:502,7:874]*ds.ro_n[7:502,7:502,7:874])*vnvp $
    + ds.ion[7:502,7:502,7:874]*ds.ro_n[7:502,7:502,7:874]*vn2

    ir = total(ir_element)*dxm*dym*dzm
    ir1=[ir1,ir]
endfor
    

;## PIP case - P2 ##
for t =0,19 do begin ;19
    print,t
    rdmpi,ds,datapath=filename2,time_step=t,var=['pr_n','pr_p','ro_n','ro_p','vx_p','vy_p','vz_p','vx_n','vy_n','vz_n']

    dxm = ds.x[1] - ds.x[0]
    dym = ds.y[1] - ds.y[0]
    dzm = ds.z[1] - ds.z[0]
        
    Tn = (5.0/3.0)*(ds.pr_n[7:502,7:502,7:874]/ds.ro_n[7:502,7:502,7:874])
    Tp = (5.0/6.0)*(ds.pr_p[7:502,7:502,7:874]/ds.ro_p[7:502,7:502,7:874])
    
    meanTp=mean(Tp)
    meanTn=mean(Tn)
    meantp2=[meantp2,meanTp]
    meantn2=[meantn2,meanTn]
        
    vD = (ds.vx_n[7:502,7:502,7:874]-ds.vx_p[7:502,7:502,7:874])^2 $
    + (ds.vy_n[7:502,7:502,7:874]-ds.vy_p[7:502,7:502,7:874])^2 $
    + (ds.vz_n[7:502,7:502,7:874]-ds.vz_p[7:502,7:502,7:874])^2
    
    alpha = 1.0*sqrt(0.5*(Tn+Tp))*sqrt(1 + (9*!pi/64.0)*(5.0/6.0)*(vD/(Tn+Tp)))

    fric_element = total(0.5*alpha*ds.ro_n[7:502,7:502,7:874]*ds.ro_p[7:502,7:502,7:874]*vD)*dxm*dym*dzm
    fric2=[fric2,fric_element]
    
    vp2 = ds.vx_p[7:502,7:502,7:874]^2 + ds.vy_p[7:502,7:502,7:874]^2 $
    + ds.vz_p[7:502,7:502,7:874]^2
    vn2 = ds.vx_n[7:502,7:502,7:874]^2 + ds.vy_n[7:502,7:502,7:874]^2 $
    + ds.vz_n[7:502,7:502,7:874]^2
    vnvp = ds.vx_n[7:502,7:502,7:874]*ds.vx_p[7:502,7:502,7:874] $
    + ds.vy_n[7:502,7:502,7:874]*ds.vy_p[7:502,7:502,7:874] $
    + ds.vz_n[7:502,7:502,7:874]*ds.vz_p[7:502,7:502,7:874]
    
    ir_element = ds.rec[7:502,7:502,7:874]*ds.ro_p[7:502,7:502,7:874]*vp2 $
    - (ds.rec[7:502,7:502,7:874]*ds.ro_p[7:502,7:502,7:874] $
    + ds.ion[7:502,7:502,7:874]*ds.ro_n[7:502,7:502,7:874])*vnvp $
    + ds.ion[7:502,7:502,7:874]*ds.ro_n[7:502,7:502,7:874]*vn2

    ir = total(ir_element)*dxm*dym*dzm
    ir2=[ir2,ir]
endfor

save,fric1,fric2,ir1,ir2,meantp1,meantn1,meantp2,meantn2,filename='fig_8_fh.sav'

endif
if popul eq 0 then begin
	restore,'fig_8_oh.sav'
	restore,'fig_8_fh.sav'
endif

file = 'fig_8_data.h5'
outvar=['ohm0','ohm1','ohm2','time0','time1','time2','fric1','fric2','ir1','ir2','meanTp0','meanTp1','meanTn1','meanTp2','meanTn2']
fid = H5F_CREATE(file)

for i=0,n_elements(outvar)-1 do begin
    if outvar(i) eq 'ohm0' then data = ohm0
    if outvar(i) eq 'ohm1' then data = ohm1
    if outvar(i) eq 'ohm2' then data = ohm2
    if outvar(i) eq 'time0' then data = time0
    if outvar(i) eq 'time1' then data = time1
    if outvar(i) eq 'time2' then data = time2
    if outvar(i) eq 'fric1' then data = fric1
    if outvar(i) eq 'ir1' then data = ir1
    if outvar(i) eq 'fric2' then data = fric2
    if outvar(i) eq 'ir2' then data = ir2
    if outvar(i) eq 'meanTp0' then data = meanTp0
    if outvar(i) eq 'meanTp1' then data = meanTp1
    if outvar(i) eq 'meanTn1' then data = meanTn1
    if outvar(i) eq 'meanTp2' then data = meanTp2
    if outvar(i) eq 'meanTn2' then data = meanTn2
    
    ;; get data type and space, needed to create the dataset
    datatype_id = H5T_IDL_CREATE(data)
    dataspace_id = H5S_CREATE_SIMPLE(size(data,/DIMENSIONS))
    dataset_id = H5D_CREATE(fid,outvar(i),datatype_id,dataspace_id)

    ;; write data to dataset
    H5D_WRITE,dataset_id,data

    ;; close all open identifiers
    H5D_CLOSE,dataset_id
    H5S_CLOSE,dataspace_id
    H5T_CLOSE,datatype_id
endfor

H5F_CLOSE,fid

END
