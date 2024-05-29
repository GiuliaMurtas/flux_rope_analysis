;Departure coefficient test
rdmpi,ds,datapath='../kink_instability_PIP_0',time_step=[13],var=['ro_p','ro_n','pr_p','pr_n']

Tp=ds.pr_p(*,*,441)/ds.ro_p(*,*,441)*5.0/6.0
T0=12855.0/Tp(0,0)

Te=Tp*T0 ;Proton/Electron temperature in Kelvin
Tn=T0*ds.pr_p(*,*,441)/ds.ro_p(*,*,441)*5.0/3.0 ;Neutral Temperature in Kelvin

T_to_eV=1.0/11606.0

F=2.6e-19/sqrt(T_to_eV*Te)
G=2.91e-14*(1.0/(0.232+13.6/Te/T_to_eV))*(13.6/Te/T_to_eV)^0.39*exp(-13.6/Te/T_to_eV)
xi_i_steady=1.0/((F/G)+1.0)

xi_i=ds.ro_p(*,*,441)/(ds.ro_p(*,*,441)+ds.ro_n(*,*,441))

file = 'departure_coef_PIP0_t13.h5'
fid = H5F_CREATE(file)
print,file

varlist=['xi_i','xi_i_steady','xgrid',$
    'ygrid']
for var=0,n_elements(varlist)-1 do begin
    varname=varlist(var)

    if varname eq 'xi_i' then data = xi_i
    if varname eq 'xi_i_steady' then data = xi_i_steady
    if varname eq 'xgrid' then data = ds.x
    if varname eq 'ygrid' then data = ds.y

    print,varname,max(data),min(data)

    datatype_id = H5T_IDL_CREATE(data)
    dataspace_id = H5S_CREATE_SIMPLE(size(data,/DIMENSIONS))
    dataset_id = H5D_CREATE(fid,$
    varname,datatype_id,dataspace_id)
    H5D_WRITE,dataset_id,data
    H5D_CLOSE,dataset_id
    H5S_CLOSE,dataspace_id
    H5T_CLOSE,datatype_id
endfor

H5F_CLOSE,fid

END

