pro rdmpi,pv,datapath=datapath,current=current,flag_double=flag_double,$
           time_step=time_step,flag_az=flag_az,flag_te=flag_te,var=var,$
            h5read=h5read,nidealt=nidealt
  Compile_Opt DEFINT32  
  if(n_elements(datapath)eq 0 ) then datapath="tmp/" else datapath=datapath+"/"
  if(n_elements(current) eq 0) then current=0
  if(n_elements(flag_double) eq 0) then flag_double=1
  if(n_elements(time_step) eq 0) then time_step=-1
  if(n_elements(flag_az) eq 0) then flag_az = 0 ;; vector potential (for 2D)
  if(n_elements(flag_te) eq 0) then flag_te = 0 ;; temperature
  if(n_elements(var) eq 0) then var = [] ;; selected variables
  if(n_elements(h5read) eq 0) then h5read=0
  if(n_elements(nidealt) eq 0) then nidealt=1

  get_param2,datapath,info
  gm=info.gm
  pv=create_struct("info",info)
  ts=info.nt
  if time_step[0] eq -1 then n_read=ts else n_read=n_elements(time_step)

;;;;; test loop to prevent memory overflow
;if n_read GT 151 then n_read=151  

;define margin------------------------------
  margin=info.margin
;-------------------------------------------
  eqs=info.eqs
  flag_mhd=0
  flag_pip=0
  flag_afr=0
  case eqs of 
     'HD':nvar=5
     'MHD':begin
        nvar=8
        flag_mhd=1
     end
     'PIP':begin
        nvar=13
        flag_mhd=1
        flag_pip=1
     end
     "AFR":begin
;        nvar=29
        nvar=13
        flag_mhd=1
        flag_pip=1
        flag_afr=1
     end
  end

;Snow
if (n_elements(var) ne 0) then nvar=n_elements(var)
  
  pv=create_struct(pv,["eqs","fl_mhd","fl_pip","fl_afr"], $
                   eqs,flag_mhd,flag_pip,flag_afr)

  if h5read eq 0 then begin  
	tfile=file_search(datapath+"t.dac.*")
	dacget0s,tfile,t,narg=time_step
  endif
  if(eqs eq "AFR") then begin
     gridfile=file_search(datapath+"region.dac.*")
     dacget2s,gridfile,grid,narg=time_step
     pv=create_struct(pv,["grid"],grid)
  endif
  ix=info.ix
  jx=info.jx
  kx=info.kx

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Load the grid

;Using dac files
if (h5read eq 0) then begin
  if (where(tag_names(info) eq "MPI_DOMAINS"))[0] ne -1 then begin
     mpi_x=info.mpi_domains[0]
     mpi_y=info.mpi_domains[1]
     mpi_z=info.mpi_domains[2]
  endif else begin
     mpi_x=(mpi_y=(mpi_z=1))  
  endelse
  ix_m=(jx_m=(kx_m=1)) 
;  xfile=file_search(datapath+"x.dac.*")
  xfile=file_search(datapath+"x.dac."+string(indgen(mpi_x),form="(i4.4)"))
  ix_m=info.ix
  ix=margin[0]*2+mpi_x*(ix_m-margin[0]*2)
  x=findgen(ix)
  for n=0,mpi_x-1 do begin
     dacget1d,xfile[n],xc
     x0=n*(ix_m-2*margin[0])
     x[x0:x0+ix_m-1]=xc     
  endfor

  pv=create_struct(pv,"t",t,"x",x)
  ndim=pv.info.ndim
  if ndim ge 2 then begin     
     yfile=file_search(datapath+"y.dac."+string(indgen(mpi_y),form="(i4.4)"))
;       yfile=file_search(datapath+"y.dac.*")
     jx_m=info.jx
     jx=margin[1]*2+mpi_y*(jx_m-margin[1]*2)
     y=findgen(jx)
     for n=0,mpi_y-1 do begin
        dacget1d,yfile[n],yc
        y0=n*(jx_m-2*margin[1])
        y[y0:y0+jx_m-1]=yc     
     endfor
     pv=create_struct(pv,"y",y)
     if ndim ge 3 then begin
        zfile=file_search(datapath+"z.dac."+string(indgen(mpi_z),form="(i4.4)"))
        kx_m=info.kx
        kx=margin[2]*2+mpi_z*(kx_m-margin[2]*2)
        z=findgen(kx)
        for n=0,mpi_z-1 do begin
           dacget1d,zfile[n],zc
           z0=n*(kx_m-2*margin[2])
           z[z0:z0+kx_m-1]=zc     
        endfor
        pv=create_struct(pv,"z",z)
     endif
  endif
  pv=create_struct(pv,"ix",n_elements(x))
  if ndim ge 2 then jx=n_elements(y) else jx=1  
  if ndim ge 3 then kx=n_elements(z) else kx=1  
  pv=create_struct(pv,"jx",jx,"kx",kx)  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Get time 
  n_cpu=mpi_x*mpi_y*mpi_z
  var_names=[]
  if time_step[0] eq -1 then begin
     time_step=indgen(n_read)
  endif
  output=dblarr(ix,jx,kx)


  if(0 eq 1 ) then begin
     add_pv,pv,["ro_amb","en_amb","mx_amb","my_amb","mz_amb",$
                "bx_amb","by_amb","bz_amb"],$
            1,var_names,mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m,$
            time_step,datapath    
     add_pv,pv,["ro_mhd","en_mhd","mx_mhd","my_mhd","mz_mhd",$
                "bx_mhd","by_mhd","bz_mhd"],$
            1,var_names,mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m,$
            time_step,datapath    
  endif  

;Snow
if (n_elements(var) eq 0) then begin
  if (flag_pip eq 1 or flag_mhd eq 0) then begin
     add_pv,pv,["ro_n","en_n","mx_n","my_n","mz_n"] ,0,var_names,$
            mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m,$
            time_step,datapath    
  endif  
  if (flag_mhd eq 1) then begin
     add_pv,pv,["ro_p","en_p","mx_p","my_p","mz_p","bx","by","bz"],$
            1,var_names,mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m,$
            time_step,datapath
  endif

endif else begin
     add_pv_var,pv,var ,0,var_names,$
            mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m,$
            time_step,datapath   
;help,pv,/str
endelse

;===============================================-
; Non-ideal terms
;===============================================
if (nidealt eq 1) then begin

  if (((info.flag_ir mod 10) ge 1) and (info.flag_pip ge 1))  then begin
     flag_ir=1
     if flag_ir eq 1 then begin
        ir_flag=info.flag_ir mod 10
        if ir_flag eq 1 then begin
           ion=dblarr(ix,jx,kx,n_read)
           rec=dblarr(ix,jx,kx,n_read)
        endif else begin
           ion=dblarr(ix,jx,kx,n_read)
           rec=dblarr(ix,jx,kx,n_read)
        endelse                   
print, 'ION'
	pv=create_struct(pv,["flag_ir"], info.flag_ir)
	pv=create_struct(pv,["flag_ir_type"], info.flag_ir_type) 
	pv=create_struct(pv,["T0"], info.T_norm) 
        for np=0,n_read-1 do begin
        files=file_search(datapath+"/"+string(time_step[np],form="(i4.4)")+"ion.dac.*")
	mpi_read,ion1,files,mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m
          ion[*,*,*,np]=ion1
        files=file_search(datapath+"/"+string(time_step[np],form="(i4.4)")+"rec.dac.*")
	mpi_read,rec1,files,mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m
          rec[*,*,*,np]=rec1
        endfor
        pv=create_struct(pv,["ion"], reform(ion))  
        pv=create_struct(pv,["rec"], reform(rec))        
     endif
  endif

  if current eq 1 then begin
	print,'Getting eta'
	eta=dblarr(ix,jx,kx)
        files=file_search(datapath+"/"+string(time_step[0],form="(i4.4)")+"et.dac.*")
        mpi_read,eta1,files,mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m
          eta[*,*,*]=eta1
        ;endfor
        pv=create_struct(pv,["eta"], reform(eta))

        print,'CURRENT ONLY WORKS IN 2D AT PRESENT'
        j_x=dblarr(ix,jx,kx,n_read)
        j_y=dblarr(ix,jx,kx,n_read)
        j_z=dblarr(ix,jx,kx,n_read)
        dx=abs(pv.x[1]-pv.x[0]) ; grid point size
        dy=abs(pv.y[1]-pv.y[0])
	dz=abs(pv.z[1]-pv.z[0])
        for np=0,n_read-1 do begin
                j_x[2:(ix-3),2:(jx-3),2:(kx-3),np]=(-pv.Bz[2:(ix-3),4:(jx-1),2:(kx-3),np]+8.*pv.Bz[2:(ix-3),3:(jx-2),2:(kx-3),np]-8.*pv.Bz[2:(ix-3),1:(jx-4),2:(kx-3),np]+pv.Bz[2:(ix-3),0:(jx-5),2:(kx-3),np])/(12.*dy) +$
			(pv.By[2:(ix-3),2:(jx-3),4:(kx-1),np]+8.*pv.By[2:(ix-3),2:(jx-3),3:(kx-2),np]-8.*pv.By[2:(ix-3),2:(jx-3),1:(kx-4),np]+pv.By[2:(ix-3),2:(jx-3),0:(kx-5),np])/(12.*dz) 

                j_y[2:(ix-3),2:(jx-3),2:(kx-3),np]=(-pv.Bz[4:(ix-1),2:(jx-3),2:(kx-3),np]+8.*pv.Bz[3:(ix-2),2:(jx-3),2:(kx-3),np]-8.*pv.Bz[1:(ix-4),2:(jx-3),2:(kx-3),np]+pv.Bz[0:(ix-5),2:(jx-3),2:(kx-3),np])/(12.*dx) -$
			(-pv.Bx[2:(ix-3),2:(jx-3),4:(kx-1),np]+8.*pv.Bx[2:(ix-3),2:(jx-3),3:(kx-2),np]-8.*pv.Bx[2:(ix-3),2:(jx-3),1:(kx-4),np]+pv.Bx[2:(ix-3),2:(jx-3),0:(kx-5),np])/(12.*dz)

                j_z[2:(ix-3),2:(jx-3),2:(kx-3),np]=(-pv.By[4:(ix-1),2:(jx-3),2:(kx-3),np]+8.*pv.By[3:(ix-2),2:(jx-3),2:(kx-3),np]-8.*pv.By[1:(ix-4),2:(jx-3),2:(kx-3),np]+pv.By[0:(ix-5),2:(jx-3),2:(kx-3),np])/(12.*dx) - $
                (-pv.Bx[2:(ix-3),4:(jx-1),2:(kx-3),np]+8.*pv.Bx[2:(ix-3),3:(jx-2),2:(kx-3),np]-8.*pv.Bx[2:(ix-3),1:(jx-4),2:(kx-3),np]+pv.Bx[2:(ix-3),0:(jx-5),2:(kx-3),np])/(12.*dy)
        endfor
        pv = create_struct(pv,["j_x"],reform(j_x))
        pv = create_struct(pv,["j_y"],reform(j_y))
        pv = create_struct(pv,["j_z"],reform(j_z))
  endif


endif

  close,/all

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
endif else begin
    RESOLVE_ROUTINE,'h5get'
;Read grid grom h5 files
    ndim=pv.info.ndim
    fpath=datapath+'/t0000.h5' 
    fpath=datapath+"/t"+string(string(time_step),form="(i4.4)")+'.h5'
print,'ONLY READING ONE TIME STEP FOR TESTS'   
    readvars=['xgrid']
    if ndim ge 2 then readvars=[readvars,'ygrid']
    if ndim ge 3 then readvars=[readvars,'zgrid']
    h5get,pv,fpath,readvars,1

    ;Read conserved variables
    if (n_elements(var) eq 0) then begin
      if (flag_pip eq 1 or flag_mhd eq 0) then begin
;         h5get,pv,fpath,["ro_n","en_n","mx_n","my_n","mz_n"],0
         h5get,pv,fpath,["ro_n","pr_n","vx_n","vy_n","vz_n"],0
      endif  
      if (flag_mhd eq 1) then begin
;         h5get,pv,fpath,["ro_p","en_p","mx_p","my_p","mz_p","bx","by","bz"],0
         h5get,pv,fpath,["ro_p","pr_p","vx_p","vy_p","vz_p","bx","by","bz"],0
      endif

    endif else begin
         h5get,pv,fpath,var,0
    ;help,pv,/str
    endelse

;   Non-ideal terms
    if n_elements(var) eq 0 then begin
        if (info.flag_rad ge 1) then begin
        edref=dblarr(ix,jx,kx,n_read)
        print, 'Rad_cooling Loaded'
        h5get,pv,fpath,["edref_m"],1
        radrho=info.radrhoref
        radt=info.rad_ts
        pv=create_struct(pv,["radt"], reform(radt))
        pv=create_struct(pv,["radrho"], reform(radrho)) 
        endif   
	if (((info.flag_ir mod 10) ge 1) and (info.flag_pip ge 1)) then begin
	 
		if (info.flag_ir_type eq 0) then begin
		print, 'Losses loading'
		h5get,pv,fpath,["aheat"],1
		h5get,pv,fpath,["ion_loss"],1
		endif  
		if (info.flag_ir eq 4) then begin
		print, 'Hydrogen levels loading'
		h5get,pv,fpath,["nexcite1","nexcite2","nexcite3","nexcite4","nexcite5","nexcite6"],1
		endif        
	endif
    endif

endelse
end
