#if NDIM==3
subroutine star_formation(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  use cooling_module, ONLY: XH=>X
  use constants, only: Myr2sec, Gyr2sec, mH, pi, rhoc, twopi
  use random
  use mpi_mod
#ifdef RT
#ifdef NTRACEGROUPS
  use rt_parameters, only: minind,rt_tracemassbins
  use SED_module,     only: get_tracer_group_from_halo
#endif
#endif 
  implicit none
#ifndef WITHOUTMPI
  integer::info,info2,dummy_io
  integer,parameter::tag=1120
#endif
  integer::ilevel
  !----------------------------------------------------------------------
  ! Description: This subroutine spawns star-particle of constant mass
  ! using a Poisson probability law if some gas condition are fulfilled.
  ! It modifies hydrodynamic variables according to mass conservation
  ! and assumes an isothermal transformation...
  ! On exit, the gas velocity and sound speed are unchanged.
  ! New star particles are synchronized with other collisionless particles.
  ! Array flag2 is used as temporary work space.
  ! Yann Rasera  10/2002-01/2003
  !----------------------------------------------------------------------
  ! local constants
  real(dp)::d0,mgas,mcell
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:twotondim,1:3)::xc
  ! other variables
  integer ::ncache,nnew,ivar,ngrid,icpu,index_star,ndebris_tot,ilun=10
  integer ::igrid,ix,iy,iz,ind,i,n,iskip,nx_loc,idim
  integer ::ntot,ntot_all,nstar_corrected,ncell
  logical ::ok_free
  real(dp)::d,x,y,z,u,v,w,e,tg,zg,zg_ox
    !shubham
  real(dp)::estar,Tkstar
  !shbham 

  real(dp)::mstar,mstar1,dstar,tstar,nISM,nCOM=0.0d0,phi_t,phi_x,theta,sigs,scrit,b_turb,zeta
  real(dp)::T2,nH,T_poly,cs2,cs2_poly,trel,t_dyn,t_ff,tdec,uvar
  real(dp)::ul,ur,fl,fr,trgv,alpha0
  real(dp)::d_gmc,mstar_nsn,scale_msun
  real(dp)::sigma2,sigma2_comp,sigma2_sole,lapld,flong,ftot,pcomp=0.3d0
  real(dp)::divv,divv2,curlv,curlva,curlvb,curlvc,curlv2
  real(dp)::birth_epoch,factG
  real(kind=8)::mlost_all,mtot_all
#ifndef WITHOUTMPI
  real(kind=8)::mlost,mtot
#endif
  real(kind=8)::RandNum,PoissMean
  real(dp)::v_kick(1:3),v_kick_mag
  real(dp),dimension(1:3)::skip_loc
  real(dp)::dx,dx_loc,scale,vol_loc,dx_min,vol_min,d1,d2,d3,d4,d5,d6
  real(dp)::mdebris
  real(dp),dimension(1:nvector)::sfr_ff, alpha_fk=0.0
  integer ,dimension(1:ncpu,1:IRandNumSize)::allseed
  integer ,dimension(1:nvector),save::ind_grid,ind_cell,ind_cell2,nstar
  integer ,dimension(1:nvector),save::ind_grid_new,ind_cell_new,ind_part
  integer ,dimension(1:nvector),save::ind_debris
  integer ,dimension(1:nvector,0:twondim)::ind_nbor
  logical ,dimension(1:nvector),save::ok,ok_new=.true.
  integer ,dimension(1:ncpu)::ntot_star_cpu,ntot_star_all
  character(LEN=80)::filename,filedir,fileloc,filedirini
  character(LEN=5)::nchar,ncharcpu
  logical::file_exist
  integer  ::idx3,idx3_new
  real(dp)::Zscale
  real(dp)::sigma_d,chi,tau_c,sss,f_h2
#ifdef SOLVERmhd
  real(dp)::bx1,bx2,by1,by2,bz1,bz2,A,B,C,emag,beta,fbeta
#endif
#if NENER>0
  integer::irad
#endif
#ifdef NCHEM
  real(dp),dimension(1:nchem) :: chem1=0.0
  integer  :: iche
#endif
#ifdef NTRACEGROUPS
  integer::dumpid
#endif
#ifdef POP3
  ! for PopIII stars
  real(dp) :: r3pc
  logical  :: pop3_ok
  integer  :: kk
  real(dp),dimension(1:1000),save              :: mass_pop3
  integer ,dimension(1:1000),save              :: idx3icell
#endif

  !tracer 
#ifdef MC_tracer
  ! MC Tracer patch
  integer :: nattach, ip, ipart
  real(dp) :: delta_m_over_m
  logical, dimension(1:nvector), save :: tok
  integer, dimension(1:nvector), save :: itracer, istar_tracer
  real(dp), dimension(1:nvector, 1:3), save :: xstar
  real(dp), dimension(1:nvector), save :: proba
  logical :: move_tracer
  ! End MC Tracer patch
#endif 
  !tracer 

  if(numbtot(1,ilevel)==0) return
  if(.not. hydro)return
  if(ndim.ne.3)return
  if(static)return
  if (sf_lmax .and. ilevel.lt.nlevelmax_current) return ! Only form stars at levelmax

  if(verbose)write(*,*)' Entering star_formation'

  if(sf_log_properties.and.ifout.gt.1) then
     call title(ifout-1,nchar)
     if(IOGROUPSIZEREP>0) then
        call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
        filedirini='output_'//TRIM(nchar)//'/'
        filedir='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/'
     else
        filedir='output_'//TRIM(nchar)//'/'
     endif
     filename=TRIM(filedir)//'stars_'//TRIM(nchar)//'.out'
     ilun=myid+103
     call title(myid,nchar)
     fileloc=TRIM(filename)//TRIM(nchar)
     ! Wait for the token
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if (mod(myid-1,IOGROUPSIZE)/=0) then
           call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
                & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
        end if
     endif
#endif

     inquire(file=fileloc,exist=file_exist)
     if((.not.file_exist).or.(abs(t-trestart).lt.dtnew(ilevel))) then
        open(ilun, file=fileloc, form='formatted')
        write(ilun,'(A24)',advance='no') '# event id  ilevel  mp  '
        do idim=1,ndim
           write(ilun,'(A2,I1,A2)',advance='no') 'xp',idim,'  '
        enddo
        do idim=1,ndim
           write(ilun,'(A2,I1,A2)',advance='no') 'vp',idim,'  '
        enddo
        do ivar=1,nvar
           if(ivar.ge.10) then
              write(ilun,'(A1,I2,A2)',advance='no') 'u',ivar,'  '
           else
              write(ilun,'(A1,I1,A2)',advance='no') 'u',ivar,'  '
           endif
        enddo
        write(ilun,'(A5)',advance='no') 'tag  '
        write(ilun,'(A1)') ' '
     else
        open(ilun, file=fileloc, status="old", position="append", action="write", form='formatted')
     endif
  endif

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_msun = scale_l**3*scale_d/1.989d33

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim
  dx_min=(0.5D0**nlevelmax)*scale
  vol_min=dx_min**ndim

  trel=sf_trelax*Myr2sec/scale_t ! relaxation timescale

  ! ISM density threshold from H/cc to code units
  nISM = n_star
  if(cosmo)then
     nCOM = del_star*omega_b*rhoc*(h0/100)**2/aexp**3*XH/mH
     nISM = MAX(nCOM,nISM)
  endif
  d0   = nISM/scale_nH

  !------------------------------------------------------
  ! Set the star particle mass from the number of SN [TK]
  !------------------------------------------------------
  if(nsn2mass>0.and.fstar_min<0)then
     ! Mesh spacing
     mstar_nsn  = (nsn2mass*M_SNII)/eta_sn/scale_msun
     ! ISM density threshold from H/cc to code units    
     mstar      = n_star/(scale_nH*aexp**3)*vol_min
     fstar_min  = mstar_nsn/mstar
     if(myid==1) write(*,*) ">>>TKNOTE: Mstar,min=",mstar_nsn*scale_msun,fstar_min,nsn2mass,M_SNII,eta_sn
  endif

  !tracer 
#ifdef MC_tracer
  ! MC Tracer patch
  if(MC_tracer) then 
    tok = .false.
    nattach = 0
  endif 
  ! End MC Tracer patch
#endif 
  !tracer 

  ! Initial star particle mass
  if(m_star < 0d0)then
     mstar=n_star/(scale_nH*aexp**3)*vol_min*fstar_min
  else
     mstar=m_star*mass_sph
  endif
  dstar=mstar/vol_loc

  factG = 1d0
  if(cosmo) factG = 3d0/4d0/twopi*omega_m*aexp

  ! Birth epoch as proper time
  if(use_proper_time)then
     birth_epoch=texp
  else
     birth_epoch=t
  endif
  d_gmc = MAX(nCOM,n_gmc)/scale_nH

  ! Cells center position relative to grid center position
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)*dx
     xc(ind,2)=(dble(iy)-0.5D0)*dx
     xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! If necessary, initialize random number generator
  if(localseed(1)==-1)then
     call rans(ncpu,iseed,allseed)
     localseed=allseed(myid,1:IRandNumSize)
  end if

  !------------------------------------------------
  ! Convert hydro variables to primitive variables
  !------------------------------------------------
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        do i=1,ngrid
           d=uold(ind_cell(i),1)
           u=uold(ind_cell(i),2)/d
           v=uold(ind_cell(i),3)/d
           w=uold(ind_cell(i),4)/d
           e=uold(ind_cell(i),5)
#ifdef SOLVERmhd
           bx1=uold(ind_cell(i),6)
           by1=uold(ind_cell(i),7)
           bz1=uold(ind_cell(i),8)
           bx2=uold(ind_cell(i),nvar+1)
           by2=uold(ind_cell(i),nvar+2)
           bz2=uold(ind_cell(i),nvar+3)
           e=e-0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)
#endif
           e=e-0.5d0*d*(u**2+v**2+w**2)
#if NENER>0
           do irad=0,nener-1
              e=e-uold(ind_cell(i),inener+irad)
           end do
#endif
           uold(ind_cell(i),1)=d
           uold(ind_cell(i),2)=u
           uold(ind_cell(i),3)=v
           uold(ind_cell(i),4)=w
           uold(ind_cell(i),5)=e/d
        end do
        do ivar=imetal,nvar
           do i=1,ngrid
              d=uold(ind_cell(i),1)
              w=uold(ind_cell(i),ivar)/d
              uold(ind_cell(i),ivar)=w
           end do
        end do
     end do
  end do

! get values of uold for density and velocities in virtual boundaries
#ifndef WITHOUTMPI
  do ivar=1,4
     call make_virtual_fine_dp(uold(1,ivar),ilevel)
  end do
#endif

  !------------------------------------------------
  ! Compute number of new stars in each cell
  !------------------------------------------------
  ntot=0
  ndebris_tot=0
  idx3=0  ! index for PopIII stars
  ! Loop over grids
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     ! Star formation criterion ---> logical array ok(i)
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        sfr_ff(:)=eps_star
        ! Flag leaf cells
        do i=1,ngrid
           ok(i)=son(ind_cell(i))==0
        end do
        ! Geometrical criterion
        if(ivar_refine>0)then
           do i=1,ngrid
              d=uold(ind_cell(i),ivar_refine)
              if(d<=var_cut_refine)ok(i)=.false.
           end do
        endif
        if(sf_virial)then
           call get_sfeff_virial()     ! ->ok, sfr_ff
        else
           call get_sfeff_non_virial() ! ->ok, sfr_ff
        endif    
        ! Calculate number of new stars in each cell using Poisson statistics
        do i=1,ngrid
           nstar(i)=0
           if(ok(i))then
              ! Compute mean number of events
              d=uold(ind_cell(i),1)
              !shubham 
              estar  = uold(ind_cell(i),5) 
              Tkstar = (gamma - 1.0d0) * estar * scale_T2  
              !shubham  
              
              if(KMT09)then  !Agertz: calculate f_H2, which enters mcell below.
                 Zscale=1
                 if(metal) then
                    Zscale=1.06d0*uold(ind_cell(i),imetal)/0.02d0
                    if(is_oxygen) &
                    Zscale=Zscale+2.09d0*uold(ind_cell(i),imetal+1)/0.02d0
                 endif
                 !Average, solar mix in gas (Asplund)
                 sigma_d=1.0d-21*Zscale !Dust coss section per H nucleus
                 chi=2.3*(1.+3.1*Zscale**0.365)/3.    !See Kuhlen for constants
                 tau_c=d*dx_loc*scale_d*scale_l*sigma_d/2.3d-24
                 sss=log(1.+0.6*chi+0.01*chi**2)/0.6/tau_c
                 f_h2=1.-3.*sss/4./(1.+0.25*sss)
                 f_h2=max(f_h2,fh2_min)  !floor
                 if(d.ge.fh2_rho/scale_nH) then  !Enforce threshold where we assume 100% molecular
                     f_h2=1.0
                  endif
              else
                 f_h2=1.0 !accounts for self-shielded gas
              endif
              mcell=d*vol_loc
              ! Free fall time of an homogeneous sphere
              tstar= .5427d0*sqrt(1.0d0/(factG*max(d,smallr)))
              ! Gas mass to be converted into stars
              mgas=dtnew(ilevel)*(sfr_ff(i)/tstar)*mcell*f_h2
#ifdef POP3
              ! for PopIII star formation
              if(metal) zg = uold(ind_cell(i),imetal)
              pop3_ok = pop3.and.(zg.lt.Zcrit_pop3)
              if(pop3_ok) then
                 if(pop3_mass>0)then
                    mstar1 = pop3_mass/scale_msun
                    PoissMean = mgas/mstar1
                 else
                    call pop3_mass_random(mstar1)
                    mstar1 = mstar1/scale_msun
                    PoissMean = mgas/(100./scale_msun) ! 100: characteristic scale, for random sampling
                 endif
              else
#endif
                 mstar1 = mstar
                 ! Poisson mean
                 PoissMean = mgas/mstar1
#ifdef POP3
              endif
#endif
              if((trel>0.).and.(.not.cosmo)) PoissMean = PoissMean*min((t/trel), 1.0d0)
              ! Compute Poisson realisation
              call poissdev(localseed,PoissMean,nstar(i))
              ! Compute depleted gas mass
              mgas=nstar(i)*mstar1
#ifdef POP3
              ! Allow 1 popIII star formation per event
              if(pop3_ok.and.(nstar(i).ge.1)) then
                 if(mcell*scale_msun.ge.1d3.and.mgas<0.9*mcell) then ! for random sampling
                    nstar(i)=1
                    idx3 = idx3 + 1
                    if(idx3.gt.10000)then
                       write(*,*)' ERROR in star formation: increase PopIII reservoir'
                       call clean_stop
                    endif
                    mass_pop3(idx3) = mstar1
                    idx3icell(idx3) = ind_cell(i)
                 else
                    nstar(i)=0
                 endif
              endif 
#endif
              ! Security to prevent more than 90% of gas depletion
              if (mgas > 0.9d0*mcell) then
                 nstar_corrected=int(0.9d0*mcell/mstar1)
                 mstar_lost=mstar_lost+(nstar(i)-nstar_corrected)*mstar1
                 nstar(i)=nstar_corrected
              endif
#ifdef POP3
              ! check if there is any active Pop III stars within 3 pc radius
              r3pc = 3.*3.08d18/(scale_l*boxlen)
              r3pc = max(r3pc,3*0.5**nlevelmax)
              if(nstar(i)>0)then ! stupid way of doing this...but..
                 do kk=1,npartmax
                    if(idp(kk)<0)then
                       if(zp(kk).lt.Zcrit_pop3)then
                          x = (xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                          y = (xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                          z = (xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
                          if(  (abs(xp(kk,1)-x).le.r3pc).and.&
                          &    (abs(xp(kk,2)-y).le.r3pc).and.&
                          &    (abs(xp(kk,3)-z).le.r3pc)) nstar(i)=0
                       endif 
                    endif
                 end do
              endif
#endif

              ! Compute new stars local statistics
              mstar_tot=mstar_tot+nstar(i)*mstar1
              if(nstar(i)>0)then
                 ntot=ntot+1
                 if(sf_log_properties) then
399                  format('SF z= ',f9.5,1x,' nH= ',f7.3, ' T= ',f7.3  ,' N= ',I5,' eps_ff= ',f9.5,' alpha= ',f9.5)
                     write(*,399) 1./aexp-1,log10(d*scale_nH), log10(Tkstar) ,nstar(i),sfr_ff(i),alpha_fk(i)
                 endif
                 if(SFdiagnostics) then
666                  format('SF z= ',f9.5,1x,' nH= ',f7.3,' N= ',I5,' eps_ff= ',f9.5,' alpha= ',f9.5)
                     write(SF_filenr,666) 1./aexp-1,log10(d*scale_nH),nstar(i),sfr_ff(i),alpha_fk(i)
                 endif
                 if(f_w>0)ndebris_tot=ndebris_tot+1
              else
                 nstar(i)=0
              endif
           endif
        enddo
        ! Store nstar in array flag2
        do i=1,ngrid
           flag2(ind_cell(i))=nstar(i)
        end do
     end do
  end do

  !---------------------------------
  ! Check for free particle memory
  !---------------------------------
  ok_free=(numbp_free-ntot-ndebris_tot)>=0
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(numbp_free,numbp_free_tot,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  numbp_free_tot=numbp_free
#endif
  if(.not. ok_free)then
     write(*,*)'No more free memory for particles'
     write(*,*)'Increase npartmax'
#ifndef WITHOUTMPI
    call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
    stop
#endif
  end if

  !---------------------------------
  ! Compute global stars statistics
  !---------------------------------
#ifndef WITHOUTMPI
  mlost=mstar_lost; mtot=mstar_tot
  call MPI_ALLREDUCE(ntot,ntot_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mtot,mtot_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mlost,mlost_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  ntot_all=ntot
  mtot_all=mstar_tot
  mlost_all=mstar_lost
#endif
  ntot_star_cpu=0; ntot_star_all=0
  ntot_star_cpu(myid)=ntot
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(ntot_star_cpu,ntot_star_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  ntot_star_cpu(1)=ntot_star_all(1)
#endif
  do icpu=2,ncpu
     ntot_star_cpu(icpu)=ntot_star_cpu(icpu-1)+ntot_star_all(icpu)
  end do
  nstar_tot=nstar_tot+ntot_all
  !shubham
  mstar_total_all = mtot_all
  !shubham
  if(myid==1)then
     if(ntot_all.gt.0)then
        write(*,'(" Level=",I6," New star=",I6," Tot=",I10," Mass=",1PE10.3," Lost=",0PF5.1,"%")')&
             & ilevel,ntot_all,nstar_tot,mtot_all*scale_msun,mlost_all/(mlost_all+mtot_all)*100.
     endif
  end if

  !------------------------------
  ! Create new star particles
  !------------------------------
  ! Starting identity number
  if(myid==1)then
     index_star=nstar_tot-ntot_all
  else
     index_star=nstar_tot-ntot_all+ntot_star_cpu(myid-1)
  end if

  ! Loop over grids
  idx3_new = 1  ! index for PopIII (ind_cell_new)
  idx3   = 1  ! index for PopIII (ind_cell)
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! Loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

        ! Flag cells with at least one new star
        do i=1,ngrid
           ok(i)=flag2(ind_cell(i))>0
        end do

        ! Gather new star arrays
        nnew=0
        do i=1,ngrid
           if (ok(i))then
              nnew=nnew+1
              ind_grid_new(nnew)=ind_grid(i)
              ind_cell_new(nnew)=ind_cell(i)
           end if
        end do

        ! Update linked list for stars
        call remove_free(ind_part,nnew)
        call add_list(ind_part,ind_grid_new,ok_new,nnew)

        ! Update linked list for debris
        if(f_w>0)then
           call remove_free(ind_debris,nnew)
           call add_list(ind_debris,ind_grid_new,ok_new,nnew)
        endif

        ! Calculate new star particle and modify gas density
        do i=1,nnew
           index_star=index_star+1

           ! Get gas variables
           n=flag2(ind_cell_new(i))
           d=uold(ind_cell_new(i),1)
           u=uold(ind_cell_new(i),2)
           v=uold(ind_cell_new(i),3)
           w=uold(ind_cell_new(i),4)
           x=(xg(ind_grid_new(i),1)+xc(ind,1)-skip_loc(1))*scale
           y=(xg(ind_grid_new(i),2)+xc(ind,2)-skip_loc(2))*scale
           z=(xg(ind_grid_new(i),3)+xc(ind,3)-skip_loc(3))*scale
           tg=uold(ind_cell_new(i),5)*(gamma-1)*scale_T2
           if(metal) then
              zg=uold(ind_cell_new(i),imetal)
              if(is_oxygen) zg_ox=uold(ind_cell_new(i),imetal+1)  ! oxygen
#ifdef NCHEM
              do iche=1,nchem
                 chem1(iche) = uold(ind_cell_new(i),ichem+iche-1)
              enddo
#endif
           endif
           mstar1 = mstar
#ifdef POP3
           if(pop3)then
              if(ind_cell_new(i).eq.idx3icell(idx3_new))then
                 mstar1 = mass_pop3(idx3_new)
                 idx3_new = idx3_new + 1
              endif
           endif
#endif

           ! Set star particle variables
           tp(ind_part(i)) = birth_epoch     ! Birth epoch
           mp(ind_part(i)) = n*mstar1        ! Mass
           mp0(ind_part(i)) = mp(ind_part(i)) ! Initial Mass
           levelp(ind_part(i)) = ilevel      ! Level
           idp(ind_part(i))    = -index_star ! Star identity
           typep(ind_part(i))%family = FAM_STAR
           typep(ind_part(i))%tag = 0
           xp(ind_part(i),1)   = x
           xp(ind_part(i),2)   = y
           xp(ind_part(i),3)   = z
           vp(ind_part(i),1)   = u
           vp(ind_part(i),2)   = v
           vp(ind_part(i),3)   = w

           ! Random kick -------------------
           if(sf_kick_kms .gt. 0d0) then
              call ranf(localseed,RandNum)
              v_kick(1) = (Randnum - 0.5) * 2.
              call ranf(localseed,RandNum)
              v_kick(2) = (Randnum - 0.5) * 2.
              call ranf(localseed,RandNum)
              v_kick(3) = (Randnum - 0.5) * 2.
              v_kick_mag = sqrt(sum((v_kick(1:3))**2))
              if(v_kick_mag .gt. 0d0) then
                 v_kick = v_kick / v_kick_mag ! Normalize to unit length
              endif
              call ranf(localseed,RandNum)
              if(sf_kick_pow < 0) then
                 v_kick_mag = (0.1**(sf_kick_pow+1)+((1.1**(sf_kick_pow+1))-(0.1**(sf_kick_pow+1)))*(Randnum+0.1))**(1./(sf_kick_pow+1))
                 v_kick_mag = v_kick_mag/((0.1**(sf_kick_pow+1)+((1.1**(sf_kick_pow+1))-(0.1**(sf_kick_pow+1)))*1.1)**(1./(sf_kick_pow+1)))
                 v_kick_mag = v_kick_mag * sf_kick_kms*1d5 / scale_v
              else
               ! 0-sf_kick km/s kick:
                 v_kick_mag = Randnum * sf_kick_kms*1d5 / scale_v
              endif
              v_kick = v_kick * v_kick_mag
              vp(ind_part(i),1)   = vp(ind_part(i),1) + v_kick(1)
              vp(ind_part(i),2)   = vp(ind_part(i),2) + v_kick(2)
              vp(ind_part(i),3)   = vp(ind_part(i),3) + v_kick(3)
           endif
           ! End random kick ---------------------------------------------

           if (metal) then
              zp(ind_part(i)) = zg  ! Initial star metallicity
              if(is_oxygen) zp_ox(ind_part(i))=zg_ox
#ifdef NCHEM
              do iche=1,nchem
                 chp(ind_part(i),iche) = chem1(iche)  ! Initial chemical abudance
              enddo
#endif
           endif
#ifdef RT
#ifdef NTRACEGROUPS
           dumpid = 0
           if (rt_tracemassbins(1).gt.-1.d0) then
              call get_tracer_group_from_halo(xp(ind_part(i),1:ndim),dumpid,minind)
           endif
           ptracegroup(ind_part(i)) = dumpid
#endif
#endif

           ! Set GMC particle variables
           if(f_w>0)then
              ! Compute GMC mass without more than 50% of gas depletion
              mdebris=min(f_w*n*mstar1,0.5d0*(d*vol_loc-n*mstar1))
              ! Add supernova ejecta
              mdebris=mdebris+eta_sn*n*mstar1
              ! Remove ejecta from the long lived star mass
              mp(ind_part(i))=n*mstar1-eta_sn*n*mstar1
              mp0(ind_part(i))=mp(ind_part(i))
              ! Set GMC particle variables
              tp(ind_debris(i))=birth_epoch  ! Birth epoch
              mp(ind_debris(i))=mdebris      ! Mass
              levelp(ind_debris(i))=ilevel   ! Level
              idp(ind_debris(i))=-n          ! Number of individual stars
              typep(ind_debris(i))%family = FAM_DEBRIS
              typep(ind_debris(i))%tag = 0
              xp(ind_debris(i),1)=x
              xp(ind_debris(i),2)=y
              xp(ind_debris(i),3)=z
              vp(ind_debris(i),1)=u
              vp(ind_debris(i),2)=v
              vp(ind_debris(i),3)=w
              ! GMC metallicity + yield from ejecta
              if(metal)zp(ind_debris(i))=zg+eta_sn*yield*(1-zg)*n*mstar1/mdebris
           endif

           if(sf_log_properties) then
              write(ilun,'(I10)',advance='no') 0
              write(ilun,'(2I10,E24.12)',advance='no') idp(ind_part(i)),ilevel,mp(ind_part(i))
              do idim=1,ndim
                 write(ilun,'(E24.12)',advance='no') xp(ind_part(i),idim)
              enddo
              do idim=1,ndim
                 write(ilun,'(E24.12)',advance='no') vp(ind_part(i),idim)
              enddo
              write(ilun,'(E24.12)',advance='no') uold(ind_cell_new(i),1)
              do ivar=2,nvar
                 if(ivar.eq.ndim+2)then
                    ! Temperature
                    uvar=(gamma-1.0d0)*(uold(ind_cell_new(i),ndim+2))*scale_T2
                 else
                    uvar=uold(ind_cell_new(i),ivar)
                 endif
                 write(ilun,'(E24.12)',advance='no') uvar
              enddo
              write(ilun,'(I10)',advance='no') typep(ind_part(i))%tag
              write(ilun,'(A1)') ' '
           endif


        end do
        ! End loop over new star particles

        ! Modify gas density according to mass depletion
        do i=1,ngrid
          if(flag2(ind_cell(i))>0)then
            n=flag2(ind_cell(i))
            d=uold(ind_cell(i),1)
            ! for PopIII stars
            mstar1 = mstar

#ifdef POP3
            if(pop3)then
              if(ind_cell(i).eq.idx3icell(idx3))then
                mstar1 = mass_pop3(idx3)
                idx3 = idx3 + 1
              endif
            endif
#endif
            dstar = mstar1/vol_loc
            uold(ind_cell(i),1)=max(d-n*dstar*(1.0d0+f_w),0.5d0*(d-n*dstar))
          endif
        end do

        !tracer    
#ifdef MC_tracer
        if(MC_tracer) then  
          dstar=mstar/vol_loc
          do i=1,nnew
            n=flag2(ind_cell_new(i))
            d=uold(ind_cell_new(i),1)
            if(mechanical_feedback==0) then
              delta_m_over_m = min(n*dstar*(1.0+f_w), 0.5_dp*d) / d
              !uold(ind_cell_new(i),1)=max(d-n*dstar*(1.0+f_w), 0.5*d)
            else
              delta_m_over_m = n*dstar / d
              !uold(ind_cell_new(i),1)=d-n*dstar
            end if

            ! Loop over particles in grid
            ipart = headp(ind_grid_new(i))

            ! Get star location
            x=xp(ind_part(i), 1)
            y=xp(ind_part(i), 2)
            z=xp(ind_part(i), 3)

            do ip = 1, numbp(ind_grid_new(i))
              ! Keep only *UNMOVED* tracer particles attached to cell containing new stars
              move_tracer = (is_gas_tracer(typep(ipart)) .and. partp(ipart) == ind_cell_new(i))
              if (move_tracer) then
                nattach = nattach + 1
                tok(nattach)          = move_tracer
                itracer(nattach)      = ipart
                istar_tracer(nattach) = ind_part(i)
                proba(nattach)        = delta_m_over_m
                ! The star tracers are moved in the move fine routine
                ! This here is just to give a non exact yet approximative location
                ! for the tracer
                xstar(nattach, 1)     = x
                xstar(nattach, 2)     = y
                xstar(nattach, 3)     = z
              end if

              if (nattach == nvector) then
                call tracer2star(itracer, proba, xstar, istar_tracer, nattach)
                nattach = 0
                tok = .false.
                itracer = 0
                istar_tracer = 0
                proba = 0
                xstar = 0
              end if
              ipart = nextp(ipart)
            end do      
          end do !nnew
        endif 
#endif 
        !tracer 
     end do
     ! End loop over cells
  end do
  ! End loop over grids

  !tracer 
#ifdef MC_tracer
  if(MC_tracer) then 
    ! Empty tracer part cache
    if (MC_tracer .and. nattach > 0) then
      call tracer2star(itracer, proba, xstar, istar_tracer, nattach)
      nattach = 0
    end if
  endif 
#endif 
  !tracer 

  !---------------------------------------------------------
  ! Convert hydro variables back to conservative variables
  !---------------------------------------------------------
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        do i=1,ngrid
           d=uold(ind_cell(i),1)
           u=uold(ind_cell(i),2)
           v=uold(ind_cell(i),3)
           w=uold(ind_cell(i),4)
           e=uold(ind_cell(i),5)*d
#ifdef SOLVERmhd
           bx1=uold(ind_cell(i),6)
           by1=uold(ind_cell(i),7)
           bz1=uold(ind_cell(i),8)
           bx2=uold(ind_cell(i),nvar+1)
           by2=uold(ind_cell(i),nvar+2)
           bz2=uold(ind_cell(i),nvar+3)
           e=e+0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)
#endif
           e=e+0.5d0*d*(u**2+v**2+w**2)
#if NENER>0
           do irad=0,nener-1
              e=e+uold(ind_cell(i),inener+irad)
           end do
#endif
           uold(ind_cell(i),1)=d
           uold(ind_cell(i),2)=d*u
           uold(ind_cell(i),3)=d*v
           uold(ind_cell(i),4)=d*w
           uold(ind_cell(i),5)=e
        end do
        do ivar=imetal,nvar
           do i=1,ngrid
              d=uold(ind_cell(i),1)
              w=uold(ind_cell(i),ivar)
              uold(ind_cell(i),ivar)=d*w
           end do
        end do
     end do
  end do

  if(sf_log_properties) close(ilun)

contains
  
!################################################################
!################################################################
!################################################################
!################################################################
subroutine get_sfeff_virial()  
  !-----------------------------------------------------------------
  ! Check if star formation allowed in given cells and, if so,
  ! calculate and return local star formation efficiency
  !-----------------------------------------------------------------
  implicit none
  ! TODO: when f2008 is obligatory - remove this and replace erfc_pre_f08 below by
  ! the f2008 intrinsic erfc() function:
  real(dp) erfc_pre_f08
  !-----------------------------------------------------------------
  do i=1,ngrid
     ! if cell is a leaf cell
     if (ok(i)) then
        ! Subgrid turbulence decay
        if(sf_tdiss.gt.0d0) then
           if(sf_compressive) then
              tdec = sf_tdiss*dx_loc/sqrt(uold(ind_cell(i),ivirial1)+uold(ind_cell(i),ivirial2))
              if(uold(ind_cell(i),ivirial1).gt.0d0) uold(ind_cell(i),ivirial1) = uold(ind_cell(i),ivirial1)*exp(-dtold(ilevel)/tdec)
              if(uold(ind_cell(i),ivirial2).gt.0d0) uold(ind_cell(i),ivirial2) = uold(ind_cell(i),ivirial2)*exp(-dtold(ilevel)/tdec)
           else
              tdec = sf_tdiss*dx_loc/sqrt(uold(ind_cell(i),ivirial1))
              if(uold(ind_cell(i),ivirial1).gt.0d0) uold(ind_cell(i),ivirial1) = uold(ind_cell(i),ivirial1)*exp(-dtold(ilevel)/tdec)
           endif
        endif
        d         = uold(ind_cell(i),1)
        ! Compute temperature in K/mu
        T2        = (gamma-1.0d0)*uold(ind_cell(i),5)*scale_T2
        ! Correct from polytrope
        T_poly    = T2_star*(uold(ind_cell(i),1)*scale_nH/nISM)**(g_star-1.0d0)
        T2        = T2-T_poly
        ! Compute sound speed squared
        cs2       = (gamma-1.0d0)*uold(ind_cell(i),5)
#if NENER>0
        !! Add radiation pressure to sound speed calculation
        !do irad = 0,nener-1
        !   cs2 = cs2 + uold(ind_cell(i),inener+irad) / d * (gamma_rad(irad+1)-1.0)
        !end do
#endif
        ! prevent numerical crash due to negative temperature
        cs2       = max(cs2,smallc**2)
        ! Correct from polytrope
        cs2_poly  = (T2_star/scale_T2)*(uold(ind_cell(i),1)*scale_nH/nISM)**(g_star-1.0d0)
        cs2       = cs2-cs2_poly
        ! We need to estimate the norm of the gradient of the velocity field in the cell (tensor of 2nd rank)
        ! i.e. || A ||^2 = trace( A A^T) where A = grad vec(v) is the tensor.
        ! So construct values of velocity field on the 6 faces of the cell using simple linear interpolation
        ! from neighbouring cell values and differentiate.
        ! Get neighbor cells if they exist, otherwise use straight injection from local cell
        ncell = 1 ! we just want the neighbors of that cell
        ind_cell2(1) = ind_cell(i)
        call getnbor(ind_cell2,ind_nbor,ncell,ilevel)
        d1           = uold(ind_nbor(1,1),1) ; d2 = uold(ind_nbor(1,2),1) ; d3 = uold(ind_nbor(1,3),1)
        d4           = uold(ind_nbor(1,4),1) ; d5 = uold(ind_nbor(1,5),1) ; d6 = uold(ind_nbor(1,6),1)
        sigma2       = 0d0 ; sigma2_comp = 0d0 ; sigma2_sole = 0d0
        trgv         = 0d0 ; divv = 0d0 ; curlva = 0d0 ; curlvb = 0d0 ; curlvc = 0d0
        flong        = 0d0
        !!!!!!!!!!!!!!!!!!
        ! Divergence terms
        !!!!!!!!!!!!!!!!!!
        ul        = (d2*uold(ind_nbor(1,2),2) + d*uold(ind_cell(i),2))/(d2+d)
        ur        = (d1*uold(ind_nbor(1,1),2) + d*uold(ind_cell(i),2))/(d1+d)
        if(sf_model.le.2) then
           fl     = (d2*f(ind_nbor(1,2),1)    + d*f(ind_cell(i),1))/(d2+d)
           fr     = (d1*f(ind_nbor(1,1),1)    + d*f(ind_cell(i),1))/(d1+d)
           flong  = flong+max((d2+d)/2*ul*fl-(d1+d)/2*ur*fr,0d0)
        endif
        sigma2_comp = sigma2_comp + (ur-ul)**2
        divv      = divv + (ur-ul)
        ul        = (d4*uold(ind_nbor(1,4),3) + d*uold(ind_cell(i),3))/(d4+d)
        ur        = (d3*uold(ind_nbor(1,3),3) + d*uold(ind_cell(i),3))/(d3+d)
        if(sf_model.le.2) then
           fl     = (d4*f(ind_nbor(1,4),2)    + d*f(ind_cell(i),2))/(d4+d)
           fr     = (d3*f(ind_nbor(1,3),2)    + d*f(ind_cell(i),2))/(d3+d)
           flong  = flong+max((d4+d)/2*ul*fl-(d3+d)/2*ur*fr,0d0)
        endif
        sigma2_comp = sigma2_comp + (ur-ul)**2
        divv      = divv + (ur-ul)
        ul        = (d6*uold(ind_nbor(1,6),4) + d*uold(ind_cell(i),4))/(d6+d)
        ur        = (d5*uold(ind_nbor(1,5),4) + d*uold(ind_cell(i),4))/(d5+d)
        if(sf_model.le.2) then
           fl     = (d6*f(ind_nbor(1,6),3)    + d*f(ind_cell(i),3))/(d6+d)
           fr     = (d5*f(ind_nbor(1,5),3)    + d*f(ind_cell(i),3))/(d5+d)
           flong  = flong+max((d6+d)/2*ul*fl-(d5+d)/2*ur*fr,0d0)
        endif
        sigma2_comp = sigma2_comp + (ur-ul)**2
        divv      = divv + (ur-ul)
        ftot      = flong
        !!!!!!!!!!!!
        ! Curl terms
        !!!!!!!!!!!!
        ul        = (d6*uold(ind_nbor(1,6),3) + d*uold(ind_cell(i),3))/(d6+d)
        ur        = (d5*uold(ind_nbor(1,5),3) + d*uold(ind_cell(i),3))/(d5+d)
        if(sf_model.le.2) then
           fl     = (d6*f(ind_nbor(1,6),2)    + d*f(ind_cell(i),2))/(d6+d)
           fr     = (d5*f(ind_nbor(1,5),2)    + d*f(ind_cell(i),2))/(d5+d)
           ftot   = ftot+abs((d6+d)/2*ul*fl-(d5+d)/2*ur*fr)
        endif
        sigma2_sole = sigma2_sole + (ur-ul)**2
        curlva    = curlva-(ur-ul)
        ul        = (d4*uold(ind_nbor(1,4),4) + d*uold(ind_cell(i),4))/(d4+d)
        ur        = (d3*uold(ind_nbor(1,3),4) + d*uold(ind_cell(i),4))/(d3+d)
        if(sf_model.le.2) then
           fl     = (d4*f(ind_nbor(1,4),3)    + d*f(ind_cell(i),3))/(d4+d)
           fr     = (d3*f(ind_nbor(1,3),3)    + d*f(ind_cell(i),3))/(d3+d)
           ftot   = ftot+abs((d4+d)/2*ul*fl-(d3+d)/2*ur*fr)
        endif
        sigma2_sole = sigma2_sole + (ur-ul)**2
        curlva    = (curlva + (ur-ul))
        ul        = (d6*uold(ind_nbor(1,6),2) + d*uold(ind_cell(i),2))/(d6+d)
        ur        = (d5*uold(ind_nbor(1,5),2) + d*uold(ind_cell(i),2))/(d5+d)
        if(sf_model.le.2) then
           fl     = (d6*f(ind_nbor(1,6),1)    + d*f(ind_cell(i),1))/(d6+d)
           fr     = (d5*f(ind_nbor(1,5),1)    + d*f(ind_cell(i),1))/(d5+d)
           ftot   = ftot+abs((d6+d)/2*ul*fl-(d5+d)/2*ur*fr)
        endif
        sigma2_sole = sigma2_sole + (ur-ul)**2
        curlvb    = curlvb+(ur-ul)
        ul        = (d2*uold(ind_nbor(1,2),4) + d*uold(ind_cell(i),4))/(d2+d)
        ur        = (d1*uold(ind_nbor(1,1),4) + d*uold(ind_cell(i),4))/(d1+d)
        if(sf_model.le.2) then
           fl     = (d2*f(ind_nbor(1,2),3)    + d*f(ind_cell(i),3))/(d2+d)
           fr     = (d1*f(ind_nbor(1,1),3)    + d*f(ind_cell(i),3))/(d1+d)
           ftot   = ftot+abs((d2+d)/2*ul*fl-(d1+d)/2*ur*fr)
        endif
        sigma2_sole = sigma2_sole + (ur-ul)**2
        curlvb    = (curlvb - (ur-ul))
        ul        = (d4*uold(ind_nbor(1,4),2) + d*uold(ind_cell(i),2))/(d4+d)
        ur        = (d3*uold(ind_nbor(1,3),2) + d*uold(ind_cell(i),2))/(d3+d)
        if(sf_model.le.2) then
           fl     = (d4*f(ind_nbor(1,4),1)    + d*f(ind_cell(i),1))/(d4+d)
           fr     = (d3*f(ind_nbor(1,3),1)    + d*f(ind_cell(i),1))/(d3+d)
           ftot   = ftot+abs((d4+d)/2*ul*fl-(d3+d)/2*ur*fr)
        endif
        sigma2_sole = sigma2_sole + (ur-ul)**2
        curlvc    = curlvc-(ur-ul)
        ul        = (d2*uold(ind_nbor(1,2),3) + d*uold(ind_cell(i),3))/(d2+d)
        ur        = (d1*uold(ind_nbor(1,1),3) + d*uold(ind_cell(i),3))/(d1+d)
        if(sf_model.le.2) then
           fl     = (d2*f(ind_nbor(1,2),2)    + d*f(ind_cell(i),2))/(d2+d)
           fr     = (d1*f(ind_nbor(1,1),2)    + d*f(ind_cell(i),2))/(d1+d)
           ftot   = ftot+abs((d2+d)/2*ul*fl-(d1+d)/2*ur*fr)
           pcomp  = flong/ftot
        endif
        sigma2_sole = sigma2_sole + (ur-ul)**2
        curlvc    = (curlvc + (ur-ul))
        sigma2    = sigma2_comp+sigma2_sole
        ! Trace of gradient velocity tensor
        trgv      = sigma2/dx_loc**2
        ! Velocity vector divergence
        divv      = divv/dx_loc
        ! Velocity vector curl
        curlv     = (curlva+curlvb+curlvc)/dx_loc
        divv2     = divv**2
        curlv2    = curlv**2
        ! Advect unresolved turbulence if a decay time is defined
        if(sf_tdiss.gt.0d0) then
           if(sf_compressive)then
              uold(ind_cell(i),ivirial1) = max(uold(ind_cell(i),ivirial1),0d0)+sigma2_comp
              uold(ind_cell(i),ivirial2) = max(uold(ind_cell(i),ivirial2),0d0)+sigma2_sole
              sigma2_comp = uold(ind_cell(i),ivirial1)
              sigma2_sole = uold(ind_cell(i),ivirial2)
              sigma2      = sigma2_sole+sigma2_comp
           else
              uold(ind_cell(i),ivirial1) = max(uold(ind_cell(i),ivirial1),0d0)+sigma2
              sigma2 = uold(ind_cell(i),ivirial1)
           endif
        else
           if(sf_compressive)then
              uold(ind_cell(i),ivirial1) = sigma2_comp
              uold(ind_cell(i),ivirial2) = sigma2_sole
           else
              uold(ind_cell(i),ivirial1) = sigma2
           endif
        endif
        ! Density criterion
        if(d<=d0) ok(i)=.false.
        if(ok(i)) then
           SELECT CASE (sf_model)
              ! Classical density threshold
           CASE (0)
              sfr_ff(i) = eps_star
              ! Multi-ff KM model
           CASE (1)
              ! Virial parameter
              alpha0    = (5.0d0*sigma2)/(pi*factG*d*dx_loc**2)
              ! Turbulent forcing parameter (Federrath 2008 & 2010)
              if(pcomp*ndim-1.0d0 == 0d0) then
                 zeta   = 0.5d0
              else
                 zeta   = ((pcomp-1.0d0)+sqrt((pcomp**2-pcomp)*(1.0d0-ndim)))/(pcomp*ndim-1.0d0)
              endif
              b_turb    = 1.0d0+(1.0d0/ndim-1.0d0)*zeta
#ifdef SOLVERmhd
              ! Best fit values to the Multi-ff KM model (MHD)
              phi_t     = 0.46d0
              phi_x     = 0.17d0
              A         = 0.5d0*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
              B         = 0.5d0*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
              C         = 0.5d0*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))
              emag      = 0.5d0*(A**2+B**2+C**2)
              beta      = uold(ind_cell(i),5)*d/max(emag,smallc**2*smallr)
              sigs      = log(1.0d0+(b_turb**2)*(sigma2/cs2)*beta/(beta+1.0d0))
              scrit     = log(((pi**2)/5)*(phi_x**2)*alpha0*(sigma2/cs2)/(1.0d0+1.0d0/beta))
#else
              ! Best fit values to the Multi-ff KM model (Hydro)
              phi_t     = 0.49d0
              phi_x     = 0.19d0
              sigs      = log(1.0d0+(b_turb**2)*(sigma2/cs2))
              scrit     = log(((pi**2)/5)*(phi_x**2)*alpha0*(sigma2/cs2))
#endif
              sfr_ff(i) = (eps_star_loc*phi_t/2.0d0)*exp(3.0d0/8.0d0*sigs)*(2.0d0-erfc_pre_f08((sigs-scrit)/sqrt(2.0d0*sigs)))
              ! Multi-ff PN model
           CASE (2)
              ! Virial parameter
              alpha0    = (5.0d0*sigma2)/(pi*factG*d*dx_loc**2)
              ! Turbulent forcing parameter (Federrath 2008 & 2010)
              if(pcomp*ndim-1.0d0 == 0d0) then
                 zeta   = 0.5d0
              else
                 zeta   = ((pcomp-1.0d0)+sqrt((pcomp**2-pcomp)*(1.0d0-ndim)))/(pcomp*ndim-1.0d0)
              endif
              b_turb    = 1.0d0+(1.0d0/ndim-1.0d0)*zeta
#ifdef SOLVERmhd
              ! Best fit values to the Multi-ff PN model (MHD)
              phi_t     = 0.47d0
              theta     = 1.00d0
              A         = (uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))/2
              B         = (uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))/2
              C         = (uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))/2
              emag      = (A**2+B**2+C**2)/2
              beta      = uold(ind_cell(i),5)*d/max(emag,smallc**2*smallr)
              fbeta     = ((1+0.925d0*beta**(-3.0d0/2.0d0))**(2.0d0/3.0d0))/((1.0d0+1.0d0/beta)**2)
              sigs      = log(1.0d0+(b_turb**2)*(sigma2/cs2)*beta/(beta+1.0d0))
              scrit     = log(0.067d0/(theta**2)*alpha0*(sigma2/cs2)*fbeta)
#else
              ! Best fit values to the Multi-ff PN model (Hydro)
              phi_t     = 0.49d0
              theta     = 0.97d0
              sigs      = log(1.0d0+(b_turb**2)*(sigma2/cs2))
              scrit     = log(0.067d0/(theta**2)*alpha0*(sigma2/cs2))
#endif
              sfr_ff(i) = (eps_star_loc*phi_t/2)*exp(3.0d0/8.0d0*sigs)*(2.0d0-erfc_pre_f08((sigs-scrit)/sqrt(2.0d0*sigs)))
              ! Virial criterion simple model
           CASE (3)
              ! Laplacian rho
              lapld     =       ((d1+d)/2)-2*d+((d2+d)/2)
              lapld     = lapld+((d3+d)/2)-2*d+((d4+d)/2)
              lapld     = lapld+((d5+d)/2)-2*d+((d6+d)/2)
              lapld     = lapld/(dx_loc/2)**2
              alpha0    = (trgv-cs2*lapld/d)/(4*pi*factG*d)
              if(alpha0<1.0.and.lapld<0.0) then
                 sfr_ff(i) = eps_star
              else
                 sfr_ff(i) = 0
                 ok(i)     = .false.
              endif
              ! Padoan 2012 "a simple SF law"
           CASE (4)
              ! Feedback efficiency
              t_dyn     = dx_loc/(2*sqrt(sigma2+cs2))
              t_ff      = 0.5427d0*sqrt(1/(factG*max(d,smallr)))
              sfr_ff(i) = eps_star*exp(-1.6d0*t_ff/t_dyn)
              ! Hopkins 2013
           CASE (5)
              alpha0    = 0.5d0*(divv2+curlv2)/(factG*d)
              if(alpha0<1.0) then
                 sfr_ff(i) = eps_star
              else
                 sfr_ff(i) = 0
                 ok(i)     = .false.
              endif
           END SELECT
        endif
     endif
  end do
  
end subroutine get_sfeff_virial

!################################################################
!################################################################
!################################################################
!################################################################
subroutine get_sfeff_non_virial()
  !-----------------------------------------------------------------
  ! Check if star formation allowed in given cells and, if so,
  ! calculate and return local star formation efficiency
  !-----------------------------------------------------------------
  implicit none
  real(dp),dimension(0:twondim) ::darr,uarr,varr,warr
  real(dp)::lamjt
  real(dp)::px_div,py_div,pz_div,Jx,Jy,Jz,rho_local
  real(dp) :: uavg,vavg,wavg,dtot
  integer::idir
  logical::isConvergent
  ! TODO: when f2008 is obligatory - remove this and replace erfc_pre_f08 below by
  ! the f2008 intrinsic erfc() function:
  real(dp) erfc_pre_f08
  !-----------------------------------------------------------------
  if (TRIM(star_maker)=='density') then
     do i=1,ngrid
        if (ok(i))then
           d=uold(ind_cell(i),1)
           ! Density criterion
           if(d<=d0)then
              ok(i)=.false.
           ! Temperature criterion
           else
              T2=uold(ind_cell(i),5)*scale_T2*(gamma-1.0d0)
              nH=max(uold(ind_cell(i),1),smallr)*scale_nH
              T_poly=T2_star*(nH/nISM)**(g_star-1.0d0)
              T2=T2-T_poly
              if(T2>T2thres_SF)then
                 ok(i)=.false.
              else
                 sfr_ff(i)=eps_star
              endif
           endif
        endif
     end do
     
  else if(TRIM(star_maker)=='federrath')then
     ! Enforce turbulence criterion + efficiency following Federrath & Klessen 2012
     do i=1,ngrid
        ! if cell is a leaf cell
        if (ok(i)) then 
           d = uold(ind_cell(i),1)
           ! set density threshold to avoid to compute tensor trace too often but otherwise not needed
           if (d <= d_gmc) then
              ok(i) = .false.
           else
              ! We need to estimate the norm of the gradient of the velocity field in the cell (tensor of 2nd rank)
              ! i.e. || A ||^2 = trace( A A^T) where A = grad vec(v) is the tensor. 
              ! So construct values of velocity field on the 6 faces of the cell using simple linear interpolation 
              ! from neighbouring cell values and differentiate. 
              ! Get neighbor cells if they exist, otherwise use straight injection from local cell
              ncell = 1 ! we just want the neighbors of that cell
              call getnbor(ind_cell(i),ind_nbor,ncell,ilevel)
              d1    = uold(ind_nbor(1,1),1)     ; d4    = uold(ind_nbor(1,4),1) 
              d2    = uold(ind_nbor(1,2),1)     ; d5    = uold(ind_nbor(1,5),1) 
              d3    = uold(ind_nbor(1,3),1)     ; d6    = uold(ind_nbor(1,6),1)  
              ul    = (d2*uold(ind_nbor(1,2),2) + d*uold(ind_cell(i),2))/(d2+d)
              ur    = (d1*uold(ind_nbor(1,1),2) + d*uold(ind_cell(i),2))/(d1+d)
              trgv  = (ur-ul)**2
              ul    = (d4*uold(ind_nbor(1,4),3) + d*uold(ind_cell(i),3))/(d4+d)
              ur    = (d3*uold(ind_nbor(1,3),3) + d*uold(ind_cell(i),3))/(d3+d)
              trgv  = trgv + (ur-ul)**2
              ul    = (d6*uold(ind_nbor(1,6),4) + d*uold(ind_cell(i),4))/(d6+d)
              ur    = (d5*uold(ind_nbor(1,5),4) + d*uold(ind_cell(i),4))/(d5+d)
              trgv  = trgv + (ur-ul)**2
              ul    = (d6*uold(ind_nbor(1,6),3) + d*uold(ind_cell(i),3))/(d6+d)
              ur    = (d5*uold(ind_nbor(1,5),3) + d*uold(ind_cell(i),3))/(d5+d)
              trgv  = trgv + (ur-ul)**2
              ul    = (d4*uold(ind_nbor(1,4),4) + d*uold(ind_cell(i),4))/(d4+d)
              ur    = (d3*uold(ind_nbor(1,3),4) + d*uold(ind_cell(i),4))/(d3+d)
              trgv  = trgv + (ur-ul)**2
              ul    = (d6*uold(ind_nbor(1,6),2) + d*uold(ind_cell(i),2))/(d6+d)
              ur    = (d5*uold(ind_nbor(1,5),2) + d*uold(ind_cell(i),2))/(d5+d)
              trgv  = trgv + (ur-ul)**2
              ul    = (d2*uold(ind_nbor(1,2),4) + d*uold(ind_cell(i),4))/(d2+d)
              ur    = (d1*uold(ind_nbor(1,1),4) + d*uold(ind_cell(i),4))/(d1+d)
              trgv  = trgv + (ur-ul)**2
              ul    = (d4*uold(ind_nbor(1,4),2) + d*uold(ind_cell(i),2))/(d4+d)
              ur    = (d3*uold(ind_nbor(1,3),2) + d*uold(ind_cell(i),2))/(d3+d)
              trgv  = trgv + (ur-ul)**2
              ul    = (d2*uold(ind_nbor(1,2),3) + d*uold(ind_cell(i),3))/(d2+d)
              ur    = (d1*uold(ind_nbor(1,1),3) + d*uold(ind_cell(i),3))/(d1+d)
              trgv  = trgv + (ur-ul)**2
              ! now compute sound speed squared = cs^2 
              cs2  = uold(ind_cell(i),5)*(gamma -1.0)
#if NENER>0
              !! Add non-thermal pressure to sound speed calculation
              !do irad = 0,nener-1
              !   cs2 = cs2 + uold(ind_cell(i),inener+irad) &
              !        / d * (gamma_rad(irad+1)-1.0)
              !end do
#endif
              ! prevent numerical crash due to negative temperature
              cs2  = max(cs2,smallc**2)
              !rho_local = d
              rho_local = rho(ind_cell(i)) ! be careful with dark matter particles
              if(sf_lam.gt.0) then ! Jeans length criterion
                 ! Calculate "turbulent" Jeans length in cell units, lamjt 
                 ! (see e.g. Chandrasekhar 51, Bonazzola et al 87, Federrath & Klessen 2012 eq 36)
                 lamjt = (pi*trgv + sqrt(pi*pi*trgv*trgv + 36.0*pi*cs2*factG*rho_local*dx_loc**2))/(6.0*factG*rho_local*dx_loc**2)
                 if (lamjt > sf_lam) ok(i) = .false. ! Jeans length resolved: gas is stable 
              endif
              if(ok(i)) then ! Jeans length not resolved --> form stars to lower density and stabilise gas
                 ! corresponding virial parameter for homogeneous sphere <= 1.5 in turbulent dominated 
                 ! limit and <= 0.5 in pressure dominated limit (in good agreement with observations,
                 ! at least for massive (>= 10^5 M_sun) clouds see Kauffmann, Pillai, Goldsmith 2013, Fig.1)  
                 alpha0 = 5d0/(pi*factG*rho_local)*(trgv + cs2)/dx_loc**2
                 ! compute star formation efficiency per free-fall time (Federrath & Klessen 2012 eq 41)
                 phi_t = 0.57 ; theta = 0.33 ! best fit values of the Padoan & Nordlund multi-scale sf model to GMC simulation data 
                 sigs  = log(1.0+0.16*trgv/cs2) ! factor 0.16 is b^2 where b=0.4 for a mixture of turbulence forcing modes
                 scrit = log(0.067/theta**2*alpha0*trgv/cs2) ! best fit from Padoan & Nordlund MS model again 
                 sfr_ff(i) = eps_star_loc/2.0*phi_t*exp(3.0/8.0*sigs)*(2.0-erfc_pre_f08((sigs-scrit)/sqrt(2.0*sigs)))
                 alpha_fk(i)=alpha0
              endif
           endif
        endif
     end do
     
  else if(TRIM(star_maker)=='FK2')then
     ! Enforce turbulence criterion + efficiency following Federrath & Klessen 2012
     do i=1,ngrid
        ! if cell is a leaf cell and high density
        d = uold(ind_cell(i),1)
        ok(i) = ok(i).and.(d>=d_gmc)
        if (ok(i)) then 
           ! We need to estimate the norm of the gradient of the velocity field in the cell (tensor of 2nd rank)
           ! i.e. || A ||^2 = trace( A A^T) where A = grad vec(v) is the tensor. 
           ! So construct values of velocity field on the 6 faces of the cell using simple linear interpolation 
           ! from neighbouring cell values and differentiate. 
           ! Get neighbor cells if they exist, otherwise use straight injection from local cell
           ncell = 1 ! we just want the neighbors of that cell
           call getnbor(ind_cell(i),ind_nbor,ncell,ilevel)
           
           darr(0) = d
           do idir=1,6
              darr(idir) = uold(ind_nbor(1,idir),1)
           end do
           
           ! local density maxima
           if(ok(i)) then
              if(maxval(darr)>d) ok(i)=.false.
           endif
           
           ! converging flow check
           if(ok(i))then
              uarr(0) = uold(ind_cell(i),2)
              varr(0) = uold(ind_cell(i),3)
              warr(0) = uold(ind_cell(i),4)
              do idir=1,6
                 uarr(idir) = uold(ind_nbor(1,idir),2) 
                 varr(idir) = uold(ind_nbor(1,idir),3) 
                 warr(idir) = uold(ind_nbor(1,idir),4) 
              end do
              divv  = (uarr(2)*darr(2)-uarr(1)*darr(1)) &
                   & + (varr(4)*darr(4)-varr(3)*darr(3)) & 
                   & + (warr(6)*darr(6)-warr(5)*darr(5))
              if(divv>0) ok(i) = .false.  ! diverging flow
           endif
           
           if(ok(i)) then
              !average velocity
              dtot  = sum(darr)
              uavg  = sum(darr*uarr)/dtot
              vavg  = sum(darr*varr)/dtot
              wavg  = sum(darr*warr)/dtot
              
              ! subtract the mean velocity field
              uarr(:) = uarr(:) - uavg
              varr(:) = varr(:) - vavg
              warr(:) = warr(:) - wavg
              
              ! subtract the symmetric divergence field                    
              ! ex)  (---->,<--): only subtract (-->,<--): result (-->,0) 
              ! ex)  (<----,-->): only subtract (<--,-->): result (<--,0)
              px_div = min( abs(darr(1)*uarr(1)),abs(darr(2)*uarr(2)))
              py_div = min( abs(darr(3)*varr(3)),abs(darr(4)*varr(4)))
              pz_div = min( abs(darr(5)*warr(5)),abs(darr(6)*warr(6)))
              
              isConvergent = darr(2)*uarr(2) - darr(1)*uarr(1) < 0 
              if (isConvergent) then
                 uarr(1) = uarr(1) - px_div/darr(1)
                 uarr(2) = uarr(2) + px_div/darr(2)
              else ! comment out if you do not want to subtract outflows
                 uarr(1) = uarr(1) + px_div/darr(1)
                 uarr(2) = uarr(2) - px_div/darr(2)
              endif
              
              isConvergent = darr(4)*varr(4) - darr(3)*varr(3) < 0
              if (isConvergent) then
                 varr(3) = varr(3) - py_div/darr(3)
                 varr(4) = varr(4) + py_div/darr(4)
              else ! comment out if you do not want to subtract outflows
                 varr(3) = varr(3) + py_div/darr(3)
                 varr(4) = varr(4) - py_div/darr(4)
              endif
              
              isConvergent = darr(6)*warr(6) - darr(5)*warr(5) < 0
              if (isConvergent) then 
                 warr(5) = warr(5) - pz_div/darr(5)
                 warr(6) = warr(6) + pz_div/darr(6)
              else ! comment out if you do not want to subtract outflows
                 warr(5) = warr(5) + pz_div/darr(5)
                 warr(6) = warr(6) - pz_div/darr(6)
              endif
              
              ! subtract the rotational velocity field (x-y) plane
              ! ^y       <-        |4|        |-u|
              ! |       |  |     |1| |2|   |-v|  |+v|
              ! --->x    ->        |3|        |+u|
              Jz  = - varr(1)*darr(1) + varr(2)*darr(2) &
                   &   + uarr(3)*darr(3) - uarr(4)*darr(4)
              Jz  = Jz / 4.0
              
              varr(1) = varr(1) + Jz/darr(1) 
              varr(2) = varr(2) - Jz/darr(2) 
              uarr(3) = uarr(3) - Jz/darr(3)
              uarr(4) = uarr(4) + Jz/darr(4)
              
              ! subtract the rotational velocity field (y-z) plane
              ! ^z       <-        |6|        |-v|  
              ! |       |  |     |3| |4|   |-w|  |+w|
              ! --->y    ->        |5|        |+v|
              Jx  = - warr(3)*darr(3) + warr(4)*darr(4) &
                   &   + varr(5)*darr(5) - varr(6)*darr(6)
              Jx  = Jx / 4.0
              
              warr(3) = warr(3) + Jx/darr(3) 
              warr(4) = warr(4) - Jx/darr(4) 
              varr(5) = varr(5) - Jx/darr(5)
              varr(6) = varr(6) + Jx/darr(6)
              
              ! subtract the rotational velocity field (x-z) plane
              ! ^z       ->        |6|        |+u|  
              ! |       |  |     |1| |2|   |+w|  |-w|
              ! --->x    <-        |5|        |-u|
              Jy  = + warr(1)*darr(1) - warr(2)*darr(2) &
                   &   - uarr(5)*darr(5) + uarr(6)*darr(6)
              Jy  = Jy / 4.0
              
              warr(1) = warr(1) - Jy/darr(1) 
              warr(2) = warr(2) + Jy/darr(2) 
              uarr(5) = uarr(5) + Jy/darr(5)
              uarr(6) = uarr(6) - Jy/darr(6)
              
              ! From this point, uarr,uarr,warr is just the turbulent velocity
              trgv  = 0.0
              
              !x-direc
              ul    = (darr(2)*uarr(2) + d*uarr(0))/(darr(2)+d)
              ur    = (darr(1)*uarr(1) + d*uarr(0))/(darr(1)+d)
              trgv  = trgv + (ur-ul)**2
              !y-direc
              ul    = (darr(4)*varr(4) + d*varr(0))/(darr(4)+d)
              ur    = (darr(3)*varr(3) + d*varr(0))/(darr(3)+d)
              trgv  = trgv + (ur-ul)**2
              !z-direc
              ul    = (darr(6)*warr(6) + d*warr(0))/(darr(6)+d)
              ur    = (darr(5)*warr(5) + d*warr(0))/(darr(5)+d)
              trgv  = trgv + (ur-ul)**2
              !z-direc; tangential component - y
              ul    = (darr(6)*varr(6) + d*varr(0))/(darr(6)+d)
              ur    = (darr(5)*varr(5) + d*varr(0))/(darr(5)+d)
              trgv  = trgv + (ur-ul)**2
              !y-direc; tangential component - z
              ul    = (darr(4)*warr(4) + d*warr(0))/(darr(4)+d)
              ur    = (darr(3)*warr(3) + d*warr(0))/(darr(3)+d)
              trgv  = trgv + (ur-ul)**2
              !z-direc; tangential component - x
              ul    = (darr(6)*uarr(6) + d*uarr(0))/(darr(6)+d)
              ur    = (darr(5)*uarr(5) + d*uarr(0))/(darr(5)+d)
              trgv  = trgv + (ur-ul)**2
              !x-direc; tangential component - z
              ul    = (darr(2)*warr(2) + d*warr(0))/(darr(2)+d)
              ur    = (darr(1)*warr(1) + d*warr(0))/(darr(1)+d)
              trgv  = trgv + (ur-ul)**2
              !y-direc; tangential component - x
              ul    = (darr(4)*uarr(4) + d*uarr(0))/(darr(4)+d)
              ur    = (darr(3)*uarr(3) + d*uarr(0))/(darr(3)+d)
              trgv  = trgv + (ur-ul)**2
              !x-direc; tangential component - y
              ul    = (darr(2)*varr(2) + d*varr(0))/(darr(2)+d)
              ur    = (darr(1)*varr(1) + d*varr(0))/(darr(1)+d)
              trgv  = trgv + (ur-ul)**2
              ! now compute sound speed squared = cs^2 
              cs2  = uold(ind_cell(i),5)*(gamma -1.0)
              ! prevent numerical crash due to negative temperature
#if NENER>0
              !! Add radiation pressure to sound speed calculation
              !do irad = 0,nener-1
              !   cs2 = cs2 + uold(ind_cell(i),inener+irad) &
              !       / d * (gamma_rad(irad+1)-1.0)
              !end do
#endif
              cs2  = max(cs2,smallc**2)
              !rho_local = d
              rho_local = rho(ind_cell(i)) ! be careful with dark matter particles
              if(sf_lam.gt.0) then ! Jeans length criterion
                 ! Calculate "turbulent" Jeans length in cell units, lamjt 
                 ! (see e.g. Chandrasekhar 51, Bonazzola et al 87, Federrath & Klessen 2012 eq 36)
                 lamjt = (pi*trgv + sqrt(pi*pi*trgv*trgv + 36.0*pi*cs2*factG*rho_local*dx_loc**2))/(6.0*factG*rho_local*dx_loc**2)
                 if (lamjt > sf_lam) ok(i) = .false. ! Jeans length resolved: gas is stable 
              endif
              if(ok(i)) then ! Jeans length not resolved --> form stars to lower density and stabilise gas
                 ! corresponding virial parameter for homogeneous sphere <= 1.5 in turbulent dominated 
                 ! limit and <= 0.5 in pressure dominated limit (in good agreement with observations,
                 ! at least for massive (>= 10^5 M_sun) clouds see Kauffmann, Pillai, Goldsmith 2013, Fig.1)  
                 alpha0 = 5d0/(pi*factG*rho_local)*(trgv + cs2)/dx_loc**2
                 alpha_fk(i) =  alpha0
                 ! compute star formation efficiency per free-fall time (Federrath & Klessen 2012 eq 41) 
                 phi_t = 0.57 ; theta = 0.33 ! best fit values of the PN multi-ff model, originally from FK2012 and updated with later results  
                 sigs  = log(1.0+0.16*trgv/cs2) ! factor 0.16 is b^2 where b=0.4 for a mixture of turbulence forcing modes
                 scrit = log(0.067/theta**2*alpha0*trgv/cs2) ! best fit from Padoan & Nordlund MS model again 
                 sfr_ff(i) = eps_star_loc/2.0*phi_t*exp(3.0/8.0*sigs)*(2.0-erfc_pre_f08((sigs-scrit)/sqrt(2.0*sigs)))
              endif
              
           endif
        endif
     end do
     
  else
     if(myid.eq.1)then 
        write(*,*) ' star_maker does not match any of the pre-defined functions'
        write(*,*) star_maker 
     endif
     stop
  endif ! star maker

end subroutine get_sfeff_non_virial

end subroutine star_formation

#endif
!################################################################
!################################################################
!################################################################
!################################################################
subroutine getnbor(ind_cell,ind_father,ncell,ilevel)
  use amr_commons
  implicit none
  integer::ncell,ilevel
  integer,dimension(1:nvector)::ind_cell
  integer,dimension(1:nvector,0:twondim)::ind_father
  !-----------------------------------------------------------------
  ! This subroutine determines the 2*ndim neighboring cells
  ! cells of the input cell (ind_cell).
  ! If for some reasons they don't exist, the routine returns
  ! the input cell.
  !-----------------------------------------------------------------
  integer::i,j,iok,ind
  integer,dimension(1:nvector),save::ind_grid_father,pos
  integer,dimension(1:nvector,0:twondim),save::igridn,igridn_ok
  integer,dimension(1:nvector,1:twondim),save::icelln_ok


  if(ilevel==1)then
     write(*,*) 'Warning: attempting to form stars on level 1 --> this is not allowed ...'
     return
  endif

  ! Get father cell
  do i=1,ncell
     ind_father(i,0)=ind_cell(i)
  end do

  ! Get father cell position in the grid
  do i=1,ncell
     pos(i)=(ind_father(i,0)-ncoarse-1)/ngridmax+1
  end do

  ! Get father grid
  do i=1,ncell
     ind_grid_father(i)=ind_father(i,0)-ncoarse-(pos(i)-1)*ngridmax
  end do

  ! Get neighboring father grids
  call getnborgrids(ind_grid_father,igridn,ncell)

  ! Loop over position
  do ind=1,twotondim

     ! Select father cells that sit at position ind
     do j=0,twondim
        iok=0
        do i=1,ncell
           if(pos(i)==ind)then
              iok=iok+1
              igridn_ok(iok,j)=igridn(i,j)
           end if
        end do
     end do

     ! Get neighboring cells for selected cells
     if(iok>0)call getnborcells(igridn_ok,ind,icelln_ok,iok)

     ! Update neighboring father cells for selected cells
     do j=1,twondim
        iok=0
        do i=1,ncell
           if(pos(i)==ind)then
              iok=iok+1
              if(icelln_ok(iok,j)>0)then
                 ind_father(i,j)=icelln_ok(iok,j)
              else
                 ind_father(i,j)=ind_cell(i)
              end if
           end if
        end do
     end do

  end do


end subroutine getnbor
!##############################################################
!##############################################################
!##############################################################
!##############################################################
function erfc_pre_f08(x)

! complementary error function
  use amr_commons, ONLY: dp
  implicit none
  real(dp) erfc_pre_f08
  real(dp) x, y
  real(kind=8) pv, ph
  real(kind=8) q0, q1, q2, q3, q4, q5, q6, q7
  real(kind=8) p0, p1, p2, p3, p4, p5, p6, p7
  parameter(pv= 1.26974899965115684d+01, ph= 6.10399733098688199d+00)
  parameter(p0= 2.96316885199227378d-01, p1= 1.81581125134637070d-01)
  parameter(p2= 6.81866451424939493d-02, p3= 1.56907543161966709d-02)
  parameter(p4= 2.21290116681517573d-03, p5= 1.91395813098742864d-04)
  parameter(p6= 9.71013284010551623d-06, p7= 1.66642447174307753d-07)
  parameter(q0= 6.12158644495538758d-02, q1= 5.50942780056002085d-01)
  parameter(q2= 1.53039662058770397d+00, q3= 2.99957952311300634d+00)
  parameter(q4= 4.95867777128246701d+00, q5= 7.41471251099335407d+00)
  parameter(q6= 1.04765104356545238d+01, q7= 1.48455557345597957d+01)

  y = x*x
  y = exp(-y)*x*(p7/(y+q7)+p6/(y+q6) + p5/(y+q5)+p4/(y+q4)+p3/(y+q3) &
       &       + p2/(y+q2)+p1/(y+q1)+p0/(y+q0))
  if (x < ph) y = y+2.0d0/(exp(pv*x)+1.0d0)
  erfc_pre_f08 = y

  return

end function erfc_pre_f08
!##############################################################
!##############################################################
!##############################################################
subroutine pop3_analytic (m,proba)
   implicit none
   real(kind=8)::m,mcha,proba
   mcha = 100d0
   proba = m**(-1.3)*exp(-(mcha/m)**(1.6))
end subroutine pop3_analytic
!##############################################################
!##############################################################
!##############################################################
subroutine pop3_mass_random (mass)
   use amr_commons, ONLY:dp
   implicit none
   real(kind=dp)::mass
   real(kind=8)::fx_cap,m,logmin,logmax,fx,y
   logical::not_found
   integer,save::iseed=-399778
   real(kind=8),external::ran1
   m = 100. ! charateristic mass
   call pop3_analytic(m, fx_cap)
   fx_cap = fx_cap * 2

   logmin = 1.0
   logmax = 3.0

   not_found=.true.
   do while(not_found)
      m = ran1(iseed)
      m = m*(logmax-logmin)+logmin
      m = 10d0**m
      call pop3_analytic(m, fx)
      y = fx_cap*ran1(iseed)
      if (y.le.fx) not_found=.false.
   end do
   mass = m
end subroutine pop3_mass_random
!##############################################################
!##############################################################
!##############################################################
!tracer 
#ifdef MC_tracer
subroutine tracer2star(ind_tracer, proba, xstar, istar, nattach)
  use amr_commons
  use random
  use pm_commons
  implicit none

  integer, intent(in) :: nattach
  integer, dimension(1:nvector), intent(in) :: ind_tracer, istar
  real(dp), dimension(1:nvector), intent(in) :: proba
  real(dp), dimension(1:nvector, 1:3), intent(in) :: xstar

  logical, dimension(1:nvector), save :: attach = .false.
  integer :: i, idim, n
  real(dp) :: r

  n = 0
  do i = 1, nattach
     call ranf(tracer_seed, r)
     attach(i) = r < proba(i)
  end do

  ! Change particle pointer and kind
  do i = 1, nattach
     if (attach(i)) then
        n = n + 1
        partp(ind_tracer(i)) = istar(i)
        ! Tag particle with star id
        partp(ind_tracer(i)) = istar(i)
        typep(ind_tracer(i))%family = FAM_TRACER_STAR
     end if
  end do

  ! Change particle location
  do idim = 1, ndim
     do i = 1, nattach
        if (attach(i)) then
           xp(ind_tracer(i), idim) = xstar(i, idim)
        end if
     end do
  end do

  ! Save state of the particle
  do i = 1, nattach
     if (attach(i)) then
        ! Set move_flag to 1 to prevent further move
        move_flag(ind_tracer(i)) = 1
     end if
  end do

end subroutine tracer2star
#endif 
!tracer 
