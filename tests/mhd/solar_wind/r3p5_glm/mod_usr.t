
module mod_usr
  use mod_mhd
  use mod_physics, only: phys_gamma
  use mod_yinyang, only: yy_basis
  implicit none

  real(8), allocatable :: B_init(:,:,:,:)
  integer           :: status,readwrite,blocksize
  integer*4         :: i, j, k, ilevel, nx_ss, ny_ss
  double precision  :: k_B,miu0,mass_H,usr_grav,SRadius,rhob,Tiso,rhob1,qs

  double precision  :: rc, Vs
  double precision  :: Vout, V_surface
  integer           :: ix^D,nth,nph,nvec,nrr
  double precision :: dth,dph,dr0,mag_divisor
  integer             :: ixg2,ixg3
  double precision    :: xlen2, xlen3
  integer, dimension(4) :: shp

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods
    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr
    usr_gravity         => gravity
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output
    usr_refine_grid     => special_refine_grid

    call set_coordinate_system("yinyang")
    call mhd_activate()
  end subroutine usr_init

  subroutine initglobaldata_usr()
    use mod_global_parameters

    character(len=100)    :: filename
    real(8)               :: dr
    integer               :: ibc

    k_B = 1.3806d-16          ! erg*K^-1
    miu0 = 4.d0*dpi           ! Gauss^2 cm^2 dyne^-1
    mass_H = 1.67262d-24      ! g
    unit_length        = 6.955d10 ! cm
    unit_temperature   = 1.d6 ! K
    unit_numberdensity = 1.d9 ! cm^-3
    unit_density       = 1.4d0*mass_H*unit_numberdensity               ! 2.341668000000000E-015 g*cm^-3
    unit_pressure      = 2.3d0*unit_numberdensity*k_B*unit_temperature ! 0.317538000000000 erg*cm^-3
    unit_magneticfield = dsqrt(miu0*unit_pressure)                     ! 1.99757357615242 Gauss
    unit_velocity      = unit_magneticfield/dsqrt(miu0*unit_density)   ! 1.16448846777562E007 cm/s = 116.45 km/s
    unit_time          = unit_length/unit_velocity                     ! 5972.5794 s = 99.543 min

    usr_grav=-2.74d4*unit_length/unit_velocity**2  ! solar gravity
    SRadius=6.955d10/unit_length                   ! Solar radius

    rc = 3.45d0
    Vs = 117.54d0 

    rhob = 2.0d8/unit_numberdensity
    Tiso = 2.01d0
    mag_divisor = 1.0d0



    filename = './off/OFF_hmi_synchronic_0502_v140_Rss_3.5_lmax10.bin'

    call read_initial_magnetic_field(filename, shp, B_init)

    nrr = shp(2)
    nth = shp(3)
    nph = shp(4)
    dth = dpi/dble(nth-1)
    dph = (2.d0*dpi)/dble(nph-1)
    dr0 = dble(xprobmax1-xprobmin1)/dble(nrr-0)

    call cal_parker_solar_wind(1.0d0,rc,vs,V_surface)

    if (mype == 0) then
       call cal_parker_solar_wind(xprobmax1, rc, vs, vout)
       print *,'number density at maxmal r',(rhob*V_surface)/(Vout*xprobmax1**2) 
       print *,'velocity ar maxmal r', vout
    endif

  end subroutine initglobaldata_usr

  subroutine read_initial_magnetic_field(filename, shp, B_init)
      use mod_global_parameters
      character(len=*), intent(in)               :: filename
      double precision, allocatable, intent(out) :: B_init(:,:,:,:)
      integer, dimension(4), intent(out)         :: shp
      integer                                    :: ios
      integer                                    :: unit

      open(newunit=unit, file=filename, form='unformatted', access='stream', status='old')
      read(unit) shp

      if (any(shp /= [3,500,180,360])) call mpistop("Unexpected May 2 OFF dimensions")
      allocate(B_init(shp(1),shp(2),shp(3),shp(4)))
      read(unit) B_init
      close(unit)
      ! Latitude reversal is applied at lookup.

      if(mype==0) then
      print *,'MHD simulation initialized by amrvac PFSS model'
      print *, '  Initial Data Shape: ', shp
      print *, '  Load Initial Magnetic field data from: ', filename
      print *, 'First element before normalization: ', B_init(1,1,1,1)
      print *, 'Last  element before normalization: ', B_init(size(B_init,1),size(B_init,2),size(B_init,3),size(B_init,4))
      endif

      B_init = B_init/unit_magneticfield/mag_divisor
      if(mype==0) then
      print*,'    min, max values of initial B (normalized):', minval(B_init), maxval(B_init)
      endif
  end subroutine read_initial_magnetic_field

  function patch_field(irr,theta,phi) result(field)
    use mod_global_parameters, only: dpi, yy_patch
    integer, intent(in) :: irr
    double precision, intent(in) :: theta,phi
    double precision :: field(3),local(3,3),global(3,3),t,p
    integer :: it,ip
    t=theta
    p=phi
    if (yy_patch==1) then
      call yy_basis(theta,phi,1,local)
      t=atan2(hypot(local(1,1),local(2,1)),local(3,1))
      p=atan2(local(2,1),local(1,1))
    end if
    ! Match the author's endpoint spacing and ceiling convention.
    it=max(1,min(nth,ceiling(t/dth)))
    ip=max(1,min(nph,ceiling(modulo(p,2.d0*dpi)/dph)))
    field=B_init(:,irr,nth+1-it,ip)
    if (yy_patch==1) then
      call yy_basis(t,p,0,global)
      field=matmul(transpose(local),matmul(global,field))
    end if
  end function patch_field

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)

    use mod_global_parameters    
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: qs
    integer          :: irr,ith,iph
    logical, save:: first=.true.

    if (first .and. mype==0) then
      print *,'Relax a solor wind model with politropic MHD'
      print *,'User Tiso(MK): ', Tiso
      print *,'V_surface(unit_velocity): ', V_surface
      first=.false.
    end if

    qs=1.0d0
    w(ixO^S,mom(1:3))=0.0d0

    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      call cal_parker_solar_wind(x(ix^D,1), rc, vs, vout)
      w(ix^D,rho_)  = (rhob*V_surface)/(Vout*x(ix^D,1)**2)
      w(ix^D,mom(1))= Vout*w(ix^D,rho_)

      irr = ceiling((x(ix1,ix2,ix3,1)-xprobmin1)/dr0)

      w(ix^D,mag(1:3))=patch_field(irr,x(ix^D,2),x(ix^D,3))
    {end do\}

    w(ixO^S,p_)   = w(ixO^S,rho_)*Tiso/(phys_gamma-1.0d0)

    if(mhd_glm) w(ixO^S,psi_)=0.d0

  end subroutine initonegrid_usr

  subroutine cal_parker_solar_wind(ra_pksw, rc_pksw, vs_pksw, v_pksw)
  use mod_global_parameters
  real*8, intent(in)  :: ra_pksw, rc_pksw, vs_pksw
  real*8, intent(out) :: v_pksw
  real*8              :: W0_pksw, Wn1_pksw, Dr_pksw, nDr_pksw, plw, L_pksw, M_pksw

   Dr_pksw = ((ra_pksw/rc_pksw)**(-4.0d0))*dexp(4.0d0*(1.0d0-(rc_pksw/ra_pksw))-1.0d0)
   nDr_pksw = -1.0d0*Dr_pksw
   if (ra_pksw .le. rc_pksw-0.05d0) then
       if (nDr_pksw .gt. -0.25) then
            W0_pksw = nDr_pksw-nDr_pksw**2.0d0 + 1.5d0*nDr_pksw*nDr_pksw*nDr_pksw
       endif
       if (nDr_pksw .le. -0.25) then
           plw=sqrt(2.0d0*(exp(1.0d0)*nDr_pksw+1.0d0))
           W0_pksw=-1.0d0+plw-(1.0d0/3.0d0)*plw**2.0d0+(11.0d0/72.0d0)*plw*plw*plw
       endif
     v_pksw = dsqrt(-1.0d0*vs_pksw*vs_pksw*W0_pksw)
   endif

    if (ra_pksw .ge. rc_pksw-0.05d0 .and. ra_pksw .le. rc_pksw+0.05d0) then
     v_pksw = vs_pksw*sqrt(3.0d0-2.0d0*rc_pksw/ra_pksw)
    endif

   if (ra_pksw .gt. rc_pksw+0.05d0) then
      L_pksw=log(abs(nDr_pksw))
      M_pksw=log(abs(L_pksw))
      Wn1_pksw = L_pksw - M_pksw + (M_pksw/L_pksw) + &
                 M_pksw*(M_pksw-2.0d0)/(2.0d0*L_pksw*L_pksw)+&
                 M_pksw*(6.0d0-9.0d0*M_pksw+2.0d0*M_pksw*M_pksw)/(6.0d0*L_pksw*L_pksw*L_pksw)
      if (abs(nDr_pksw) .ge. 0.28796090d0) then
         plw=-sqrt(2.0d0*(exp(1.0d0)*nDr_pksw+1.0d0))
         Wn1_pksw=-1.0d0+plw-(1.0d0/3.0d0)*plw*plw+(11.0d0/72.0d0)*plw*plw*plw
      endif
     v_pksw = dsqrt(-1.0d0*vs_pksw*vs_pksw*Wn1_pksw)
   endif

     v_pksw = v_pksw*1.0d5/unit_velocity  ! Km/s -> normalized velocity

  end subroutine

  subroutine specialbound_usr(qdt,qt,ixI^L,ixO^L,iB,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qdt,qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: Qp(ixI^S),xlen^D
    double precision, allocatable :: magCC(:,:,:),xex(:,:,:)
    integer :: ix^D,ixOs^L,jxO^L,idir,ixbc^D,ixA^L

    integer :: ixIM^L, ith, iph
    double precision :: tmp1(ixI^S),tmp2(ixI^S),pth(ixI^S)
    double precision :: coeffrho
    double precision :: q,q1,a0,a1,a2,a3,b0,b1,b2,c1,c2,c3,c4

    q  = 1.0d0
    q1 = 1.0d0/q
    a0 = -1.0d0-1.0d0/(1.0d0+q)-1/(1.0d0+q+q**2)
    a1 = 1.0d0+1.0d0/q**2+1.0d0/q**2
    a2 = -(1.0d0+q+q**2)/(q**3+q**4)
    a3 = 1.0d0/(q**3+q**4+q**5)
    b0 = -(2.0d0+q1)/(1.0d0+q1)
    b1 = 1.0d0+1.0d0/q1
    b2 = -1.0d0/(1.0d0+q1**2)

    c1 = -24.0d0*q1**6/(1.0d0-2.0d0*q1-3.0d0*q1**2+2.0d0*q1**3+8.0d0*q1**4+12.0d0*q1**5-24.0d0*q1**6)
    c2 = 6.0d0*q1**3/(2.0d0-7.0d0*q1+2.0d0*q1**2+14.0d0*q1**3-12.0d0*q1**4)
    c3 = 8.0d0*q1/(-6.0d0+17.0d0*q1+6.0d0*q1**2-51.0d0*q1**3+36.0d0*q1**4)
    c4 = 3.0d0/(3.0d0-4.0d0*q1-6.0d0*q1**2-4.0d0*q1**3+16.0d0*q1**4+24.0d0*q1**5-32.0d0*q1**6)

    select case(iB)
     case(1)

       w(ixO^S,mom(:))=0.0d0

       {do ix^DB=ixOmin^DB,ixOmax^DB\}
         call cal_parker_solar_wind(x(ix^D,1), rc, vs, vout)
         w(ix^D,rho_)  = (rhob*V_surface)/(Vout*x(ix^D,1)**2)
         w(ix^D,mom(1))= vout*w(ix^D,rho_)
         w(ix^D,p_)    = w(ix^D,rho_)*Tiso/(phys_gamma-1.0d0)
       {end do\}
       w(ixO^S,mom(2)) = zero
       w(ixO^S,mom(3)) = zero

       do ix2=ixOmin2,ixOmax2
         do ix3 = ixOmin3,ixOmax3

           w(ixOmax1,ix2,ix3,mag(1:3))=patch_field(1,x(ixOmax1,ix2,ix3,2),x(ixOmax1,ix2,ix3,3))
         end do
       end do

       do ix1=ixOmax1-1,ixOmin1,-1
         w(ix1^%1ixO^S,mag(1:3))=(a0-a1)/a0*w(ix1+1^%1ixO^S,mag(1:3))+&
                                 (a1-a2)/a0*w(ix1+2^%1ixO^S,mag(1:3))+&
                                 (a2-a3)/a0*w(ix1+3^%1ixO^S,mag(1:3))+&
                                 (a3- 0)/a0*w(ix1+4^%1ixO^S,mag(1:3))
       end do

     case(2)

       do ix1 = ixOmin1,ixOmax1

         w(ix1^%1ixO^S,rho_) = (-b1/b0)*w(ix1-1^%1ixO^S,rho_)+&
                               (-b2/b0)*w(ix1-2^%1ixO^S,rho_)
         w(ix1^%1ixO^S,p_)   = (-b1/b0)*w(ix1-1^%1ixO^S,p_)+&
                               (-b2/b0)*w(ix1-2^%1ixO^S,p_)

         w(ix1^%1ixO^S,mom(2:3)) = zero
         w(ix1^%1ixO^S,mag(2:3)) = zero

         w(ix1^%1ixO^S,mag(1)) = ((-b1/b0)*w(ix1-1^%1ixO^S,mag(1))*x(ix1-1^%1ixO^S,1)**2+&
                                  (-b2/b0)*w(ix1-2^%1ixO^S,mag(1))*x(ix1-2^%1ixO^S,1)**2)/x(ix1^%1ixO^S,1)**2 

       end do

       do ix1 = ixOmin1,ixOmax1
         w(ix1^%1ixO^S,mom(1)) = w(ixOmin1-1^%1ixO^S,mom(1))
         where(w(ix1^%1ixO^S,mom(1))<0.d0)
           w(ix1^%1ixO^S,mom(1))=0.d0
         end where
       end do

     case(3)

     case(4)

     case(5)

     case(6)

     case default
       call mpistop("Special boundary is not defined for this region")
    end select

  end subroutine specialbound_usr

  subroutine special_refine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    refine = -1
    coarsen = -1

  end subroutine special_refine_grid

  subroutine getggrav(ggrid,ixI^L,ixO^L,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: ggrid(ixI^S)

    ggrid(ixO^S)=usr_grav*SRadius**2/(x(ixO^S,1))**2
  end subroutine

  subroutine gravity(ixI^L,ixO^L,wCT,x,gravity_field)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(out)   :: gravity_field(ixI^S,ndim)
    double precision                :: ggrid(ixI^S)

    gravity_field=0.d0
    call getggrav(ggrid,ixI^L,ixO^L,x)
    gravity_field(ixO^S,1)=ggrid(ixO^S)
  end subroutine gravity

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision                   :: w1(ixI^S,nw+nwauxio)
    double precision                   :: qvec(ixI^S,1:ndim)
    double precision                   :: tmp(ixI^S)
    double precision                   :: current(ixI^S,1:ndim), Lorentz(ixI^S,1:ndim)
    double precision                   :: xlen^D
    integer :: ix^D,ixbc^D,idirmin,idir,jdir,kdir
    double precision                   :: divb(ixI^S),B2(ixI^S),dip(ixI^S)
    double precision                   :: divv(ixI^S), div_free_fi(ixI^S)
    double precision                   :: vels(ixI^S,1:ndim)

    w(ixO^S,nw+1:nw+nwauxio)=zero
    w(ixO^S,nw+1 )=w(ixO^S,mag(1))*unit_magneticfield
    w(ixO^S,nw+2 )=w(ixO^S,mag(2))*unit_magneticfield
    w(ixO^S,nw+3 )=w(ixO^S,mag(3))*unit_magneticfield

    call get_current(w,ixI^L,ixO^L,idirmin,current)
    B2(ixO^S)=sum(w(ixO^S,mag(1:ndir))**2,dim=ndim+1)
    w(ixO^S,nw+4 )=current(ixO^S,1)
    w(ixO^S,nw+5 )=current(ixO^S,2)
    w(ixO^S,nw+6 )=current(ixO^S,3)

    call get_divb(w,ixI^L,ixO^L,tmp)
    w(ixO^S,nw+7 )=tmp(ixO^S)

    call mhd_get_pthermal(w,x,ixI^L,ixO^L,tmp)
    w(ixO^S,nw+8 )=tmp(ixO^S)

    w(ixO^S,nw+9 )=w(ixO^S,mom(1))/w(ixO^S,rho_)*unit_velocity/1.d5
    vels(ixI^S, 1)=w(ixI^S,mom(1))/w(ixI^S,rho_)
    vels(ixI^S, 2)=w(ixI^S,mom(2))/w(ixI^S,rho_)
    vels(ixI^S, 3)=w(ixI^S,mom(3))/w(ixI^S,rho_)

    call divvector(vels(ixI^S,1:ndir),ixI^L,ixO^L,divv)
    w(ixO^S,nw+10)=divv(ixO^S)
    call get_normalized_divb(w,ixI^L,ixO^L,div_free_fi)
    w(ixO^S,nw+11)=div_free_fi
    w(ixO^S,nw+12)=current(ixO^S,2)*w(ixO^S,mag(3))-current(ixO^S,3)*w(ixO^S,mag(2))

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)

    character(len=*) :: varnames
    varnames='Br Bt Bp Jr Jt Jp divB pth Vr divv fi Lr'

  end subroutine specialvarnames_output

end module mod_usr
