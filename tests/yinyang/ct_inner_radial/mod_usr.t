module mod_usr
  use mod_mhd
  use mod_physics, only: phys_get_pthermal
  use mod_functions_bfield, only: get_divb
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  logical  ::  native=.false.
  double precision  ::  drive=0.d0
  double precision,parameter  ::  temp=2.01d0,rc=3.49584d0,rhob=0.05d0,bref=0.02d0
contains
  subroutine usr_init()
    integer  ::  n
    namelist /radial_test/ native,drive
    do n=1,size(par_files)
      open(unitpar,file=trim(par_files(n)),status='old')
      read(unitpar,radial_test,end=10)
10    close(unitpar)
    end do
    if(native)then
      call set_coordinate_system('spherical')
    else
      call set_coordinate_system('yinyang')
    end if
    usr_init_one_grid=>initial
    usr_special_bc=>inner
    usr_set_electric_field=>boundary_emf
    usr_gravity=>gravity
    usr_process_global=>check
    usr_process_adv_global=>advanced
    call mhd_activate()
  end subroutine

  ! Subsonic isothermal Parker branch on this shell, in code units.
  pure double precision function speed(r)
    double precision,intent(in) :: r
    double precision :: rhs,lo,hi,y
    integer :: n
    rhs=4*(log(r/rc)+rc/r)-3
    lo=-rhs;hi=0.d0
    do n=1,54
      y=(lo+hi)/2
      if(exp(y)-y>rhs)then
        lo=y
      else
        hi=y
      end if
    end do
    speed=sqrt(temp)*exp((lo+hi)/4)
  end function

  subroutine plasma(ixI^L,ixO^L,w,x)
    integer,intent(in) :: ixI^L,ixO^L
    double precision,intent(in) :: x(ixI^S,ndim)
    double precision,intent(inout) :: w(ixI^S,nw)
    integer :: i
    double precision :: r,v,rho
    do i=ixOmin1,ixOmax1
      r=x(i,ixOmin2,ixOmin3,1);v=speed(r)
      rho=rhob*speed(1.d0)/(r*r*v)
      w(i,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)=rho
      w(i,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(1))=rho*v
      w(i,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(2:3))=0.d0
      w(i,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)=rho*temp*eos%inv_gamma_minus_1
    end do
  end subroutine

  subroutine initial(ixI^L,ixO^L,w,x)
    integer,intent(in) :: ixI^L,ixO^L
    double precision,intent(in) :: x(ixI^S,ndim)
    double precision,intent(inout) :: w(ixI^S,nw)
    integer :: i
    double precision :: r
    if(.not.mhd_internal_e)call mpistop('radial test uses internal energy')
    block%ws=0.d0
    do i=lbound(block%ws,1),ubound(block%ws,1)
      r=x(ixImin1,ixImin2,ixImin3,1)+(i-ixImin1+0.5d0)*dxlevel(1)
      block%ws(i,:,:,1)=bref/r**2
    end do
    call mhd_face_to_center(ixO^L,block)
    call plasma(ixI^L,ixO^L,w,x)
  end subroutine

  subroutine inner(qdt,qt,ixI^L,ixO^L,iB,w,x)
    integer,intent(in) :: ixI^L,ixO^L,iB
    double precision,intent(in) :: qdt,qt,x(ixI^S,ndim)
    double precision,intent(inout) :: w(ixI^S,nw)
    if(iB/=1)call mpistop('unexpected special boundary')
    call plasma(ixI^L,ixO^L,w,x)
  end subroutine

  subroutine gravity(ixI^L,ixO^L,w,x,grav)
    integer,intent(in) :: ixI^L,ixO^L
    double precision,intent(in) :: w(ixI^S,nw),x(ixI^S,ndim)
    double precision,intent(out) :: grav(ixI^S,ndim)
    grav=0.d0
    grav(ixO^S,1)=-2*temp*rc/x(ixO^S,1)**2
  end subroutine

  ! Integrate E=drive*(z_hat cross r_hat) over the actual inner-face edges at r=1.
  subroutine boundary_emf(ixI^L,ixO^L,qt,qdt,fE,s)
    integer,intent(in) :: ixI^L,ixO^L
    double precision,intent(in) :: qt,qdt
    type(state) :: s
    double precision,intent(inout) :: fE(ixI^S,sdim:3)
    integer :: i,j,k
    double precision :: t,p,h,dp
    if(.not.s%is_physical_boundary(1))return
    i=ixOmin1-1;h=dxlevel(2);dp=dxlevel(3)
    do k=ixOmin3-1,ixOmax3
      do j=ixOmin2-1,ixOmax2
        t=s%x(i,j,k,2)+h/2;p=s%x(i,j,k,3)
        if(.not.native.and.yy_patch==1)then
          fE(i,j,k,2)=qdt*drive*cos(p+dp/2)*h
          fE(i,j,k,3)=-qdt*drive*sin(t)*cos(t)*sin(p)*2*sin(dp/2)
        else
          fE(i,j,k,2)=0.d0
          fE(i,j,k,3)=qdt*drive*sin(t)**2*dp
        end if
      end do
    end do
  end subroutine

  subroutine advanced(it,qt)
    integer,intent(in) :: it
    double precision,intent(in) :: qt
    call check(it+1,qt+dt)
  end subroutine

  subroutine check(it,qt)
    integer,intent(in) :: it
    double precision,intent(in) :: qt
    integer :: n,g,i,j,k,rank,counts(2),total(2)
    double precision :: pmin,pworld,pressure(ixG^T),divb(ixG^T),err(4),worst(4),t,p,h,dp,zint,area,ref,r,rho
    err=0;counts=0;pmin=huge(1.d0)
    do n=1,igridstail
      g=igrids(n)
      if(.not.ps(g)%is_physical_boundary(1))cycle
      block=>ps(g)
      call get_divb(block%w,ixG^LL,ixM^LL,divb)
      call phys_get_pthermal(block%w,block%x,ixG^LL,ixM^LL,pressure)
      i=ixMlo1-1
      do k=ixMlo3,ixMhi3
        do j=ixMlo2,ixMhi2
          t=block%x(i,j,k,2);p=block%x(i,j,k,3)
          h=block%dx(i,j,k,2);dp=block%dx(i,j,k,3)
          area=block%surfaceC(i,j,k,1)
          zint=0.5d0*(sin(t+h/2)**2-sin(t-h/2)**2)*dp
          if(.not.native.and.yy_patch==1) zint=sin(p)*sin(dp/2)*(h-cos(2*t)*sin(h))
          ref=bref*area-2*drive*qt*zint
          if(.not.ieee_is_finite(block%ws(i,j,k,1)+divb(i+1,j,k)))call mpistop('nonfinite CT state')
          if(.not.ieee_is_finite(pressure(i+1,j,k)))call mpistop('nonfinite inner pressure')
          pmin=min(pmin,pressure(i+1,j,k))
          err(1)=max(err(1),abs(block%ws(i,j,k,1)*area-ref))
          err(2)=max(err(2),abs(divb(i+1,j,k)))
          counts(1)=counts(1)+1
        end do
      end do
      ! Include all radial-angular ghost corners of the overset patches.
      do i=ixGlo1,ixMlo1-1
        r=block%x(i,ixMlo2,ixMlo3,1)
        rho=rhob*speed(1.d0)/(r*r*speed(r))
        do k=ixGlo3,ixGhi3
          do j=ixGlo2,ixGhi2
            if(native.and.(j<ixMlo2.or.j>ixMhi2.or.k<ixMlo3.or.k>ixMhi3))cycle
            err(3)=max(err(3),abs(block%w(i,j,k,rho_)-rho)/rho)
            err(4)=max(err(4),abs(block%w(i,j,k,e_)/(rho*temp*eos%inv_gamma_minus_1)-1))
            counts(2)=counts(2)+1
          end do
        end do
      end do
    end do
    call MPI_ALLREDUCE(err,worst,4,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
    call MPI_ALLREDUCE(counts,total,2,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierrmpi)
    call MPI_ALLREDUCE(pmin,pworld,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierrmpi)
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierrmpi)
    if(rank==0)write(*,'(a,i5,es16.7,2i9,5es17.8)')'RADIAL_CT ',it,qt,total,worst,pworld
    if(pworld<=0.or.any(total<=0)) call mpistop('invalid pressure or empty sample')
    if(worst(1)>1.d-12) call mpistop('incorrect radial face flux')
    if(maxval(worst(2:4))>1.d-10) call mpistop('invalid divergence or fluid ghosts')
  end subroutine
end module
