module mod_usr
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use mod_mhd
  use mod_physics, only: phys_to_primitive,phys_to_conserved
  use mod_yinyang, only: yy_basis
  implicit none
  private
  public :: usr_init
  logical :: native=.false.,initial_refine=.false.
  integer :: last_check=-1
  integer,allocatable :: previous_owner(:,:,:,:)
contains
  subroutine usr_init()
    integer :: n
    namelist /amr_test/ native,initial_refine
    do n=1,size(par_files)
      open(unitpar,file=trim(par_files(n)),status='old')
      read(unitpar,amr_test,end=10)
10    close(unitpar)
    end do
    if(native) then
      call set_coordinate_system('spherical')
    else
      call set_coordinate_system('yinyang')
    end if
    call mhd_activate()
    usr_init_one_grid=>initial
    usr_refine_grid=>test_refinement_flags
    usr_process_adv_grid=>seed_step
    usr_after_refine=>check_after_refine
    usr_before_main_loop=>check_state
  end subroutine

  ! Patch-dependent scalars and Cartesian vectors that are bilinear in native
  ! donor angles. This allows a continuous-coordinate reference for the four
  ! donor samples without calling the interpolation code under test.
  pure subroutine profile(radius,theta,phi,patch,rho,p,psi,v,b)
    double precision,intent(in) :: radius,theta,phi
    integer,intent(in) :: patch
    double precision,intent(out) :: rho,p,psi,v(3),b(3)
    rho=1.d0+0.05d0*patch+0.1d0*radius+0.01d0*theta*phi
    p=1.d0+0.02d0*radius-0.02d0*theta+0.004d0*phi
    psi=0.001d0*(theta+phi+patch)
    v=[0.02d0*phi,0.01d0*theta,0.02d0+0.003d0*patch]
    b=[0.02d0+0.003d0*theta*phi,-0.01d0*radius,0.01d0*phi]
  end subroutine

  subroutine initial(ixI^L,ixO^L,w,x)
    integer,intent(in) :: ixI^L,ixO^L
    double precision,intent(in) :: x(ixI^S,ndim)
    double precision,intent(inout) :: w(ixI^S,nw)
    integer :: i,j,k
    double precision :: frame(3,3),r,p,psi,v(3),b(3)
    do k=ixOmin3,ixOmax3
      do j=ixOmin2,ixOmax2
        call yy_basis(x(ixOmin1,j,k,2),x(ixOmin1,j,k,3),yy_patch,frame)
        do i=ixOmin1,ixOmax1
          call profile(x(i,j,k,1),x(i,j,k,2),x(i,j,k,3),yy_patch,r,p,psi,v,b)
          w(i,j,k,rho_)=r
          w(i,j,k,p_)=p
          w(i,j,k,mom(:))=matmul(transpose(frame),v)
          w(i,j,k,mag(:))=matmul(transpose(frame),b)
          if(psi_>0) w(i,j,k,psi_)=psi
        end do
      end do
    end do
    call phys_to_conserved(ixI^L,ixO^L,w,x)
  end subroutine

  subroutine seed_step(igrid,level,ixI^L,ixO^L,qt,w,x)
    integer,intent(in) :: igrid,level,ixI^L,ixO^L
    double precision,intent(in) :: qt,x(ixI^S,ndim)
    double precision,intent(inout) :: w(ixI^S,nw)
    call initial(ixI^L,ixO^L,w,x)
  end subroutine

  subroutine test_refinement_flags(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    integer,intent(in) :: igrid,level,ixI^L,ixO^L
    double precision,intent(in) :: qt,w(ixI^S,nw),x(ixI^S,ndim)
    integer,intent(inout) :: refine,coarsen
    if(level==1) then
      ! Select one interior root outside the receiver and donor stencils.
      if(yy_patch==0.and.(it==1.or.(it==0.and.initial_refine))) then
        if(node(pig1_,igrid)==1.and.node(pig2_,igrid)==2.and.node(pig3_,igrid)==5) refine=1
      end if
    else if(it>=3) then
      coarsen=1
    end if
  end subroutine

  subroutine check_after_refine(n_coarsen,n_refine)
    integer,intent(in) :: n_coarsen,n_refine
    call check_state()
  end subroutine

  subroutine check_state()
    double precision :: w(ixG^T,nw),ref(nw),frame(3,3),ang(2),inplane(3),r,p,psi,v(3),b(3)
    double precision :: err,maxerr,dr,dtheta,dphi,radial,theta,phi
    integer :: n,g,i,j,k,gj,gk,gi,bi,bj,bk,a,c,received,expected,actual,fine,counts(3),count_all(3),moved
    integer,allocatable :: owner(:,:,:,:),donor(:,:,:,:),local_owner(:,:,:,:),local_donor(:,:,:,:)

    if(it==last_check) return
    last_check=it
    fine=0
    do n=1,igridstail
      if(node(plevel_,igrids(n))==2) fine=fine+1
    end do
    expected=0
    if(refine_max_level==2.and.yy_patch==0) then
      if(it==1.or.it==2.or.(it==0.and.initial_refine)) expected=8
    end if
    call MPI_ALLREDUCE(fine,actual,1,MPI_INTEGER,MPI_SUM,icomm,ierrmpi)
    if(actual/=expected) call mpistop('GLM AMR test: unexpected number of fine leaves')
    if(native) then
      if(mype==0) write(*,'(a,2(1x,i0))') 'NATIVE AMR',it,actual
      return
    end if

    allocate(owner(ng1(1),ng2(1),ng3(1),0:1),donor(ng1(1),ng2(1),ng3(1),0:1))
    allocate(local_owner(ng1(1),ng2(1),ng3(1),0:1),local_donor(ng1(1),ng2(1),ng3(1),0:1))
    local_owner=-1;local_donor=0;err=0.d0;received=0
    dr=(xprobmax1-xprobmin1)/domain_nx1
    dtheta=(xprobmax2-xprobmin2)/domain_nx2
    dphi=(xprobmax3-xprobmin3)/domain_nx3
    do n=1,igridstail
      g=igrids(n)
      if(node(plevel_,g)/=1) then
        if(node(pig2_,g)==1.or.node(pig2_,g)==ng2(2).or. &
            node(pig3_,g)==1.or.node(pig3_,g)==ng3(2)) &
          call mpistop('GLM AMR test: interior fixture created an angular receiver above level one')
        cycle
      end if
      local_owner(node(pig1_,g),node(pig2_,g),node(pig3_,g),yy_patch)=yy_world_rank
      if(node(pig2_,g)/=1.and.node(pig2_,g)/=ng2(1).and. &
          node(pig3_,g)/=1.and.node(pig3_,g)/=ng3(1)) cycle
      block=>ps(g)
      w=ps(g)%w
      call phys_to_primitive(ixG^LL,ixG^LL,w,ps(g)%x)
      do k=ixGlo3,ixGhi3
        gk=(node(pig3_,g)-1)*block_nx3+k-ixMlo3
        do j=ixGlo2,ixGhi2
          gj=(node(pig2_,g)-1)*block_nx2+j-ixMlo2
          if(gj>=0.and.gj<domain_nx2.and.gk>=0.and.gk<domain_nx3) cycle
          theta=ps(g)%x(ixMlo1,j,k,2);phi=ps(g)%x(ixMlo1,j,k,3)
          inplane=[-sin(theta)*cos(phi),cos(theta),sin(theta)*sin(phi)]
          ang=[acos(max(-1.d0,min(1.d0,inplane(3)))),atan2(inplane(2),inplane(1))]
          call yy_basis(theta,phi,yy_patch,frame)
          bj=floor((ang(1)-xprobmin2)/dtheta-0.5d0)
          bk=floor((ang(2)-xprobmin3)/dphi-0.5d0)
          do i=ixGlo1,ixGhi1
            gi=max(0,min(domain_nx1-1,(node(pig1_,g)-1)*block_nx1+i-ixMlo1))
            bi=1+gi/block_nx1
            radial=xprobmin1+(gi+0.5d0)*dr
            do c=0,1
              do a=0,1
                local_donor(bi,1+(bj+a)/block_nx2,1+(bk+c)/block_nx3,1-yy_patch)=1
              end do
            end do
            call profile(radial,ang(1),ang(2),1-yy_patch,r,p,psi,v,b)
            ref(rho_)=r;ref(p_)=p;ref(psi_)=psi
            ref(mom(:))=matmul(transpose(frame),v)
            ref(mag(:))=matmul(transpose(frame),b)
            if(any(.not.ieee_is_finite(w(i,j,k,:)))) call mpistop('GLM AMR test: non-finite ghost state')
            err=max(err,maxval(abs(w(i,j,k,:)-ref)))
            received=received+1
          end do
        end do
      end do
    end do
    call MPI_ALLREDUCE(local_owner,owner,size(owner),MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
    call MPI_ALLREDUCE(local_donor,donor,size(donor),MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
    if(any(donor==1.and.owner<0)) call mpistop('GLM AMR test: reference stencil has a missing level-one donor')
    call MPI_ALLREDUCE(err,maxerr,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
    counts=[received,merge(actual,0,yy_world_rank==0),merge(actual,0,yy_world_rank==yy_world_size/2)]
    call MPI_ALLREDUCE(counts,count_all,3,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierrmpi)
    if(count_all(1)==0.or.maxerr>2.d-13) call mpistop('GLM AMR test: wrong overset ghost values')
    moved=0
    if(allocated(previous_owner)) then
      moved=count(donor(:,:,:,0)==1.and.owner(:,:,:,0)/=previous_owner(:,:,:,0))
      if(any(donor(:,:,:,1)==1.and.owner(:,:,:,1)/=previous_owner(:,:,:,1))) &
        call mpistop('GLM AMR test: Yang donors unexpectedly migrated')
    end if
    previous_owner=owner
    if(yy_world_rank==0) write(*,'(a,5(1x,i0),1x,es14.6)') 'YY AMR',it,count_all,moved,maxerr
  end subroutine
end module
