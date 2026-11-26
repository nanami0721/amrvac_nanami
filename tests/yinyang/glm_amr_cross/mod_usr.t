module mod_usr
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use mod_mhd
  use mod_eos, only: eos
  use mod_mhd_phys, only: mhd_internal_e
  use mod_physics, only: phys_to_primitive,phys_to_conserved
  use mod_ghostcells_update, only: getbc
  use mod_yinyang, only: yy_basis
  implicit none
  private
  public :: usr_init
  logical :: initial_refine=.false., uniform_state=.false., special_inner=.false.
  integer :: last_check=-1, local_special_corners=0, previous_level
  integer,allocatable :: previous_cover(:,:,:,:)

contains

  subroutine usr_init()
    integer :: n
    namelist /cross_test/ initial_refine,uniform_state,special_inner
    do n=1,size(par_files)
      open(unitpar,file=trim(par_files(n)),status='old')
      read(unitpar,cross_test,end=10)
10    close(unitpar)
    end do
    call set_coordinate_system('yinyang')
    call mhd_activate()
    usr_init_one_grid=>cross_initial
    usr_refine_grid=>cross_flags
    usr_after_refine=>cross_after_refine
    usr_before_main_loop=>cross_check
    if(special_inner) usr_special_bc=>cross_special_bc
  end subroutine

  subroutine profile(r,x2,x3,patch,w)
    double precision,intent(in) :: r,x2,x3
    integer,intent(in) :: patch
    double precision,intent(out) :: w(nw)
    w(rho_)=1.d0+0.01d0*patch
    w(p_)=1.d0
    w(mom(:))=[0.012d0,-0.014d0,0.015d0]
    w(mag(:))=[0.011d0,0.012d0,-0.013d0]
    w(psi_)=0.001d0
    if(.not.uniform_state) then
      w(rho_)=w(rho_)+0.03d0*r+0.01d0*x2-0.01d0*x3
      w(p_)=w(p_)+0.02d0*r-0.01d0*x2+0.02d0*x3
      w(mom(:))=w(mom(:))+0.005d0*[x2,x3,r]
      w(mag(:))=w(mag(:))+0.003d0*[x3,r,x2]
      w(psi_)=w(psi_)+0.001d0*x2
    end if
    if(special_inner) w(mag(:))=0.d0
  end subroutine

  subroutine cross_initial(ixI^L,ixO^L,w,x)
    integer,intent(in) :: ixI^L,ixO^L
    double precision,intent(in) :: x(ixI^S,ndim)
    double precision,intent(inout) :: w(ixI^S,nw)
    integer :: i,j,k
    double precision :: frame(3,3),v(nw)
    do k=ixOmin3,ixOmax3
      do j=ixOmin2,ixOmax2
        call yy_basis(x(ixOmin1,j,k,2),x(ixOmin1,j,k,3),yy_patch,frame)
        do i=ixOmin1,ixOmax1
          call profile(x(i,j,k,1),x(i,j,k,2),x(i,j,k,3),yy_patch,v)
          v(mom(:))=matmul(transpose(frame),v(mom(:)))
          v(mag(:))=matmul(transpose(frame),v(mag(:)))
          w(i,j,k,:)=v
        end do
      end do
    end do
    call phys_to_conserved(ixI^L,ixO^L,w,x)
  end subroutine

  subroutine special_state(theta,phi,patch,w)
    double precision,intent(in) :: theta,phi
    integer,intent(in) :: patch
    double precision,intent(out) :: w(nw)
    double precision :: frame(3,3)
    call profile(1.d0,theta,phi,patch,w)
    call yy_basis(theta,phi,patch,frame)
    w(mag(:))=0.02d0*frame(:,1)
    w(rho_)=w(rho_)+0.1d0*(1.d0+patch)
  end subroutine

  subroutine cross_special_bc(qdt,qt,ixI^L,ixO^L,iB,w,x)
    integer,intent(in) :: ixI^L,ixO^L,iB
    double precision,intent(in) :: qdt,qt,x(ixI^S,ndim)
    double precision,intent(inout) :: w(ixI^S,nw)
    integer :: i,j,k
    double precision :: frame(3,3),v(nw)
    if(iB/=1) call mpistop('cross-test: unexpected special boundary')
    if(ixOmin2==ixGlo2.and.ixOmax2==ixGhi2.and. &
      ixOmin3==ixGlo3.and.ixOmax3==ixGhi3) local_special_corners=local_special_corners+1
    do k=ixOmin3,ixOmax3
      do j=ixOmin2,ixOmax2
        call yy_basis(x(ixOmin1,j,k,2),x(ixOmin1,j,k,3),yy_patch,frame)
        do i=ixOmin1,ixOmax1
          call special_state(x(i,j,k,2),x(i,j,k,3),yy_patch,v)
          v(mom(:))=matmul(transpose(frame),v(mom(:)))
          v(mag(:))=matmul(transpose(frame),v(mag(:)))
          w(i,j,k,:)=v
        end do
      end do
    end do
    call phys_to_conserved(ixI^L,ixO^L,w,x)
  end subroutine

  subroutine cross_flags(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    integer,intent(in) :: igrid,level,ixI^L,ixO^L
    double precision,intent(in) :: qt,w(ixI^S,nw),x(ixI^S,ndim)
    integer,intent(inout) :: refine,coarsen
    double precision,parameter :: phi_anchor=0.015d0
    if(it>=4) then
      if(level>1) coarsen=1
      return
    end if
    if(it==0.and..not.initial_refine) return
    if(yy_patch/=0.or.level>=refine_max_level) return
    if(node(pig1_,igrid)/=1.or.node(pig2_,igrid)/=1) return
    if(rnode(rpxmin3_,igrid)<=phi_anchor.and.rnode(rpxmax3_,igrid)>phi_anchor) refine=1
  end subroutine

  subroutine cross_after_refine(n_coarsen,n_refine)
    integer,intent(in) :: n_coarsen,n_refine
    ! Startup can pass through several different layouts at the same time.
    last_check=-1
    call cross_check()
  end subroutine

  subroutine cross_check()
    integer :: local_high,high,n,g,factor,lo(3),hi(3),width(3),bn(3)
    integer,allocatable :: map(:,:,:,:),owned(:,:,:,:)
    double precision,allocatable :: edges(:,:)
    integer :: l,nr,e,edge_section
    double precision :: dr,error,max_error
    integer :: counts(8),total(8)

    if(it==last_check) return
    last_check=it
    local_high=levmax
    call MPI_ALLREDUCE(local_high,high,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
    bn=[block_nx1,block_nx2,block_nx3]
    width=[ng1(high),ng2(high),ng3(high)]
    allocate(owned(width(1),width(2),width(3),0:1),map(width(1),width(2),width(3),0:1))
    owned=0
    do n=1,igridstail
      g=igrids(n)
      block=>ps(g)
      call cross_initial(ixG^LL,ixM^LL,ps(g)%w,ps(g)%x)
      l=node(plevel_,g)
      factor=2**(high-l)
      lo=([node(pig1_,g),node(pig2_,g),node(pig3_,g)]-1)*factor+1
      hi=lo+factor-1
      owned(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),yy_patch)=l*(yy_world_size+1)+yy_world_rank+1
    end do
    call MPI_ALLREDUCE(owned,map,size(map),MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
    if(minval(map)==0) call mpistop('cross-test: missing leaf in reference map')
    allocate(edges(0:width(1)*bn(1),high))
    do l=1,high
      nr=ng1(l)*bn(1)
      edges(0,l)=xprobmin1
      if(stretch_type(1)==stretch_symm) edge_section=nstretchedblocks(l,1)*bn(1)/2
      do e=1,nr
        select case(stretch_type(1))
        case(stretch_uni)
          dr=dxfirst(l,1)*qstretch(l,1)**(e-1)
        case(stretch_symm)
          if(e<=edge_section) then
            dr=dxfirst(l,1)*qstretch(l,1)**(edge_section-e)
          else if(e>nr-edge_section) then
            dr=dxfirst(l,1)*qstretch(l,1)**(e-nr+edge_section-1)
          else
            dr=dxmid(l,1)
          end if
        case default
          dr=dx(1,l)
        end select
        edges(e,l)=edges(e-1,l)+dr
      end do
    end do
    local_special_corners=0
    call getbc(global_time,0.d0,ps,iwstart,nwgc)
    error=0.d0
    counts=0
    do n=1,igridstail
      g=igrids(n)
      call compare_block(g,high,map,edges,error,counts)
    end do
    previous_cover=map
    previous_level=high
    counts(6)=local_special_corners
    call MPI_ALLREDUCE(counts,total,size(total),MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierrmpi)
    call MPI_ALLREDUCE(error,max_error,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
    if(yy_world_rank==0) write(*,'(a,10(1x,i0),1x,es13.5)') 'CROSS AMR',it,high,total,max_error
    if(total(1)==0.or.max_error>2.d-12) call mpistop('cross-test: ghost comparison failed')
  end subroutine

  subroutine compare_block(g,high,map,edges,error,counts)
    integer,intent(in) :: g,high,map(:,:,:,0:)
    double precision,intent(in) :: edges(0:,:)
    double precision,intent(inout) :: error
    integer,intent(inout) :: counts(:)
    double precision :: actual(ixG^T,nw),expected(nw),value(nw),theta,phi,dq(2),other(3),ang(2)
    double precision :: frac(2),frame(3,3)
    integer :: i,j,k,a,b,l,root(3),bn(3),nn(3),v(3),idx(3),base(2),factor,donors(0:1,0:1),ri
    integer :: cell(3),current,previous
    logical :: at_wall,is_mixed
    l=node(plevel_,g)
    bn=[block_nx1,block_nx2,block_nx3]
    nn=[ng1(l),ng2(l),ng3(l)]*bn
    root=([node(pig1_,g),node(pig2_,g),node(pig3_,g)]-1)*bn
    if(root(2)>0.and.root(2)+bn(2)<nn(2).and.root(3)>0.and.root(3)+bn(3)<nn(3)) return
    dq=[xprobmax2-xprobmin2,xprobmax3-xprobmin3]/nn(2:3)
    factor=2**(high-l)
    block=>ps(g)
    actual=ps(g)%w
    call phys_to_primitive(ixG^LL,ixG^LL,actual,ps(g)%x)
    do k=ixGlo3,ixGhi3
      do j=ixGlo2,ixGhi2
        v(2:3)=root(2:3)+[j-ixMlo2,k-ixMlo3]
        if(all(v(2:3)>=0).and.all(v(2:3)<nn(2:3))) cycle
        theta=ps(g)%x(ixMlo1,j,k,2)
        phi=ps(g)%x(ixMlo1,j,k,3)
        other=[-sin(theta)*cos(phi),cos(theta),sin(theta)*sin(phi)]
        ang=[acos(other(3)),atan2(other(2),other(1))]
        frac=(ang-[xprobmin2,xprobmin3])/dq-0.5d0
        base=floor(frac)
        frac=frac-base
        call yy_basis(theta,phi,yy_patch,frame)
        do i=ixGlo1,ixGhi1
          v(1)=root(1)+i-ixMlo1
          at_wall=v(1)<0.or.v(1)>=nn(1)
          do b=0,1
            do a=0,1
              cell=[max(0,min(nn(1)-1,v(1))),base(1)+a,base(2)+b]
              idx=(cell*factor)/bn+1
              current=map(idx(1),idx(2),idx(3),1-yy_patch)
              donors(a,b)=current/(yy_world_size+1)
              if(allocated(previous_cover).and..not.at_wall) then
                if(previous_level>=l) then
                  idx=cell*2**(previous_level-l)/bn+1
                else
                  idx=cell/(bn*2**(l-previous_level))+1
                end if
                previous=previous_cover(idx(1),idx(2),idx(3),1-yy_patch)
                if(previous/(yy_world_size+1)==donors(a,b).and.previous/=current) counts(8)=counts(8)+1
              end if
            end do
          end do
          if(any(abs(donors-l)>1)) call mpistop('cross-test: actual stencil violates 2:1')
          is_mixed=any(donors/=l)
          expected=0.d0
          if(at_wall.and.is_mixed) then
            counts(5)=counts(5)+1
            if(special_inner.and.v(1)<0) then
              call special_state(theta,phi,yy_patch,expected)
              ! Native spherical GLM closure removes the seeded radial field
              ! from its stencil. It leaves the outermost angular row intact.
              if(j>ixGlo2.and.j<ixGhi2.and.k>ixGlo3.and.k<ixGhi3) then
                expected(mag(:))=0.d0
                counts(7)=counts(7)+1
              end if
            else
              ri=merge(ixMlo1,ixMhi1,v(1)<0)
              expected=actual(ri,j,k,:)
              expected(mom(:))=matmul(frame,expected(mom(:)))
              expected(mag(:))=matmul(frame,expected(mag(:)))
            end if
          else
            do b=0,1
              do a=0,1
                idx=[max(0,min(nn(1)-1,v(1))),base(1)+a,base(2)+b]
                if(at_wall.and.special_inner.and.v(1)<0) then
                  ang=[xprobmin2,xprobmin3]+(idx(2:3)+0.5d0)*dq
                  call special_state(ang(1),ang(2),1-yy_patch,value)
                  value(mag(:))=0.d0
                else
                  call expected_donor(l,idx,donors(a,b),edges,value)
                end if
                expected=expected+value*merge(frac(1),1.d0-frac(1),a==1)*merge(frac(2),1.d0-frac(2),b==1)
                counts(donors(a,b)-l+3)=counts(donors(a,b)-l+3)+1
              end do
            end do
          end if
          expected(mom(:))=matmul(transpose(frame),expected(mom(:)))
          expected(mag(:))=matmul(transpose(frame),expected(mag(:)))
          if(any(.not.ieee_is_finite(actual(i,j,k,:)))) call mpistop('cross-test: invalid angular ghost')
          error=max(error,maxval(abs(expected-actual(i,j,k,:))))
          counts(1)=counts(1)+1
        end do
      end do
    end do
  end subroutine

  subroutine expected_donor(l,idx,donor_level,edges,value)
    integer,intent(in) :: l,idx(3),donor_level
    double precision,intent(in) :: edges(0:,:)
    double precision,intent(out) :: value(nw)
    double precision :: sample(nw),dv,vsum,dq(2),rlo,rhi,theta,phi
    integer :: i,j,k,q(3),n
    if(donor_level<=l) then
      call sample_at(l,idx,edges,value)
      return
    end if

    ! Analytic native child averages in the common Cartesian basis.
    value=0.d0
    vsum=0.d0
    n=l+1
    dq=[xprobmax2-xprobmin2,xprobmax3-xprobmin3]/([ng2(n),ng3(n)]*[block_nx2,block_nx3])
    do k=0,1
      do j=0,1
        do i=0,1
          q=2*idx+[i,j,k]
          call sample_at(n,q,edges,sample)
          theta=xprobmin2+q(2)*dq(1)
          rlo=edges(q(1),n)
          rhi=edges(q(1)+1,n)
          dv=(rhi-rlo)*(rhi*rhi+rlo*rhi+rlo*rlo)*(cos(theta)-cos(theta+dq(1)))
          sample(p_)=sample(p_)/(eos%gamma-1.d0)
          if(.not.mhd_internal_e) sample(p_)=sample(p_)+ &
            0.5d0*(sample(rho_)*sum(sample(mom(:))**2)+sum(sample(mag(:))**2))
          sample(mom(:))=sample(rho_)*sample(mom(:))
          value=value+dv*sample
          vsum=vsum+dv
        end do
      end do
    end do
    value=value/vsum
    value(mom(:))=value(mom(:))/value(rho_)
    if(.not.mhd_internal_e) value(p_)=value(p_)-0.5d0*(value(rho_)*sum(value(mom(:))**2)+sum(value(mag(:))**2))
    value(p_)=value(p_)*(eos%gamma-1.d0)
  end subroutine

  subroutine sample_at(level,cell,edges,w)
    integer,intent(in) :: level,cell(3)
    double precision,intent(in) :: edges(0:,:)
    double precision,intent(out) :: w(nw)
    double precision :: r,theta,phi
    r=0.5d0*(edges(cell(1),level)+edges(cell(1)+1,level))
    theta=xprobmin2+(cell(2)+0.5d0)*(xprobmax2-xprobmin2)/(ng2(level)*block_nx2)
    phi=xprobmin3+(cell(3)+0.5d0)*(xprobmax3-xprobmin3)/(ng3(level)*block_nx3)
    call profile(r,theta,phi,1-yy_patch,w)
  end subroutine
end module mod_usr
