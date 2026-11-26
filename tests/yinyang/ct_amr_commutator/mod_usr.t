module mod_usr
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use mod_mhd
  use mod_physics, only: phys_apply_ct_edges,phys_to_conserved
  use mod_functions_bfield, only: get_divb
  use mod_yinyang, only: yy_basis,yy_ct_stage_plan,yy_ct_ghost_plan, &
    yy_prepare_ct_edges,yy_exchange_ct_edges,yy_prepare_ct_faces,yy_exchange_ct_faces
  implicit none
  private
  public :: usr_init

  ! One Faraday update on the leaf mesh, independent of time integration, so
  ! the error is exactly the spatial commutator and has no RK or Riemann error.
  double precision, parameter :: probe_dt=0.013d0
  double precision, parameter :: r_min=1.18d0,r_max=1.32d0
  double precision, parameter :: phi_min_turn=0.165d0,phi_max_turn=0.225d0

contains

  subroutine usr_init()
    call set_coordinate_system('yinyang')
    call mhd_activate()
    usr_init_one_grid=>initial
    usr_refine_grid=>refine_yang_strip
    usr_before_main_loop=>measure
  end subroutine

  subroutine initial(ixI^L,ixO^L,w,x)
    integer,intent(in) :: ixI^L,ixO^L
    double precision,intent(in) :: x(ixI^S,ndim)
    double precision,intent(inout) :: w(ixI^S,nw)

    w(ixO^S,rho_)=1.d0
    w(ixO^S,p_)=1.d0
    w(ixO^S,mom(:))=0.d0
    w(ixO^S,mag(:))=0.d0
    block%ws=0.d0
    call phys_to_conserved(ixI^L,ixO^L,w,x)
  end subroutine

  subroutine refine_yang_strip(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    integer,intent(in) :: igrid,level,ixI^L,ixO^L
    double precision,intent(in) :: qt,w(ixI^S,nw),x(ixI^S,ndim)
    integer,intent(inout) :: refine,coarsen

    refine=-1
    coarsen=-1
    if(level/=1.or.yy_patch/=1.or.node(pig1_,igrid)/=1) return
    if(rnode(rpxmin2_,igrid)>0.15d0*2.d0*dpi) return
    if(rnode(rpxmax3_,igrid)<=0.131d0*2.d0*dpi) return
    refine=1
  end subroutine

  pure function electric_field(pos) result(value)
    double precision,intent(in) :: pos(3)
    double precision :: value(3),x,y,z
    x=pos(1)
    y=pos(2)
    z=pos(3)
    value=[0.2d0*y*y+0.1d0*z*x, -0.3d0*z*z+0.15d0*x*y, 0.4d0*x*x-0.12d0*y*z]
  end function

  pure function minus_curl(pos) result(value)
    double precision,intent(in) :: pos(3)
    double precision :: value(3)
    value=[-0.48d0*pos(3),0.7d0*pos(1),0.25d0*pos(2)]
  end function

  subroutine fill_native_edges(s,fE)
    type(state),intent(in) :: s
    double precision,intent(out) :: fE(ixG^T,3)
    integer :: lo(3),hi(3),c,i,j,k
    double precision :: offset(3),coord(3),frame(3,3),pos(3)

    fE=0.d0
    hi=[ixMhi1,ixMhi2,ixMhi3]
    do c=1,3
      lo=[ixMlo1,ixMlo2,ixMlo3]-1
      lo(c)=lo(c)+1
      offset=0.5d0
      offset(c)=0.d0
      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            coord=s%x(i,j,k,:)+offset*s%dx(i,j,k,:)
            call yy_basis(coord(2),coord(3),yy_patch,frame)
            pos=coord(1)*frame(:,1)
            fE(i,j,k,c)=probe_dt*s%dsC(i,j,k,c)*dot_product(electric_field(pos),frame(:,c))
          end do
        end do
      end do
    end do
  end subroutine

  ! Both moments use the physical surface weight r sin(theta).
  subroutine reference_theta_moments(center,delta,patch,mean,rms)
    double precision,intent(in) :: center(3),delta(3)
    integer,intent(in) :: patch
    double precision,intent(out) :: mean,rms
    double precision :: frame(3,3),coord(3),pos(3),value,w,total
    integer :: ir,ip
    double precision,parameter :: nodes(5)=[-0.906179845938664d0,-0.538469310105683d0, &
      0.d0,0.538469310105683d0,0.906179845938664d0]
    double precision,parameter :: weights(5)=[0.236926885056189d0,0.478628670499366d0, &
      0.568888888888889d0,0.478628670499366d0,0.236926885056189d0]

    mean=0.d0
    rms=0.d0
    total=0.d0
    do ip=1,5
      do ir=1,5
        coord=center
        coord(1)=center(1)+0.5d0*delta(1)*nodes(ir)
        coord(3)=center(3)+0.5d0*delta(3)*nodes(ip)
        call yy_basis(coord(2),coord(3),patch,frame)
        pos=coord(1)*frame(:,1)
        w=weights(ir)*weights(ip)*coord(1)*sin(coord(2))
        value=dot_product(minus_curl(pos),frame(:,2))
        mean=mean+w*value
        rms=rms+w*value**2
        total=total+w
      end do
    end do
    mean=mean/total
    rms=sqrt(rms/total)
  end subroutine

  double precision function target_area(s,i,j,k) result(area)
    type(state),intent(in) :: s
    integer,intent(in) :: i,j,k
    double precision :: dr,dphi,r,theta,phi,lo,hi,phi_width

    r=s%x(i,j,k,1)
    dr=s%dx(i,j,k,1)
    theta=s%x(i,j,k,2)+0.5d0*s%dx(i,j,k,2)
    phi=s%x(i,j,k,3)
    dphi=s%dx(i,j,k,3)
    lo=max(r_min,r-0.5d0*dr)
    hi=max(lo,min(r_max,r+0.5d0*dr))
    phi_width=max(0.d0,min(phi_max_turn*2.d0*dpi,phi+0.5d0*dphi) &
      -max(phi_min_turn*2.d0*dpi,phi-0.5d0*dphi))
    area=0.5d0*(hi-lo)*(hi+lo)*sin(theta)*phi_width
  end function

  subroutine select_targets(faces,edges)
    integer,allocatable,intent(out) :: faces(:,:),edges(:,:)
    integer :: count,pass,n,g,i,j,k,idx(3),level

    do pass=1,2
      count=0
      do n=1,igridstail
        g=igrids(n)
        if(node(pig2_,g)/=1) cycle
        j=ixMlo2-1
        level=node(plevel_,g)
        do k=ixMlo3,ixMhi3
          do i=ixMlo1,ixMhi1
            if(target_area(ps(g),i,j,k)<=0.d0) cycle
            if(level/=1+yy_patch) call mpistop('CT commutator receiver has wrong AMR level')
            count=count+1
            if(pass==1) cycle
            idx=[i,j,k]
            faces(:,count)=[g,idx,2]
            edges(:,4*count-3)=[g,idx,1]
            edges(:,4*count-2)=[g,i,j,k-1,1]
            edges(:,4*count-1)=[g,idx,3]
            edges(:,4*count)=[g,i-1,j,k,3]
          end do
        end do
      end do
      if(pass==1) allocate(faces(5,count),edges(5,4*count))
    end do
  end subroutine

  subroutine measure()
    type(ct_stage_flux),allocatable :: fluxes(:)
    type(state),allocatable,target :: probe(:)
    type(state),pointer :: probe_block
    type(yy_ct_ghost_plan) :: face_plans(2)
    type(yy_ct_stage_plan) :: edge_plans(2)
    integer,allocatable :: faces(:,:),edges(:,:)
    logical :: needed(max_blocks)
    double precision,allocatable :: face_first(:),native(:)
    integer :: n,g,p,i,j,k,record,idx(3),count(0:1),total_count(0:1)
    double precision :: area,exact,edge_first,error(4),sums(6,0:1),total(6,0:1)
    double precision :: local_max(5,0:1),global_max(5,0:1),divb(ixG^T)
    double precision :: edge_change(0:1),center(3),delta(3),reference_rms,expected_area
    character(len=80) :: filename

    if(.not.associated(phys_apply_ct_edges)) call mpistop('CT commutator needs the native UCT-HLL curl')
    call select_targets(faces,edges)
    call prepare_plans(faces,edges,face_plans,edge_plans,needed)
    allocate(fluxes(max_blocks),probe(max_blocks),native(size(faces,2)))
    do n=1,igridstail
      g=igrids(n)
      if(.not.needed(g)) cycle
      allocate(fluxes(g)%edge(ixG^T,3))
      allocate(probe(g)%ws,mold=ps(g)%ws)
      probe(g)%surfaceC=>ps(g)%surfaceC
      probe(g)%dvolume=>ps(g)%dvolume
      call fill_native_edges(ps(g),fluxes(g)%edge)
      probe(g)%ws=0.d0
      call phys_apply_ct_edges(ixG^LL,ixM^LL,fluxes(g)%edge,probe(g))
    end do

    do n=1,size(faces,2)
      g=faces(1,n)
      idx=faces(2:4,n)
      native(n)=probe(g)%ws(idx(1),idx(2),idx(3),2)
    end do

    call exchange_curl_then_transfer(face_plans,probe,faces,native,face_first)
    call exchange_transfer_then_curl(edge_plans,fluxes,edges,edge_change)

    sums=0.d0
    local_max=0.d0
    count=0
    count(yy_patch)=size(faces,2)
    do n=1,igridstail
      g=igrids(n)
      if(.not.needed(g)) cycle
      probe(g)%ws=0.d0
      call phys_apply_ct_edges(ixG^LL,ixM^LL,fluxes(g)%edge,probe(g))
      probe_block=>probe(g)
      block=>probe_block
      call get_divb(ps(g)%w,ixG^LL,ixM^LL,divb)
      block=>ps(g)
      local_max(5,yy_patch)=max(local_max(5,yy_patch),maxval(abs(divb(ixM^T)))/probe_dt)
    end do

    write(filename,'(a,i0,a)') 'samples-rank-',yy_world_rank,'.txt'
    open(newunit=record,file=trim(filename),status='replace')
    write(record,'(2a)') '# patch level block i j k r theta phi dr dphi ', &
      'clipped_area commutator_rate'
    do n=1,size(faces,2)
      g=faces(1,n)
      i=faces(2,n)
      j=faces(3,n)
      k=faces(4,n)
      area=target_area(ps(g),i,j,k)
      center=ps(g)%x(i,j,k,:)
      delta=ps(g)%dx(i,j,k,:)
      center(2)=center(2)+0.5d0*delta(2)
      call reference_theta_moments(center,delta,yy_patch,exact,reference_rms)
      edge_first=probe(g)%ws(i,j,k,2)
      ! Both CT fields are densities. Multiplying by A recovers flux, and
      ! dividing its increment by A and dt yields the measured rate.
      error=[edge_first-face_first(n),edge_first-probe_dt*exact, &
        face_first(n)-probe_dt*exact,native(n)-probe_dt*exact]/probe_dt
      if(.not.all(ieee_is_finite(error)).or.area<=0.d0) call mpistop('invalid CT commutator sample')
      write(record,'(6i8,7es24.15)') yy_patch,node(plevel_,g),g,i,j,k, &
        center,delta(1),delta(3),area,error(1)
      sums(1,yy_patch)=sums(1,yy_patch)+area
      sums(2:5,yy_patch)=sums(2:5,yy_patch)+area*error**2
      sums(6,yy_patch)=sums(6,yy_patch)+area*exact**2
      local_max(1:4,yy_patch)=max(local_max(1:4,yy_patch),abs(error))
    end do

    close(record)
    call MPI_ALLREDUCE(count,total_count,2,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierrmpi)
    call MPI_ALLREDUCE(sums,total,12,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierrmpi)
    call MPI_ALLREDUCE(local_max,global_max,10,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
    if(any(total_count==0).or.any(total(1,:)<=0.d0)) call mpistop('missing CT commutator direction')
    center=[0.5d0*(r_min+r_max),xprobmin2,(phi_min_turn+phi_max_turn)*dpi]
    delta=[r_max-r_min,0.d0,(phi_max_turn-phi_min_turn)*2.d0*dpi]
    expected_area=0.5d0*(r_max**2-r_min**2)*sin(center(2))*delta(3)
    if(maxval(abs(total(1,:)-expected_area))>1.d-12*expected_area) &
      call mpistop('CT commutator support area does not match analytic region')
    if(yy_world_rank==0) then
      do p=0,1
        call reference_theta_moments(center,delta,p,exact,reference_rms)
        write(*,'(a,2i8,13es24.15)') 'CT_COMMUTATOR ',p,total_count(p), &
          total(1,p),sqrt(total(2:6,p)/total(1,p)),global_max(:,p),edge_change(p),reference_rms
      end do
    end if
  end subroutine

  subroutine prepare_plans(faces,edges,face_plans,edge_plans,needed)
    integer,intent(in) :: faces(:,:),edges(:,:)
    type(yy_ct_ghost_plan),intent(out) :: face_plans(2)
    type(yy_ct_stage_plan),intent(out) :: edge_plans(2)
    logical,intent(out) :: needed(max_blocks)
    integer :: level,nfaces,nedges

    needed=.false.
    call mark_blocks(faces(1,:),needed)
    do level=1,2
      nfaces=0
      nedges=0
      if(yy_patch==level-1) then
        nfaces=size(faces,2)
        nedges=size(edges,2)
      end if
      call yy_prepare_ct_faces(face_plans(level),ps,faces(:,:nfaces),level=level)
      call yy_prepare_ct_edges(edge_plans(level),ps,edges(:,:nedges),level=level)
      if(size(face_plans(level)%donor_level_range,2)/=nfaces.or. &
        any(face_plans(level)%donor_level_range/=3-level)) &
        call mpistop('CT commutator face stencil is not entirely cross-level')
      if(size(edge_plans(level)%donor_level_range,2)/=nedges.or. &
        any(edge_plans(level)%donor_level_range/=3-level)) &
        call mpistop('CT commutator edge stencil is not entirely cross-level')
      call mark_blocks(face_plans(level)%faces%radial%request(1,:),needed)
      call mark_blocks(face_plans(level)%faces%tangent%request(1,:),needed)
      call mark_blocks(edge_plans(level)%edge%samples%request(1,:),needed)
    end do
  end subroutine

  subroutine mark_blocks(ids,needed)
    integer,intent(in) :: ids(:)
    logical,intent(inout) :: needed(max_blocks)
    integer :: n

    do n=1,size(ids)
      needed(ids(n))=.true.
    end do
  end subroutine

  subroutine exchange_curl_then_transfer(plans,probe,faces,native,values)
    type(yy_ct_ghost_plan),intent(inout) :: plans(2)
    type(state),intent(inout) :: probe(max_blocks)
    integer,intent(in) :: faces(:,:)
    double precision,intent(in) :: native(:)
    double precision,allocatable,intent(out) :: values(:)
    integer :: level,nt,n,g,i,j,k

    allocate(values(size(faces,2)))
    do level=1,2
      nt=size(plans(level)%receivers,2)
      call yy_exchange_ct_faces(plans(level),probe)
      ! Preserve original donor data for the other collective direction.
      do n=1,nt
        g=faces(1,n)
        i=faces(2,n)
        j=faces(3,n)
        k=faces(4,n)
        values(n)=probe(g)%ws(i,j,k,2)
        probe(g)%ws(i,j,k,2)=native(n)
      end do
    end do
  end subroutine

  subroutine exchange_transfer_then_curl(plans,fluxes,edges,change)
    type(yy_ct_stage_plan),intent(inout) :: plans(2)
    type(ct_stage_flux),intent(inout) :: fluxes(:)
    integer,intent(in) :: edges(:,:)
    double precision,intent(out) :: change(0:1)
    double precision,allocatable :: values(:),native(:)
    double precision :: local(0:1)
    integer :: level,nt,n,g,i,j,k,c

    allocate(values(size(edges,2)),native(size(edges,2)))
    do n=1,size(edges,2)
      g=edges(1,n)
      i=edges(2,n)
      j=edges(3,n)
      k=edges(4,n)
      c=edges(5,n)
      native(n)=fluxes(g)%edge(i,j,k,c)
    end do
    do level=1,2
      nt=size(plans(level)%receivers,2)
      call yy_exchange_ct_edges(plans(level),fluxes)
      do n=1,nt
        g=edges(1,n)
        i=edges(2,n)
        j=edges(3,n)
        k=edges(4,n)
        c=edges(5,n)
        values(n)=fluxes(g)%edge(i,j,k,c)
      end do
      do n=1,nt
        g=edges(1,n)
        i=edges(2,n)
        j=edges(3,n)
        k=edges(4,n)
        c=edges(5,n)
        fluxes(g)%edge(i,j,k,c)=native(n)
      end do
    end do

    local=0.d0
    do n=1,size(edges,2)
      g=edges(1,n)
      i=edges(2,n)
      j=edges(3,n)
      k=edges(4,n)
      c=edges(5,n)
      fluxes(g)%edge(i,j,k,c)=values(n)
      local(yy_patch)=max(local(yy_patch),abs(values(n)-native(n))/ &
        (probe_dt*ps(g)%dsC(i,j,k,c)))
    end do
    call MPI_ALLREDUCE(local,change,2,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
  end subroutine

end module
