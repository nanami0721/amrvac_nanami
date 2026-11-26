module mod_usr
  use mod_mhd
  use mod_functions_bfield, only: get_divb
  use mod_physics, only: phys_get_pthermal,phys_to_conserved
  use mod_yinyang, only: yy_basis,yy_ct_stage_plan,yy_runtime_ct_edges,yy_prepare_ct_edges, &
    yy_exchange_ct_edges,yy_ct_edge_is_owner,yy_ct_virtual_face_sources
  use mod_yinyang_staggered, only: yy_staggered_table,yy_staggered_max_sources
  use mod_connectivity, only: neighbor_type,neighbor_fine
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: usr_init
  logical :: use_spherical=.false.,check_wall=.false.
  double precision,parameter :: wall_component(3)=[0.d0,0.003d0,-0.004d0]
  integer :: layout=0, rank, stage_calls=0
contains
  subroutine usr_init()
    integer :: n
    namelist /ct_amr_test/ use_spherical,layout,check_wall
    do n=1,size(par_files)
      open(unitpar,file=trim(par_files(n)),status='old')
      read(unitpar,ct_amr_test,end=10)
10    close(unitpar)
    end do
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierrmpi)
    if(use_spherical) then
      call set_coordinate_system('spherical')
    else
      call set_coordinate_system('yinyang')
    end if
    call mhd_activate()
    usr_init_one_grid=>initial
    usr_init_vector_potential=>potential
    usr_refine_grid=>flags
    usr_after_refine=>after_regrid
    usr_process_adv_global=>after_advance
    usr_ct_stage=>check_edges
    if(check_wall) then
      usr_special_bc=>wall_fluid
      usr_set_electric_field=>wall_emf
    end if
  end subroutine

  subroutine initial(ixI^L,ixO^L,w,x)
    integer,intent(in) :: ixI^L,ixO^L
    double precision,intent(inout) :: w(ixI^S,nw)
    double precision,intent(in) :: x(ixI^S,ndim)
    integer :: ix^D
    double precision :: frame(3,3),v(3)
    w(ixO^S,rho_)=1.d0
    w(ixO^S,p_)=1.d0
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      call yy_basis(x(ix^D,2),x(ix^D,3),yy_patch,frame)
      v=[0.01d0,-0.02d0,0.03d0]
      w(ix^D,mom(:))=matmul(transpose(frame),v)
    {end do\}
    call b_from_vector_potential(block%ixGs^L,ixI^L,ixO^L,block%ws,x)
    call mhd_face_to_center(ixO^L,block)
    call phys_to_conserved(ixI^L,ixO^L,w,x)
  end subroutine

  subroutine potential(ixI^L,ixO^L,x,a,idir)
    integer,intent(in) :: ixI^L,ixO^L,idir
    double precision,intent(in) :: x(ixI^S,ndim)
    double precision,intent(out) :: a(ixI^S)
    integer :: ix^D
    double precision :: frame(3,3),r(3),cart(3)
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      call yy_basis(x(ix^D,2),x(ix^D,3),yy_patch,frame)
      r=x(ix^D,1)*frame(:,1)
      cart=0.5d0*[-0.01d0*r(3)-0.03d0*r(2), &
        0.03d0*r(1)-0.02d0*r(3),0.02d0*r(2)+0.01d0*r(1)]
      a(ix^D)=dot_product(cart,frame(:,idir))
    {end do\}
  end subroutine

  subroutine flags(g,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    integer,intent(in) :: g,level,ixI^L,ixO^L
    double precision,intent(in) :: qt,w(ixI^S,nw),x(ixI^S,ndim)
    integer,intent(inout) :: refine,coarsen
    refine=-1
    coarsen=-1
    if((layout==0.or.layout==3).and.it>=4) then
      if(level>1) coarsen=1
      return
    end if
    if(level>=refine_max_level) return
    select case(layout)
    case(0)
      if(it<1.or.yy_patch/=0) return
      if(node(pig1_,g)==2.and.node(pig2_,g)==2.and.node(pig3_,g)==6) refine=1
    case(1)
      refine=1
    case(2)
      if(yy_patch==1) refine=1
    case(3)
      if(it<1.or.yy_patch/=0) return
      if(node(pig1_,g)/=1.or.node(pig2_,g)/=1) return
      if(rnode(rpxmin3_,g)<=0.015d0.and.rnode(rpxmax3_,g)>0.015d0) refine=1
    end select
  end subroutine

  subroutine check_edges(qt,qdt,sct,s,fluxes)
    double precision,intent(in) :: qt,qdt
    type(state),intent(in) :: sct(max_blocks),s(max_blocks)
    type(ct_stage_flux),intent(inout) :: fluxes(:)
    integer :: n,g
    stage_calls=stage_calls+1
    do n=1,igridstail_active
      g=igrids_active(n)
      if(.not.all(ieee_is_finite(fluxes(g)%edge))) call mpistop('CT AMR test: nonfinite RK edge integral')
    end do
    if(check_wall.and.it==2.and.modulo(stage_calls-1,3)==0) call inspect_wall(qdt,sct,fluxes)
  end subroutine

  subroutine wall_fluid(qdt,qt,ixI^L,ixO^L,iB,w,x)
    integer,intent(in) :: ixI^L,ixO^L,iB
    double precision,intent(in) :: qdt,qt,x(ixI^S,ndim)
    double precision,intent(inout) :: w(ixI^S,nw)
    if(iB/=1) call mpistop('CT AMR test: unexpected special boundary')
    w(ixO^S,rho_)=1.d0
    w(ixO^S,mom(:))=0.d0
    w(ixO^S,p_)=1.d0
    call phys_to_conserved(ixI^L,ixO^L,w,x)
  end subroutine

  subroutine wall_emf(ixI^L,ixO^L,qt,qdt,edge,s)
    integer,intent(in) :: ixI^L,ixO^L
    double precision,intent(in) :: qt,qdt
    double precision,intent(inout) :: edge(ixI^S,sdim:3)
    type(state) :: s
    integer :: i,j,k
    if(.not.s%is_physical_boundary(1)) return
    i=ixOmin1-1
    do k=ixOmin3-1,ixOmax3
      do j=ixOmin2,ixOmax2
        edge(i,j,k,2)=qdt*s%dsC(i,j,k,2)*wall_component(2)
      end do
    end do
    do k=ixOmin3,ixOmax3
      do j=ixOmin2-1,ixOmax2
        edge(i,j,k,3)=qdt*s%dsC(i,j,k,3)*wall_component(3)
      end do
    end do
  end subroutine

  subroutine inspect_wall(qdt,s,fluxes)
    double precision,intent(in) :: qdt
    type(state),intent(in) :: s(max_blocks)
    type(ct_stage_flux),intent(in) :: fluxes(:)
    type(ct_stage_flux),allocatable :: copy(:)
    type(yy_ct_stage_plan) :: plan
    integer,allocatable :: receivers(:,:),local(:,:),leaves(:,:),counts(:),offsets(:)
    integer :: n,g,i,j,k,l,nt,pass,nn(3),bn(3),orig(3),lo(3),hi(3),pos(3),p,level_max
    integer :: checked,all_checked
    double precision :: diff,diff_max,error,error_max

    ! This preview is observational. The production stage still receives the
    ! untouched original workspace and performs the identical coupling once.
    allocate(copy(max_blocks),local(5,igridstail),counts(yy_world_size),offsets(yy_world_size))
    bn=[block_nx1,block_nx2,block_nx3]
    lo=[ixMlo1,ixMlo2,ixMlo3]
    hi=[ixMhi1,ixMhi2,ixMhi3]
    do n=1,igridstail
      g=igrids(n)
      copy(g)%edge=fluxes(g)%edge
      local(:,n)=[node(plevel_,g),node(pig1_,g),node(pig2_,g),node(pig3_,g),yy_patch]
    end do
    call MPI_ALLGATHER(5*igridstail,1,MPI_INTEGER,counts,1,MPI_INTEGER,MPI_COMM_WORLD,ierrmpi)
    offsets(1)=0
    do p=2,yy_world_size
      offsets(p)=offsets(p-1)+counts(p-1)
    end do
    allocate(leaves(5,sum(counts)/5))
    call MPI_ALLGATHERV(local,5*igridstail,MPI_INTEGER,leaves,counts,offsets,MPI_INTEGER,MPI_COMM_WORLD,ierrmpi)

    call yy_runtime_ct_edges(s,copy)
    call inspect_virtual_flux(s)
    diff=0.d0
    do n=1,igridstail
      g=igrids(n)
      if(node(pig1_,g)/=1) cycle
      i=ixMlo1-1
      do k=ixMlo3-1,ixMhi3
        do j=ixMlo2,ixMhi2
          diff=max(diff,abs(copy(g)%edge(i,j,k,2)-fluxes(g)%edge(i,j,k,2)))
        end do
      end do
      do k=ixMlo3,ixMhi3
        do j=ixMlo2-1,ixMhi2
          diff=max(diff,abs(copy(g)%edge(i,j,k,3)-fluxes(g)%edge(i,j,k,3)))
        end do
      end do
    end do
    call MPI_ALLREDUCE(diff,diff_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
    if(diff_max>1.d-18) call mpistop('CT AMR test: prescribed wall EMF changed at the seam')

    checked=0
    error=0.d0
    level_max=maxval(leaves(1,:))
    do l=1,level_max
      nn=bn*[ng1(l),ng2(l),ng3(l)]
      do pass=1,2
        nt=0
        do n=1,igridstail
          g=igrids(n)
          if(node(plevel_,g)/=l.or.node(pig1_,g)/=1) cycle
          orig=bn*([node(pig1_,g),node(pig2_,g),node(pig3_,g)]-1)
          do k=lo(3)-1,hi(3)
            do j=lo(2)-1,hi(2)
              pos=orig+[lo(1),j,k]-lo
              if(pos(2)/=-1.and.pos(2)/=nn(2)-1.and.pos(3)/=-1.and.pos(3)/=nn(3)-1) cycle
              if(.not.yy_ct_edge_is_owner(g,[lo(1),j,k],1)) cycle
              nt=nt+1
              if(pass==2) receivers(:,nt)=[g,lo(1),j,k,1]
            end do
          end do
        end do
        if(pass==1) allocate(receivers(5,nt))
      end do
      call yy_prepare_ct_edges(plan,s,receivers,radial_extrapolate=.true.,level=l)
      call yy_exchange_ct_edges(plan,copy)
      call inspect_wall_donors(plan%edge%samples,l,qdt,leaves,bn,checked,error)
      deallocate(receivers)
    end do
    call MPI_ALLREDUCE(checked,all_checked,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierrmpi)
    call MPI_ALLREDUCE(error,error_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
    if(rank==0) write(*,'(a,1x,i0,2(1x,es13.5))') 'CT AMR wall pairs',all_checked,diff_max,error_max
    if(all_checked==0.or.error_max>1.d-14) call mpistop('CT AMR test: wall donor fine-edge sum failed')
  end subroutine

  subroutine inspect_wall_donors(table,l,qdt,leaves,bn,checked,error)
    type(yy_staggered_table),intent(in) :: table
    integer,intent(in) :: l,leaves(:,:),bn(3)
    double precision,intent(in) :: qdt
    integer,intent(inout) :: checked
    double precision,intent(inout) :: error
    integer,allocatable :: tag_in(:,:),tag_out(:,:)
    integer :: n,g,t,q,term,slots(2),first(5),second(5),c,coarse(3),low(3),high(3),delta(3)
    double precision :: length,integral
    logical :: touches_coarse
    allocate(tag_in(5,size(table%request,2)),tag_out(5,size(table%recv)))
    do n=1,size(table%request,2)
      g=table%request(1,n)
      tag_in(:,n)=[node(plevel_,g),bn*([node(pig1_,g),node(pig2_,g),node(pig3_,g)]-1)+ &
        table%request(2:4,n),table%request(5,n)]
    end do
    call MPI_ALLTOALLV(tag_in,5*table%rc,5*table%rd,MPI_INTEGER,tag_out,5*table%sc,5*table%sd, &
      MPI_INTEGER,MPI_COMM_WORLD,ierrmpi)
    do t=1,size(table%covered)
      do q=1,table%stencil
        term=table%term_end(q-1,t)+1
        if(table%term_end(q,t)/=term+1) cycle
        if(maxval(abs(table%coefficient(:,q,t)))<=1.d-13) cycle
        slots=table%term_slot(term:term+1)
        first=tag_out(:,slots(1))
        second=tag_out(:,slots(2))
        c=first(5)
        if(c<2.or.first(2)/=-1.or.second(2)/=-1.or.second(5)/=c) cycle
        if(first(1)/=l+1.or.second(1)/=l+1) cycle
        delta=0
        delta(c)=1
        if(any(second(2:4)-first(2:4)/=delta)) cycle
        coarse=(first(2:4)-1)/2
        coarse(c)=first(c+1)/2
        touches_coarse=.false.
        do n=1,size(leaves,2)
          if(leaves(1,n)/=l.or.leaves(5,n)==yy_patch) cycle
          low=bn*(leaves(2:4,n)-1)
          high=low+bn-1
          low=low-1
          low(c)=low(c)+1
          if(all(coarse>=low).and.all(coarse<=high)) then
            touches_coarse=.true.
            exit
          end if
        end do
        if(.not.touches_coarse) cycle
        length=xprobmin1*dx(c,l)
        if(c==3) length=length*sin(xprobmin2+(coarse(2)+1)*dx(2,l))
        integral=sum(table%recv(slots))
        error=max(error,abs(integral-qdt*length*wall_component(c)))
        error=max(error,length*maxval(abs(table%term_weight(term:term+1)*length-1.d0)))
        checked=checked+1
      end do
    end do
  end subroutine

  ! At each native coarse/fine face, check all eight virtual children of
  ! the coarse boundary cell. The incident existing fine block keeps its own
  ! active faces. This check isolates the coarse donor's prolongation map.
  subroutine inspect_virtual_flux(s)
    type(state),intent(in) :: s(max_blocks)
    integer :: n,g,d,side,l,a,b,c,j,lo(3),hi(3),origin(3),idx(3),child(3),axis(3),bn(3)
    integer :: checked,total
    double precision :: parent,net,error,maximum,upper,lower
    logical :: local
    lo=[ixMlo1,ixMlo2,ixMlo3]
    hi=[ixMhi1,ixMhi2,ixMhi3]
    bn=[block_nx1,block_nx2,block_nx3]
    checked=0
    error=0.d0
    do n=1,igridstail
      g=igrids(n)
      l=node(plevel_,g)
      origin=bn*([node(pig1_,g),node(pig2_,g),node(pig3_,g)]-1)
      do d=1,3
        do side=-1,1,2
          axis=0
          axis(d)=side
          if(neighbor_type(axis(1),axis(2),axis(3),g)/=neighbor_fine) cycle
          idx=(lo+hi)/2
          idx(d)=merge(lo(d),hi(d),side==-1)
          parent=0.d0
          do j=1,3
            axis=0
            axis(j)=1
            parent=parent+s(g)%ws(idx(1),idx(2),idx(3),j)*s(g)%surfaceC(idx(1),idx(2),idx(3),j) &
              -s(g)%ws(idx(1)-axis(1),idx(2)-axis(2),idx(3)-axis(3),j)* &
              s(g)%surfaceC(idx(1)-axis(1),idx(2)-axis(2),idx(3)-axis(3),j)
          end do
          do c=0,1
            do b=0,1
              do a=0,1
                child=2*(idx-lo+origin)+[a,b,c]
                net=0.d0
                local=.true.
                do j=1,3
                  axis=0
                  axis(j)=1
                  upper=virtual_flux(s,l+1,child,j,local)
                  lower=virtual_flux(s,l+1,child-axis,j,local)
                  net=net+upper-lower
                end do
                if(.not.local) cycle
                checked=checked+1
                error=max(error,abs(net-parent/8.d0))
              end do
            end do
          end do
        end do
      end do
    end do
    call MPI_ALLREDUCE(checked,total,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierrmpi)
    call MPI_ALLREDUCE(error,maximum,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
    if(rank==0) write(*,'(a,1x,i0,1x,es13.5)') 'CT AMR virtual child fluxes',total,maximum
    if(total==0.or.maximum>2.d-14) call mpistop('CT AMR test: donor-aligned child flux closure failed')
  end subroutine

  double precision function virtual_flux(s,l,idx,c,local) result(flux)
    type(state),intent(in) :: s(max_blocks)
    integer,intent(in) :: l,idx(3),c
    logical,intent(inout) :: local
    integer,parameter :: capacity=yy_staggered_max_sources
    integer :: ranks(capacity),requests(5,capacity),n,m,g,pt(3)
    double precision :: weights(capacity),r0,r1,t0,t1,area
    flux=0.d0
    if(.not.local) return
    call yy_ct_virtual_face_sources(l,yy_patch,c,idx,ranks,requests,weights,n)
    if(n==0) call mpistop('CT AMR test: missing parent flux sources')
    if(any(ranks(1:n)/=yy_world_rank+1)) then
      local=.false.
      return
    end if
    do m=1,n
      g=requests(1,m)
      pt=requests(2:4,m)+[ixMlo1,ixMlo2,ixMlo3]
      flux=flux+weights(m)*s(g)%ws(pt(1),pt(2),pt(3),c)
    end do
    r0=xprobmin1+idx(1)*dx(1,l)
    r1=r0+dx(1,l)
    t0=xprobmin2+idx(2)*dx(2,l)
    t1=t0+dx(2,l)
    select case(c)
    case(1);area=r1**2*(cos(t0)-cos(t1))*dx(3,l)
    case(2);area=0.5d0*(r1**2-r0**2)*sin(t1)*dx(3,l)
    case(3);area=0.5d0*(r1**2-r0**2)*dx(2,l)
    end select
    flux=flux*area
  end function

  subroutine after_regrid(ncoarse,nfine)
    integer,intent(in) :: ncoarse,nfine
    call examine('regrid',it)
  end subroutine

  subroutine after_advance(iit,qt)
    integer,intent(in) :: iit
    double precision,intent(in) :: qt
    call examine('advance',iit+1)
  end subroutine

  subroutine examine(event,iit)
    character(len=*),intent(in) :: event
    integer,intent(in) :: iit
    integer :: n,g,counts(nlevelshi,0:1),total(nlevelshi,0:1)
    double precision :: divb(ixG^T),pressure(ixG^T),maxdiv,global_div,minp,global_p,ghost_div,global_ghost
    counts=0;maxdiv=0.d0;ghost_div=0.d0;minp=huge(1.d0)
    do n=1,igridstail
      g=igrids(n);block=>ps(g)
      counts(block%level,yy_patch)=counts(block%level,yy_patch)+1
      call get_divb(block%w,ixG^LL,ixG^LL,divb)
      call phys_get_pthermal(block%w,block%x,ixG^LL,ixG^LL,pressure)
      if(.not.all(ieee_is_finite(divb))) call mpistop('CT AMR test: nonfinite face divergence')
      if(.not.all(ieee_is_finite(pressure))) call mpistop('CT AMR test: nonfinite pressure')
      maxdiv=max(maxdiv,maxval(abs(divb(ixM^T))))
      ghost_div=max(ghost_div,maxval(abs(divb)))
      minp=min(minp,minval(pressure))
    end do
    call MPI_ALLREDUCE(counts,total,size(counts),MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierrmpi)
    call MPI_ALLREDUCE(maxdiv,global_div,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
    call MPI_ALLREDUCE(ghost_div,global_ghost,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
    call MPI_ALLREDUCE(minp,global_p,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierrmpi)
    if(rank==0) write(*,'(a,1x,a,*(1x,g0))') 'CT AMR',event,iit,stage_calls, &
      total(1:refine_max_level,0),total(1:refine_max_level,1),global_div,global_ghost,global_p
    if(global_div>2.d-12.or.global_p<=0.d0) call mpistop('CT AMR test: active divergence or pressure bound failed')
  end subroutine
end module
