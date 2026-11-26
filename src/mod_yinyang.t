!> Yin–Yang overset coupling on native spherical meshes.
module mod_yinyang
  use mod_yinyang_staggered, only: yy_staggered_table, yy_staggered_build, yy_face_location, yy_edge_location, &
    yy_staggered_offset, yy_staggered_coordinates, yy_staggered_bracket
  use mod_yinyang_ct_coupling, only: yy_ct_edge_plan, yy_ct_edge_build, yy_ct_edge_exchange, &
    yy_ct_face_plan, yy_ct_face_build, yy_ct_face_exchange, yy_ct_close_ghosts
  use mod_yinyang_patch_sync, only: yy_patch_copy_plan, yy_patch_copy_build, yy_patch_copy_exchange, &
    yy_sync_edges, yy_sync_faces, yy_sync_extended_faces
  use mod_global_parameters
  use, intrinsic :: iso_fortran_env, only: int64
  {^IFTHREED
  use mod_yinyang_glm_amr, only: radial_edge, leaf_directory, build_directory, find_leaf, covering_leaf, &
    domain_cells, key_level, key_index, key_rank, key_grid, check_level_pair
  }
  use mod_yinyang_glm_amr, only: yy_glm_invalidate,yy_glm_report,yy_basis=>yy_glm_basis
  use mod_physics, only: phys_to_primitive, phys_to_conserved, physics_type, &
    phys_ct_sync_center,phys_get_ct_edges,phys_apply_ct_edges
  use mod_comm_lib, only: mpistop
  use mod_forest, only: particle_mesh_epoch
  implicit none
  private
  ! Owned here: only input parsing and Yin-Yang initialization use this option.
  logical,public :: yy_ct_face_moments=.false.
  public :: yy_validate, yy_exchange, yy_basis, yy_position, yy_checkpoint, yy_checkpoint_begin, yy_report
  public :: yy_invalidate_after_regrid, yy_refresh_glm_mode
  logical,save :: glm_mode_known=.false., glm_mode_mixed=.false.
  public :: yy_build_staggered_table, yy_prepare_ct_edges, yy_exchange_ct_edges, yy_ct_stage_plan
  public :: yy_ct_ghost_plan, yy_prepare_ct_ghosts, yy_fill_ct_ghosts, yy_prepare_ct_faces, yy_exchange_ct_faces
  type yy_ct_ghost_plan
    type(yy_ct_face_plan) :: faces
    integer,allocatable :: receivers(:,:),donor_level_range(:,:)
    double precision,allocatable :: radial(:),tangent(:),values(:)
    integer(kind=8) :: mesh_epoch=-1_8
  end type
  public :: yy_ct_patch_plan, yy_prepare_ct_patch, yy_sync_ct_edges, yy_ordered_ct_edges
  public :: yy_sync_ct_faces, yy_ordered_ct_ghosts, yy_ct_edge_is_owner
  public :: yy_complete_ct_faces, yy_sync_ct_state, yy_ct_radial_faces
  type yy_ct_patch_plan
    type(yy_patch_copy_plan) :: edge,interior_faces,extended_faces
    integer(kind=8) :: mesh_epoch=-1_8
  end type
  type yy_ct_stage_plan
    type(yy_ct_edge_plan) :: edge
    integer,allocatable :: receivers(:,:),donor_level_range(:,:)
    double precision,allocatable :: packed(:), received(:)
    integer(kind=8) :: mesh_epoch=-1_8
  end type
  logical,save :: special_inner_ct=.false.
  public :: yy_runtime_ct_edges, yy_runtime_ct_ghosts, yy_ct_radial_cells
  type yy_runtime_ct_plan
    type(yy_ct_patch_plan) :: patch
    type(yy_ct_stage_plan) :: edges
    type(yy_ct_ghost_plan) :: ghosts
  end type
  type(yy_runtime_ct_plan),allocatable,save :: runtime_ct(:)
  logical,save :: sync_edges_before_exchange=.true.
  integer(kind=8),save :: runtime_ct_epoch=-1_8
  {^IFTHREED
  public :: yy_ct_balance, yy_ct_virtual_face_sources
  type(leaf_directory),save :: ct_directory
  integer(kind=8),save :: ct_directory_epoch=-1_8
  integer,save :: ct_build_level,ct_uniform_level(0:1)
  }
  integer(kind=8) :: epoch=-1_8
  integer :: nr=0, ns=0, nt=0, raw_nr=0
  integer(kind=8) :: exchange_calls=0_8
  integer, allocatable :: rec(:,:), slot(:,:), request(:,:), sc(:),rc(:),sd(:),rd(:)
  integer, allocatable :: handles(:), statuses(:,:)
  integer, allocatable :: shead(:),snext(:),rhead(:),rnext(:)
  double precision, allocatable :: weights(:,:), rbasis(:,:,:), sbasis(:,:,:)
  double precision, allocatable :: sendbuf(:,:), recvbuf(:,:), scratch(:^D&,:)
  double precision :: pack_seconds=0.d0, wait_seconds=0.d0, unpack_seconds=0.d0, build_seconds=0.d0
contains
  ! Build an explicitly owned table from the real fixed-grid AMRVAC layout.
  ! Caller owns its lifetime and must rebuild after mesh/ownership changes.
  ! This does not activate CT or write state/ghost faces.
  subroutine yy_build_staggered_table(table,s,family,points)
    type(yy_staggered_table),intent(inout) :: table
    type(state),intent(in) :: s(max_blocks)
    integer,intent(in) :: family
    double precision,intent(in) :: points(:,:)
    integer :: nn(3),bn(3)
    integer,allocatable :: directory(:,:)
    double precision,allocatable :: radial(:)
    double precision :: angular(2,2)
    call yy_staggered_mesh(s,bn,nn,radial,angular,directory)
    call yy_staggered_build(table,MPI_COMM_WORLD,yy_patch,family,bn,nn,radial,angular,directory,points)
  end subroutine

  subroutine yy_staggered_mesh(s,bn,nn,rmax,angular,directory,level)
    type(state),intent(in) :: s(max_blocks)
    integer,intent(out) :: bn(3),nn(3)
    double precision,allocatable,intent(out) :: rmax(:)
    double precision,intent(out) :: angular(2,2)
    integer,allocatable,intent(out) :: directory(:,:)
    integer,optional,intent(in) :: level
    {^IFTHREED
    integer :: n,g,p,i,gi,l,entries,nb(3)
    integer,allocatable :: local(:,:),counts(:),displs(:)
    double precision,allocatable :: lower(:),upper(:),rmin(:)
    double precision :: radius
    if(.not.yinyang_enabled.or.lb_automatic) call mpistop('invalid staggered Yin-Yang layout')
    l=1
    if(present(level)) then
      l=level
    else if(refine_max_level/=1) then
      call mpistop('explicit level required for adaptive staggered Yin-Yang layout')
    end if
    if(stretched_dim(2).or.stretched_dim(3)) call mpistop('staggered table needs uniform angular mesh')
    bn=[block_nx1,block_nx2,block_nx3]
    nb=[ng1(l),ng2(l),ng3(l)]
    nn=nb*bn
    allocate(local(6,igridstail),counts(yy_world_size),displs(yy_world_size))
    allocate(rmax(0:nn(1)))
    if(refine_max_level==1) then
      allocate(lower(0:nn(1)),upper(0:nn(1)),rmin(0:nn(1)))
    end if
    if(allocated(lower)) then
      lower=huge(1.d0);upper=-huge(1.d0)
    end if
    entries=0
    do n=1,igridstail
      g=igrids(n)
      if(node(plevel_,g)/=l) cycle
      entries=entries+1
      local(:,entries)=[node(pig1_,g)-1,node(pig2_,g)-1,node(pig3_,g)-1,g,yy_world_rank,yy_patch]
      if(.not.allocated(lower)) cycle
      do i=ixMlo1-1,ixMhi1
        gi=(node(pig1_,g)-1)*block_nx1+i-ixMlo1+1
        radius=s(g)%x(i,ixMlo2,ixMlo3,1)+0.5d0*s(g)%dx(i,ixMlo2,ixMlo3,1)
        lower(gi)=min(lower(gi),radius);upper(gi)=max(upper(gi),radius)
      end do
    end do
    if(allocated(lower)) then
      call MPI_ALLREDUCE(lower,rmin,nn(1)+1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierrmpi)
      call MPI_ALLREDUCE(upper,rmax,nn(1)+1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
      if(any(abs(rmax-rmin)>1.d-12*max(1.d0,abs(rmax)))) &
        call mpistop('staggered table radial meshes disagree or are incomplete')
    else
      do i=0,nn(1)
        rmax(i)=radial_edge(l,int(i,int64))
      end do
    end if
    call MPI_ALLGATHER(6*entries,1,MPI_INTEGER,counts,1,MPI_INTEGER,MPI_COMM_WORLD,ierrmpi)
    displs(1)=0
    do p=2,yy_world_size
      displs(p)=displs(p-1)+counts(p-1)
    end do
    allocate(directory(6,sum(counts)/6))
    call MPI_ALLGATHERV(local,6*entries,MPI_INTEGER,directory,counts,displs,MPI_INTEGER,MPI_COMM_WORLD,ierrmpi)
    angular(:,1)=[xprobmin2,xprobmin3];angular(:,2)=[xprobmax2,xprobmax3]
    }
    {^IFONED call mpistop('staggered Yin-Yang table requires 3D') }
    {^IFTWOD call mpistop('staggered Yin-Yang table requires 3D') }
  end subroutine

  {^IFTHREED
  ! Check proposed leaves before committing a regrid. The projection and the
  ! staggered lattice brackets are shared with the exchange plan construction.
  subroutine yy_ct_balance(dir,want)
    type(leaf_directory),intent(in) :: dir
    integer,intent(inout) :: want(:)
    integer :: n,l,nn(3),bn(3),origin(3),i,j,k,c,boxlo(3),boxhi(3),idx(3),lo(3)
    integer :: previous_level
    double precision,allocatable :: radial(:)
    double precision :: angular(2,2),step(2)
    logical :: sides(2,2)
    logical :: mask(-nghostcells-1:block_nx1+nghostcells-1, &
                    -nghostcells-1:block_nx2+nghostcells-1, &
                    -nghostcells-1:block_nx3+nghostcells-1,3)
    bn=[block_nx1,block_nx2,block_nx3]
    angular(:,1)=[xprobmin2,xprobmin3]
    angular(:,2)=[xprobmax2,xprobmax3]
    previous_level=0
    do n=1,size(dir%leaf,2)
      if(dir%leaf(key_rank,n)/=yy_world_rank) cycle
      l=dir%leaf(key_level,n)
      nn=int(domain_cells(l))
      origin=(dir%leaf(key_index:key_index+2,n)-1)*bn
      call ct_face_bounds(l,origin,bn,boxlo,boxhi,sides)
      if(.not.any(sides)) cycle
      if(previous_level/=l) then
        if(allocated(radial)) deallocate(radial)
        allocate(radial(0:nn(1)))
        do i=0,nn(1)
          radial(i)=radial_edge(l,int(i,int64))
        end do
        step=(angular(:,2)-angular(:,1))/nn(2:3)
        previous_level=l
      end if
      call ct_face_footprint(boxlo,boxhi,sides,mask)
      do c=1,3
        do k=-nghostcells-1,bn(3)+nghostcells-1
          do j=-nghostcells-1,bn(2)+nghostcells-1
            do i=-nghostcells-1,bn(1)+nghostcells-1
              if(mask(i,j,k,c)) call check_target(yy_face_location,c,origin+[i,j,k])
            end do
          end do
        end do
        lo=-1
        lo(c)=0
        do k=lo(3),bn(3)-1
          do j=lo(2),bn(2)-1
            do i=lo(1),bn(1)-1
              idx=origin+[i,j,k]
              if(.not.ct_seam_edge(idx,c,nn)) cycle
              if(special_inner_ct.and.c/=1.and.idx(1)==-1) cycle
              call check_target(yy_edge_location,c,idx)
            end do
          end do
        end do
      end do
    end do
  contains
    subroutine check_target(family,c,receiver)
      integer,intent(in) :: family,c,receiver(3)
      integer :: component,ncomponents,base(3),idx(3),lower(3),upper(3),a,b,h
      double precision :: offset(3),point(3),sph(3),xyz(3),fraction(3)
      logical :: valid,moments,edge,shared(3)
      edge=family==yy_edge_location
      moments=.not.edge.and.yy_ct_face_moments
      ncomponents=3
      if(.not.edge.and.c==1.and..not.moments) ncomponents=1
      call yy_staggered_offset(family,c,offset)
      point(1)=radial(receiver(1)+1)
      if(offset(1)==0.5d0) point(1)=0.5d0*(point(1)+radial(receiver(1)))
      point(2:3)=angular(:,1)+(receiver(2:3)+offset(2:3))*step
      call yy_staggered_coordinates(point,yy_patch,sph,xyz)
      do component=1,ncomponents
        call yy_staggered_offset(family,component,offset)
        call yy_staggered_bracket(sph,offset,nn,radial,angular,step,edge.and.component==1,moments,base,fraction,valid)
        if(.not.valid) call mpistop('projected CT target lacks angular or radial support')
        shared=edge
        shared(component)=.not.edge
        do h=0,1
          do b=0,1
            do a=0,1
              idx=base+[a,b,h]
              lower=max(0,idx)
              upper=lower
              ! An edge/face on a block plane also has a native copy in
              ! the adjacent block. Fine collinear segments can own EMF.
              where(shared.and.modulo(idx+1,bn)==0.and.idx>=0) upper=upper+1
              if(edge) then
                ! A virtual coarse edge may sum fine segments whose owner
                ! is on the opposite side of a fine block plane.
                where(shared.and.modulo(2*(idx+1),bn)==0.and.idx>=0) upper=lower+1
              end if
              upper=min(upper,nn-1)
              call check_sites(lower,upper)
            end do
          end do
        end do
      end do
    end subroutine

    subroutine check_sites(lower,upper)
      integer,intent(in) :: lower(3),upper(3)
      integer :: i,j,k,other
      do k=lower(3),upper(3)
        do j=lower(2),upper(2)
          do i=lower(1),upper(1)
            other=covering_leaf(dir,l,int([i,j,k],int64),1-yy_patch)
            call check_level_pair(dir,n,other,want)
          end do
        end do
      end do
    end subroutine
  end subroutine

  pure logical function ct_seam_edge(idx,c,nn) result(seam)
    integer,intent(in) :: idx(3),c,nn(3)
    seam=(c/=2.and.(idx(2)==-1.or.idx(2)==nn(2)-1)).or. &
         (c/=3.and.(idx(3)==-1.or.idx(3)==nn(3)-1))
  end function

  ! Active radial cells plus all native halo cells inside the physical shell.
  ! At mixed junctions the whole angular-outside footprint needs an overset
  ! donor, including halos alongside a coarse/fine radial or angular neighbor.
  subroutine ct_face_bounds(level,origin,bn,lo,hi,sides)
    integer,intent(in) :: level,origin(3),bn(3)
    integer,intent(out) :: lo(3),hi(3)
    logical,intent(out) :: sides(2,2)
    integer :: nn(3)
    nn=int(domain_cells(level))
    lo=max(-nghostcells,-origin)
    hi=min(bn+nghostcells-1,nn-1-origin)
    sides(1,:)=origin(2:3)==0
    sides(2,:)=origin(2:3)+bn(2:3)==nn(2:3)
  end subroutine

  subroutine ct_face_footprint(lo,hi,sides,mask)
    integer,intent(in) :: lo(3),hi(3)
    logical,intent(in) :: sides(2,2)
    logical,intent(out) :: mask(-nghostcells-1:,-nghostcells-1:,-nghostcells-1:,:)
    integer :: i,j,k,c,side,low(3),high(3),idx(3),fl(3)
    mask=.false.
    if(.not.any(sides)) return
    low=lo
    high=hi
    where(sides(1,:)) low(2:3)=lo(2:3)-nghostcells
    where(sides(2,:)) high(2:3)=hi(2:3)+nghostcells
    do k=low(3),high(3)
      do j=low(2),high(2)
        if(j>=lo(2).and.j<=hi(2).and.k>=lo(3).and.k<=hi(3)) cycle
        do i=lo(1),hi(1)
          do c=1,3
            fl=lo
            fl(c)=fl(c)-1
            do side=0,1
              idx=[i,j,k]
              idx(c)=idx(c)-side
              if(all(idx>=fl).and.all(idx<=hi)) cycle
              mask(idx(1),idx(2),idx(3),c)=.true.
            end do
          end do
        end do
      end do
    end do
  end subroutine

  ! Planner callbacks run collectively outside block loops. The cached sparse
  ! directory lists only leaves; virtual staging never allocates a finest mesh.
  subroutine yy_ct_prepare_directory(level)
    integer,optional,intent(in) :: level
    integer :: n,patch,leaf_level,lo(0:1),hi(0:1)
    if(.not.present(level)) call mpistop('adaptive CT plan requires a receiver level')
    if(ct_directory_epoch/=particle_mesh_epoch) then
      call build_directory(ct_directory)
      lo=huge(1)
      hi=0
      do n=1,size(ct_directory%leaf,2)
        patch=ct_directory%leaf(key_rank,n)/(yy_world_size/2)
        leaf_level=ct_directory%leaf(key_level,n)
        lo(patch)=min(lo(patch),leaf_level)
        hi(patch)=max(hi(patch),leaf_level)
      end do
      ct_uniform_level=merge(lo,0,lo==hi)
      ct_directory_epoch=particle_mesh_epoch
    end if
    ct_build_level=level
  end subroutine

  subroutine yy_ct_patch_owner(coords,idx,c,edge,rank,grid)
    integer,intent(inout) :: coords(3)
    integer,intent(in) :: idx(3),c
    logical,intent(in) :: edge
    integer,intent(out) :: rank,grid
    integer :: leaf,bn(3)
    if(edge) then
      bn=[block_nx1,block_nx2,block_nx3]
      call ct_same_level_entity(ct_build_level,yy_patch,bn,idx,c,.true.,coords,leaf)
    else
      leaf=find_leaf(ct_directory,ct_build_level,coords+1,yy_patch)
    end if
    rank=-1
    grid=-1
    if(leaf==0) return
    rank=ct_directory%leaf(key_rank,leaf)
    grid=ct_directory%leaf(key_grid,leaf)
  end subroutine

  ! A native edge can have up to four same-level copies. Pick the lower
  ! geometric owner among the actual leaves even at a mixed-level junction.
  subroutine ct_same_level_entity(level,patch,bn,idx,c,edge,coords,leaf)
    integer,intent(in) :: level,patch,bn(3),idx(3),c
    logical,intent(in) :: edge
    integer,intent(out) :: coords(3),leaf
    integer :: upper(3),lower(3),i,j,k
    logical :: shared_axis(3)
    lower=max(0,idx/bn)
    upper=lower
    shared_axis=edge
    shared_axis(c)=.not.edge
    where(shared_axis.and.modulo(idx+1,bn)==0.and.idx>=0) upper=upper+1
    upper=min(upper,[ng1(level),ng2(level),ng3(level)]-1)
    leaf=0
    do k=lower(3),upper(3)
      do j=lower(2),upper(2)
        do i=lower(1),upper(1)
          coords=[i,j,k]
          leaf=find_leaf(ct_directory,level,coords+1,patch)
          if(leaf/=0) return
        end do
      end do
    end do
    coords=lower
  end subroutine

  subroutine yy_ct_donor_sources(family,c,idx,ranks,requests,weights,levels,n)
    integer,intent(in) :: family,c,idx(3)
    integer,intent(out) :: ranks(:),requests(:,:),levels(:),n
    double precision,intent(out) :: weights(:)
    call ct_donor_sources(ct_build_level,1-yy_patch,family,c,idx,ranks,requests,weights,levels,n)
  end subroutine

  ! Inspect a coarse parent’s donor-aligned prolongation, including shared faces.
  ! The production plan must have prepared the current sparse leaf directory.
  subroutine yy_ct_virtual_face_sources(level,patch,c,idx,ranks,requests,weights,n)
    integer,intent(in) :: level,patch,c,idx(3)
    integer,intent(out) :: ranks(:),requests(:,:),n
    double precision,intent(out) :: weights(:)
    integer :: levels(size(ranks))
    if(ct_directory_epoch/=particle_mesh_epoch) call mpistop('adaptive CT leaf directory is stale')
    if(level<2.or.patch<0.or.patch>1.or.c<1.or.c>3) call mpistop('invalid virtual CT child face')
    call ct_donor_sources(level,patch,yy_face_location,c,idx,ranks,requests,weights,levels,n,source_level=level-1)
  end subroutine

  subroutine ct_donor_sources(l,patch,family,c,idx,ranks,requests,weights,levels,n,source_level)
    integer,intent(in) :: l,patch,family,c,idx(3)
    integer,intent(out) :: ranks(:),requests(:,:),levels(:),n
    double precision,intent(out) :: weights(:)
    integer,optional,intent(in) :: source_level
    integer :: ld,leaf,d,term(4,3),count(3),site(3)
    integer :: bn(3),dims(3)
    double precision :: coef(4,3),measure,weight
    logical :: integrated(3),ok,wider
    bn=[block_nx1,block_nx2,block_nx3]
    if(present(source_level)) then
      ld=source_level
    else
      leaf=covering_leaf(ct_directory,l,int(max(idx,0),int64),patch)
      ld=ct_directory%leaf(key_level,leaf)
    end if
    n=0
    if(abs(ld-l)>1) return
    if(ld==l) then
      weight=1.d0
      if(family==yy_edge_location) weight=1.d0/ct_measure(l,family,c,idx)
      call append(ld,idx,weight,ok)
      return
    end if
    integrated=(family==yy_face_location)
    integrated(c)=.not.integrated(c)
    measure=ct_measure(l,family,c,idx)
    if(ld>l) then
      ! Restrict actual fine subface fluxes or fine segment line integrals.
      do d=1,3
        if(integrated(d)) then
          count(d)=2
          term(1:2,d)=2*idx(d)+[0,1]
        else
          count(d)=1
          term(1,d)=2*idx(d)+1
        end if
        coef(:,d)=1.d0
      end do
      call tensor_terms(ld,ok)
      if(.not.ok) n=0
    else
      dims=int(domain_cells(ld))
      ! The integral and point maps commute with the incidence difference.
      ! Use one tensor pair throughout the donor patch. A face-dependent
      ! fallback would break flux closure between faces of the same child.
      wider=ct_uniform_level(patch)==ld
      do d=1,3
        call ct_prolong_axis(idx(d),dims(d),integrated(d),wider,term(:,d),coef(:,d),count(d))
      end do
      call tensor_terms(ld,ok)
      if(.not.ok) n=0
    end if
  contains
    subroutine tensor_terms(level,complete)
      integer,intent(in) :: level
      logical,intent(out) :: complete
      integer :: t1,t2,t3
      complete=.false.
      do t3=1,count(3)
        do t2=1,count(2)
          do t1=1,count(1)
            site=[term(t1,1),term(t2,2),term(t3,3)]
            weight=coef(t1,1)*coef(t2,2)*coef(t3,3)/measure
            if(family==yy_face_location) weight=weight*ct_measure(level,family,c,site)
            call append(level,site,weight,complete)
            if(.not.complete) return
          end do
        end do
      end do
    end subroutine

    subroutine append(level,point,factor,found)
      integer,intent(in) :: level,point(3)
      double precision,intent(in) :: factor
      logical,intent(out) :: found
      integer :: g,coords(3),segment,children(3,2),child_coords(3,2),child_leaf(2)
      coords=max(0,point/bn)
      if(family==yy_edge_location) then
        call ct_same_level_entity(level,patch,bn,point,c,.true.,coords,g)
        if(level<ct_directory%deepest) then
          ! Native edge reflux uses fine collinear segments where they exist.
          ! Use the same authority for coarse donors before cross-patch reads.
          do segment=1,2
            children(:,segment)=2*point+1
            children(c,segment)=2*point(c)+segment-1
            call ct_same_level_entity(level+1,patch,bn,children(:,segment),c,.true., &
              child_coords(:,segment),child_leaf(segment))
          end do
          if(all(child_leaf/=0)) then
            found=.true.
            do segment=1,2
              call record(level+1,children(:,segment),child_coords(:,segment),child_leaf(segment),factor)
            end do
            return
          end if
        end if
      else
        call ct_same_level_entity(level,patch,bn,point,c,.false.,coords,g)
      end if
      found=g/=0
      if(found) call record(level,point,coords,g,factor)
    end subroutine

    subroutine record(level,point,coords,g,factor)
      integer,intent(in) :: level,point(3),coords(3),g
      double precision,intent(in) :: factor
      if(abs(level-l)>1) call mpistop('adaptive CT actual donor exceeds the one-level stencil bound')
      n=n+1
      if(n>size(ranks)) call mpistop('adaptive CT donor stencil exceeds storage capacity')
      ranks(n)=ct_directory%leaf(key_rank,g)+1
      requests(:,n)=[ct_directory%leaf(key_grid,g),point-coords*bn,c]
      levels(n)=level
      weights(n)=factor
    end subroutine
  end subroutine

  pure subroutine ct_prolong_axis(index,ncoarse,integrated,wider,sites,coeff,count)
    integer,intent(in) :: index,ncoarse
    logical,intent(in) :: integrated,wider
    integer,intent(out) :: sites(4),count
    double precision,intent(out) :: coeff(4)
    integer :: parent,vertex,odd
    if(integrated) then
      parent=index/2
      count=1
      sites(1)=parent
      coeff(1)=0.5d0
      if(.not.wider) return
      odd=modulo(index,2)
      if(parent==0) then
        count=2
        sites(1:2)=[parent,parent+1]
        coeff(1:2)=[0.5d0+(1-2*odd)*0.125d0,-(1-2*odd)*0.125d0]
      else if(parent==ncoarse-1) then
        count=2
        sites(1:2)=[parent,parent-1]
        coeff(1:2)=[0.5d0-(1-2*odd)*0.125d0,(1-2*odd)*0.125d0]
      else
        count=3
        sites(2:3)=[parent-1,parent+1]
        coeff(2:3)=[(1-2*odd)*0.0625d0,-(1-2*odd)*0.0625d0]
      end if
    else
      vertex=index+1
      parent=vertex/2-1
      count=1
      sites(1)=parent
      coeff(1)=1.d0
      if(modulo(vertex,2)==0) return
      count=2
      sites(2)=parent+1
      coeff(1:2)=0.5d0
      if(.not.wider) return
      if(parent==-1) then
        count=3
        sites(3)=1
        coeff(1:3)=[0.375d0,0.75d0,-0.125d0]
      else if(parent==ncoarse-2) then
        count=3
        sites(3)=parent-1
        coeff(1:3)=[0.75d0,0.375d0,-0.125d0]
      else
        count=4
        sites(3:4)=[parent-1,parent+2]
        coeff(1:4)=[0.5625d0,0.5625d0,-0.0625d0,-0.0625d0]
      end if
    end if
  end subroutine

  ! Native spherical metric for the positive face normal or edge tangent.
  pure double precision function ct_measure(level,family,c,idx) result(value)
    integer,intent(in) :: level,family,c,idx(3)
    double precision :: r0,r1,t1,dteta,dphi
    r1=radial_edge(level,int(idx(1),int64)+1_int64)
    dteta=dx(2,level)
    dphi=dx(3,level)
    t1=xprobmin2+(idx(2)+1)*dteta
    if(family==yy_edge_location) then
      select case(c)
      case(1);value=r1-radial_edge(level,int(idx(1),int64))
      case(2);value=r1*dteta
      case(3);value=r1*sin(t1)*dphi
      end select
    else
      select case(c)
      case(1)
        value=2.d0*r1*r1*sin(t1-0.5d0*dteta)*sin(0.5d0*dteta)*dphi
      case(2)
        r0=radial_edge(level,int(idx(1),int64))
        value=0.5d0*(r1-r0)*(r1+r0)*sin(t1)*dphi
      case(3)
        r0=radial_edge(level,int(idx(1),int64))
        value=0.5d0*(r1-r0)*(r1+r0)*dteta
      end select
    end if
  end function
  }

  ! Receiver rows are local block ID, native i,j,k and direction. Select only
  ! desired interface edges; uncomputed UCT workspace may never be a donor.
  subroutine yy_prepare_ct_edges(plan,s,receivers,radial_extrapolate,native_midpoint,level)
    type(yy_ct_stage_plan),intent(out) :: plan
    type(state),intent(in) :: s(max_blocks)
    integer,intent(in) :: receivers(:,:)
    logical,optional,intent(in) :: radial_extrapolate,native_midpoint
    integer,optional,intent(in) :: level
    {^IFTHREED
    logical :: midpoint
    integer :: nn(3),bn(3),n,g,c,idx(3),lo(3),hi(3)
    integer,allocatable :: directory(:,:)
    double precision,allocatable :: radial(:),centers(:,:),widths(:)
    double precision :: angular(2,2),offset(3)
    if(size(receivers,1)/=5) call mpistop('CT edge receiver shape must be (5,n)')
    call yy_staggered_mesh(s,bn,nn,radial,angular,directory,level)
    allocate(centers(3,size(receivers,2)),widths(size(receivers,2)))
    do n=1,size(receivers,2)
      g=receivers(1,n);c=receivers(5,n);idx=receivers(2:4,n)
      if(c<1.or.c>3) call mpistop('invalid receiving edge direction')
      if(.not.any(igrids(1:igridstail)==g)) call mpistop('receiving CT block is not local')
      lo=[ixMlo1-1,ixMlo2-1,ixMlo3-1];hi=[ixMhi1,ixMhi2,ixMhi3];lo(c)=lo(c)+1
      if(any(idx<lo).or.any(idx>hi)) call mpistop('receiver outside valid UCT stencil')
      offset=0.5d0;offset(c)=0.d0
      centers(:,n)=s(g)%x(idx(1),idx(2),idx(3),:)+offset*s(g)%dx(idx(1),idx(2),idx(3),:)
      widths(n)=s(g)%dx(idx(1),idx(2),idx(3),c)
    end do
    midpoint=.true.;if(present(native_midpoint))midpoint=native_midpoint
    if(refine_max_level>1) then
      call yy_ct_prepare_directory(level)
      call yy_ct_edge_build(plan%edge,MPI_COMM_WORLD,yy_patch,bn,nn,radial,angular,directory, &
        centers,widths,receivers(5,:),radial_extrapolate,native_midpoint=midpoint,sources=yy_ct_donor_sources)
      if(midpoint) then
        plan%donor_level_range=plan%edge%samples%donor_level_range
      else
        allocate(plan%donor_level_range(2,size(receivers,2)))
        do n=1,size(receivers,2)
          plan%donor_level_range(:,n)=[minval(plan%edge%samples%donor_level_range(1,2*n-1:2*n)), &
            maxval(plan%edge%samples%donor_level_range(2,2*n-1:2*n))]
        end do
      end if
    else
      call yy_ct_edge_build(plan%edge,MPI_COMM_WORLD,yy_patch,bn,nn,radial,angular,directory, &
        centers,widths,receivers(5,:),radial_extrapolate,native_midpoint=midpoint)
      allocate(plan%donor_level_range(2,size(receivers,2)))
      plan%donor_level_range=1
    end if
    if(any(.not.plan%edge%covered)) call mpistop('CT edge lacks complete donor stencil')
    do n=1,size(plan%edge%samples%request,2)
      idx=plan%edge%samples%request(2:4,n)+[ixMlo1,ixMlo2,ixMlo3]
      c=plan%edge%samples%request(5,n)
      lo=[ixMlo1-1,ixMlo2-1,ixMlo3-1];lo(c)=lo(c)+1
      if(any(idx<lo).or.any(idx>[ixMhi1,ixMhi2,ixMhi3])) &
        call mpistop('donor outside computed UCT edge stencil')
    end do
    plan%receivers=receivers
    allocate(plan%packed(size(plan%edge%samples%request,2)),plan%received(size(receivers,2)))
    plan%mesh_epoch=particle_mesh_epoch
    }
    {^IFONED call mpistop('CT edge exchange requires 3D') }
    {^IFTWOD call mpistop('CT edge exchange requires 3D') }
  end subroutine

  ! Invoke at the all-block RK-stage boundary, before any edge is overwritten.
  ! All requested donor integrals are packed before the collective returns.
  subroutine yy_exchange_ct_edges(plan,fluxes,defer_write)
    type(yy_ct_stage_plan),intent(inout) :: plan
    type(ct_stage_flux),intent(inout) :: fluxes(:)
    logical,optional,intent(in) :: defer_write
    {^IFTHREED
    integer :: n,g,c,i,j,k
    if(.not.allocated(plan%packed)) call mpistop('CT edge plan has not been built')
    if(size(fluxes)<max_blocks) call mpistop('CT stage block array is too small')
    if(plan%mesh_epoch/=particle_mesh_epoch) call mpistop('CT edge plan is stale')
    do n=1,size(plan%packed)
      g=plan%edge%samples%request(1,n);c=plan%edge%samples%request(5,n)
      i=plan%edge%samples%request(2,n)+ixMlo1
      j=plan%edge%samples%request(3,n)+ixMlo2
      k=plan%edge%samples%request(4,n)+ixMlo3
      if(.not.allocated(fluxes(g)%edge)) call mpistop('CT donor stage workspace missing')
      plan%packed(n)=fluxes(g)%edge(i,j,k,c)
    end do
    call yy_ct_edge_exchange(plan%edge,plan%packed,plan%received)
    if(present(defer_write)) then
      if(defer_write) return
    end if
    call yy_write_ct_edges(plan,fluxes)
    }
    {^IFONED call mpistop('CT edge exchange requires 3D') }
    {^IFTWOD call mpistop('CT edge exchange requires 3D') }
  end subroutine

  subroutine yy_write_ct_edges(plan,fluxes)
    type(yy_ct_stage_plan),intent(in) :: plan
    type(ct_stage_flux),intent(inout) :: fluxes(:)
    {^IFTHREED
    integer :: n,g,c,i,j,k
    do n=1,size(plan%received)
      if(.not.plan%edge%covered(n)) cycle
      g=plan%receivers(1,n);i=plan%receivers(2,n);j=plan%receivers(3,n)
      k=plan%receivers(4,n);c=plan%receivers(5,n)
      if(.not.allocated(fluxes(g)%edge)) call mpistop('CT receiving workspace missing')
      fluxes(g)%edge(i,j,k,c)=plan%received(n)
    end do
    }
    {^IFONED call mpistop('CT edge exchange requires 3D') }
    {^IFTWOD call mpistop('CT edge exchange requires 3D') }
  end subroutine

  ! Enumerate the union of faces bounding external angular ghost cells. Active
  ! faces, including active-boundary anchors, are excluded before interpolation.
  ! Build the same projection used for ghost faces on explicit receiver faces.
  ! Diagnostic callers may select active faces in an independent scratch state.
  subroutine yy_prepare_ct_faces(plan,s,receivers,level,face_moments)
    type(yy_ct_ghost_plan),intent(out) :: plan
    type(state),intent(in) :: s(max_blocks)
    integer,intent(in) :: receivers(:,:)
    integer,optional,intent(in) :: level
    logical,optional,intent(in) :: face_moments
    {^IFTHREED
    integer :: bn(3),nn(3),n,g,c,idx(3),above(3)
    integer,allocatable :: directory(:,:)
    logical :: moments
    double precision,allocatable :: radial(:),centers(:,:),widths(:,:)
    double precision :: angular(2,2),offset(3)
    if(size(receivers,1)/=5) call mpistop('CT face receiver shape must be (5,n)')
    call yy_staggered_mesh(s,bn,nn,radial,angular,directory,level)
    allocate(centers(3,size(receivers,2)),widths(3,size(receivers,2)))
    allocate(plan%values(size(receivers,2)),plan%donor_level_range(2,size(receivers,2)))
    plan%receivers=receivers
    do n=1,size(receivers,2)
      g=receivers(1,n);idx=receivers(2:4,n);c=receivers(5,n)
      if(c<1.or.c>3) call mpistop('invalid receiving face direction')
      if(.not.any(igrids(1:igridstail)==g)) call mpistop('receiving CT block is not local')
      above=max(idx,[ixGlo1,ixGlo2,ixGlo3])
      offset=0.d0;offset(c)=0.5d0
      if(idx(c)<above(c)) offset(c)=-0.5d0
      widths(:,n)=s(g)%dx(above(1),above(2),above(3),:)
      centers(:,n)=s(g)%x(above(1),above(2),above(3),:)+offset*widths(:,n)
    end do
    moments=.false.
    if(present(face_moments)) moments=face_moments
    if(refine_max_level>1) then
      call yy_ct_prepare_directory(level)
      if(moments) then
        call yy_ct_face_build(plan%faces,MPI_COMM_WORLD,yy_patch,bn,nn,radial,angular,directory, &
          centers,receivers(5,:),widths,sources=yy_ct_donor_sources)
      else
        call yy_ct_face_build(plan%faces,MPI_COMM_WORLD,yy_patch,bn,nn,radial,angular,directory, &
          centers,receivers(5,:),sources=yy_ct_donor_sources)
      end if
      plan%donor_level_range(:,plan%faces%rslot)=plan%faces%radial%donor_level_range
      plan%donor_level_range(:,plan%faces%tslot)=plan%faces%tangent%donor_level_range
    else
      if(moments) then
        call yy_ct_face_build(plan%faces,MPI_COMM_WORLD,yy_patch,bn,nn,radial,angular,directory, &
          centers,receivers(5,:),widths)
      else
        call yy_ct_face_build(plan%faces,MPI_COMM_WORLD,yy_patch,bn,nn,radial,angular,directory,centers,receivers(5,:))
      end if
      plan%donor_level_range=1
    end if
    if(any(.not.plan%faces%covered)) call mpistop('CT ghost face lacks complete donor stencil')
    allocate(plan%radial(size(plan%faces%radial%request,2)),plan%tangent(size(plan%faces%tangent%request,2)))
    plan%mesh_epoch=particle_mesh_epoch
    }
    {^IFONED call mpistop('magnetic face exchange requires 3D') }
    {^IFTWOD call mpistop('magnetic face exchange requires 3D') }
  end subroutine

  subroutine yy_prepare_ct_ghosts(plan,s,face_moments,level)
    type(yy_ct_ghost_plan),intent(out) :: plan
    type(state),intent(in) :: s(max_blocks)
    logical,optional,intent(in) :: face_moments
    integer,optional,intent(in) :: level
    {^IFTHREED
    integer :: n,g,i,j,k,c,side,pass,t,idx(3),lo(3),hi(3),fl(3),pig(3),nb(3),l
    integer :: low(3),high(3),bn(3)
    logical :: sides(2,2),selected,mask(ixGlo1-1:ixGhi1,ixGlo2-1:ixGhi2,ixGlo3-1:ixGhi3,3)
    integer,allocatable :: receivers(:,:)
    l=1
    if(present(level)) l=level
    lo=[ixMlo1,ixMlo2,ixMlo3];hi=[ixMhi1,ixMhi2,ixMhi3];nb=[ng1(l),ng2(l),ng3(l)]
    bn=[block_nx1,block_nx2,block_nx3]
    do pass=1,2
      t=0
      do n=1,igridstail
        g=igrids(n)
        if(node(plevel_,g)/=l) cycle
        pig=[node(pig1_,g),node(pig2_,g),node(pig3_,g)]
        if(refine_max_level>1) then
          call ct_face_bounds(l,(pig-1)*bn,bn,low,high,sides)
          call ct_face_footprint(low,high,sides,mask)
        else
          sides(1,:)=(pig(2:3)==1);sides(2,:)=(pig(2:3)==nb(2:3));mask=.false.
          do k=ixGlo3,ixGhi3
            do j=ixGlo2,ixGhi2
              selected=((j<lo(2).and.sides(1,1)).or.(j>hi(2).and.sides(2,1))) &
                .and.k>=lo(3).and.k<=hi(3)
              if((k<lo(3).and.sides(1,2)).or.(k>hi(3).and.sides(2,2))) &
                selected=(j>=lo(2).and.j<=hi(2)).or.(j<lo(2).and.sides(1,1)).or.(j>hi(2).and.sides(2,1))
              if(.not.selected) cycle
              do i=lo(1),hi(1)
                do c=1,3
                  do side=0,1
                    idx=[i,j,k];idx(c)=idx(c)-side
                    fl=lo;fl(c)=fl(c)-1
                    if(all(idx>=fl).and.all(idx<=hi)) cycle
                    mask(idx(1),idx(2),idx(3),c)=.true.
                  end do
                end do
              end do
            end do
          end do
        end if
        do c=1,3
          do k=ixGlo3-1,ixGhi3
            do j=ixGlo2-1,ixGhi2
              do i=ixGlo1-1,ixGhi1
                if(.not.mask(i,j,k,c)) cycle
                t=t+1
                if(pass==1) cycle
                receivers(:,t)=[g,i,j,k,c]
              end do
            end do
          end do
        end do
      end do
      if(pass==1) allocate(receivers(5,t))
    end do
    call yy_prepare_ct_faces(plan,s,receivers,level,face_moments)
    }
    {^IFONED call mpistop('magnetic ghost fill requires 3D') }
    {^IFTWOD call mpistop('magnetic ghost fill requires 3D') }
  end subroutine

  ! Fill face B in angular halos, optionally followed by outward flux closure.
  ! Both donor families are packed before any writes, so a ghost never becomes
  ! a donor through ordering. Caller handles physical radial and same-patch BCs.
  subroutine yy_exchange_ct_faces(plan,s)
    type(yy_ct_ghost_plan),intent(inout) :: plan
    type(state),intent(inout) :: s(max_blocks)
    {^IFTHREED
    integer :: n,g,i,j,k,c,idx(3),lo(3)
    if(.not.allocated(plan%values)) call mpistop('magnetic face plan has not been built')
    if(plan%mesh_epoch/=particle_mesh_epoch) call mpistop('magnetic face plan is stale')
    lo=[ixMlo1,ixMlo2,ixMlo3]
    do n=1,size(plan%radial)
      g=plan%faces%radial%request(1,n);idx=plan%faces%radial%request(2:4,n)+lo
      c=plan%faces%radial%request(5,n)
      plan%radial(n)=s(g)%ws(idx(1),idx(2),idx(3),c)
    end do
    do n=1,size(plan%tangent)
      g=plan%faces%tangent%request(1,n);idx=plan%faces%tangent%request(2:4,n)+lo
      c=plan%faces%tangent%request(5,n)
      plan%tangent(n)=s(g)%ws(idx(1),idx(2),idx(3),c)
    end do
    call yy_ct_face_exchange(plan%faces,plan%radial,plan%tangent,plan%values)
    do n=1,size(plan%values)
      g=plan%receivers(1,n);i=plan%receivers(2,n);j=plan%receivers(3,n)
      k=plan%receivers(4,n);c=plan%receivers(5,n)
      s(g)%ws(i,j,k,c)=plan%values(n)
    end do
    }
    {^IFONED call mpistop('magnetic face exchange requires 3D') }
    {^IFTWOD call mpistop('magnetic face exchange requires 3D') }
  end subroutine

  subroutine yy_fill_ct_ghosts(plan,s,constrain,level)
    type(yy_ct_ghost_plan),intent(inout) :: plan
    type(state),intent(inout) :: s(max_blocks)
    logical,intent(in) :: constrain
    integer,optional,intent(in) :: level
    {^IFTHREED
    integer :: n,g,lo(3),hi(3),status,pig(3),nb(3),l,low(3),high(3),bn(3)
    logical :: sides(2,2)
    call yy_exchange_ct_faces(plan,s)
    l=1
    if(present(level)) l=level
    lo=[ixMlo1,ixMlo2,ixMlo3];hi=[ixMhi1,ixMhi2,ixMhi3];nb=[ng1(l),ng2(l),ng3(l)]
    if(constrain) then
      do n=1,igridstail
        g=igrids(n)
        if(node(plevel_,g)/=l) cycle
        pig=[node(pig1_,g),node(pig2_,g),node(pig3_,g)]
        low=lo
        high=hi
        if(refine_max_level>1) then
          bn=[block_nx1,block_nx2,block_nx3]
          call ct_face_bounds(l,(pig-1)*bn,bn,low,high,sides)
          low=low+lo
          high=high+lo
        else
          sides(1,:)=(pig(2:3)==1);sides(2,:)=(pig(2:3)==nb(2:3))
        end if
        call yy_ct_close_ghosts([ixGlo1-1,ixGlo2-1,ixGlo3-1],low,high, &
          [ixGlo1,ixGlo2,ixGlo3],[ixGhi1,ixGhi2,ixGhi3],sides,s(g)%ws,s(g)%surfaceC,status)
        if(status/=0) call mpistop('invalid geometry or field in magnetic ghost closure')
      end do
    end if
    }
    {^IFONED call mpistop('magnetic ghost fill requires 3D') }
    {^IFTWOD call mpistop('magnetic ghost fill requires 3D') }
  end subroutine

  subroutine yy_prepare_ct_patch(plan,s,level)
    type(yy_ct_patch_plan),intent(out) :: plan
    type(state),intent(in) :: s(max_blocks)
    integer,optional,intent(in) :: level
    integer :: bn(3),nn(3)
    integer,allocatable :: directory(:,:)
    double precision,allocatable :: radial(:)
    double precision :: angular(2,2)
    call yy_staggered_mesh(s,bn,nn,radial,angular,directory,level)
    {^IFTHREED
    if(refine_max_level>1) then
      call yy_ct_prepare_directory(level)
      call yy_patch_copy_build(plan%edge,MPI_COMM_WORLD,yy_patch,yy_sync_edges,bn,nn,nghostcells,directory, &
        owner_lookup=yy_ct_patch_owner)
      call yy_patch_copy_build(plan%interior_faces,MPI_COMM_WORLD,yy_patch,yy_sync_faces,bn,nn,nghostcells,directory, &
        owner_lookup=yy_ct_patch_owner)
      call yy_patch_copy_build(plan%extended_faces,MPI_COMM_WORLD,yy_patch,yy_sync_extended_faces, &
        bn,nn,nghostcells,directory,owner_lookup=yy_ct_patch_owner)
    else
      call yy_patch_copy_build(plan%edge,MPI_COMM_WORLD,yy_patch,yy_sync_edges,bn,nn,nghostcells,directory)
      call yy_patch_copy_build(plan%interior_faces,MPI_COMM_WORLD,yy_patch,yy_sync_faces,bn,nn,nghostcells,directory)
      call yy_patch_copy_build(plan%extended_faces,MPI_COMM_WORLD,yy_patch,yy_sync_extended_faces, &
        bn,nn,nghostcells,directory)
    end if
    }
    plan%mesh_epoch=particle_mesh_epoch
  end subroutine

  logical function yy_ct_edge_is_owner(g,idx,c) result(owner)
    integer,intent(in) :: g,idx(3),c
    {^IFTHREED
    integer :: d,lo(3),pig(3),bn(3),coords(3),leaf
    lo=[ixMlo1,ixMlo2,ixMlo3];pig=[node(pig1_,g),node(pig2_,g),node(pig3_,g)]
    if(ct_directory_epoch==particle_mesh_epoch) then
      bn=[block_nx1,block_nx2,block_nx3]
      call ct_same_level_entity(node(plevel_,g),yy_patch,bn,(pig-1)*bn+idx-lo,c,.true.,coords,leaf)
      owner=all(coords==pig-1)
      return
    end if
    owner=.true.
    do d=1,3
      if(d==c) cycle
      if(idx(d)==lo(d)-1.and.pig(d)>1) owner=.false.
    end do
    }
    {^IFONED owner=.false. }
    {^IFTWOD owner=.false. }
  end function

  subroutine yy_sync_ct_edges(plan,fluxes)
    type(yy_ct_patch_plan),intent(inout) :: plan
    type(ct_stage_flux),intent(inout) :: fluxes(:)
    {^IFTHREED
    integer :: n,g,c,idx(3),lo(3)
    if(plan%mesh_epoch/=particle_mesh_epoch) call mpistop('CT patch plan is stale')
    if(.not.allocated(plan%edge%send)) call mpistop('CT patch plan has not been built')
    if(size(fluxes)<max_blocks) call mpistop('CT stage array is too small')
    lo=[ixMlo1,ixMlo2,ixMlo3]
    do n=1,size(plan%edge%send)
      g=plan%edge%request(1,n);idx=plan%edge%request(2:4,n)+lo;c=plan%edge%request(5,n)
      if(.not.allocated(fluxes(g)%edge)) call mpistop('CT patch donor workspace missing')
      plan%edge%send(n)=fluxes(g)%edge(idx(1),idx(2),idx(3),c)
    end do
    call yy_patch_copy_exchange(plan%edge)
    do n=1,size(plan%edge%recv)
      g=plan%edge%receivers(1,n);idx=plan%edge%receivers(2:4,n)+lo;c=plan%edge%receivers(5,n)
      if(.not.allocated(fluxes(g)%edge)) call mpistop('CT patch receiver workspace missing')
      fluxes(g)%edge(idx(1),idx(2),idx(3),c)=plan%edge%recv(n)
    end do
    }
    {^IFONED call mpistop('CT patch sync requires 3D') }
    {^IFTWOD call mpistop('CT patch sync requires 3D') }
  end subroutine

  ! The overset receiver list must select geometric owners. Otherwise the final
  ! same-patch broadcast would erase an update applied only to a non-owner copy.
  subroutine yy_ordered_ct_edges(patch_plan,cross_plan,fluxes)
    type(yy_ct_patch_plan),intent(inout) :: patch_plan
    type(yy_ct_stage_plan),intent(inout) :: cross_plan
    type(ct_stage_flux),intent(inout) :: fluxes(:)
    integer :: n
    if(.not.allocated(cross_plan%receivers)) call mpistop('CT cross-patch plan has not been built')
    do n=1,size(cross_plan%receivers,2)
      if(.not.yy_ct_edge_is_owner(cross_plan%receivers(1,n),cross_plan%receivers(2:4,n), &
        cross_plan%receivers(5,n))) call mpistop('CT overset receiver must be the geometric owner')
    end do
    call yy_sync_ct_edges(patch_plan,fluxes)
    call yy_exchange_ct_edges(cross_plan,fluxes)
    call yy_sync_ct_edges(patch_plan,fluxes)
  end subroutine

  subroutine yy_sync_ct_faces(plan,s,extended)
    type(yy_ct_patch_plan),target,intent(inout) :: plan
    type(state),intent(inout) :: s(max_blocks)
    logical,intent(in) :: extended
    {^IFTHREED
    type(yy_patch_copy_plan),pointer :: table
    integer :: n,g,c,idx(3),lo(3)
    if(plan%mesh_epoch/=particle_mesh_epoch) call mpistop('CT patch plan is stale')
    if(extended) then
      table=>plan%extended_faces
    else
      table=>plan%interior_faces
    end if
    if(.not.allocated(table%send)) call mpistop('CT face sync plan has not been built')
    lo=[ixMlo1,ixMlo2,ixMlo3]
    do n=1,size(table%send)
      g=table%request(1,n);idx=table%request(2:4,n)+lo;c=table%request(5,n)
      if(.not.allocated(s(g)%ws)) call mpistop('CT face donor state missing')
      table%send(n)=s(g)%ws(idx(1),idx(2),idx(3),c)
    end do
    call yy_patch_copy_exchange(table)
    do n=1,size(table%recv)
      g=table%receivers(1,n);idx=table%receivers(2:4,n)+lo;c=table%receivers(5,n)
      if(.not.allocated(s(g)%ws)) call mpistop('CT face receiver state missing')
      s(g)%ws(idx(1),idx(2),idx(3),c)=table%recv(n)
    end do
    }
    {^IFONED call mpistop('CT face sync requires 3D') }
    {^IFTWOD call mpistop('CT face sync requires 3D') }
  end subroutine

  subroutine yy_ordered_ct_ghosts(patch_plan,cross_plan,s)
    type(yy_ct_patch_plan),intent(inout) :: patch_plan
    type(yy_ct_ghost_plan),intent(inout) :: cross_plan
    type(state),intent(inout) :: s(max_blocks)
    call yy_sync_ct_faces(patch_plan,s,.false.)
    call yy_fill_ct_ghosts(cross_plan,s,.true.)
    call yy_sync_ct_faces(patch_plan,s,.true.)
  end subroutine

  ! First finish angular faces, then radial physical halos across the entire
  ! angular extent, including corners. Cell thermodynamic halos must already be
  ! populated (ordinary/overset exchange) before this explicit routine is called.
  subroutine yy_complete_ct_faces(patch_plan,cross_plan,s)
    type(yy_ct_patch_plan),intent(inout) :: patch_plan
    type(yy_ct_ghost_plan),intent(inout) :: cross_plan
    type(state),intent(inout) :: s(max_blocks)
    call yy_ordered_ct_ghosts(patch_plan,cross_plan,s)
    call yy_ct_radial_faces(s)
  end subroutine

  ! Explicit radial face policy: constant tangential B,
  ! zero finite-volume divergence in radial ghost cells. Boundary Br is anchored.
  ! This is not a characteristic outflow or a line-tied velocity boundary.
  subroutine yy_ct_radial_faces(s)
    type(state),intent(inout) :: s(max_blocks)
    {^IFTHREED
    integer :: n,g,side,anchor,first,last,stride,i,j,k,c,jlo,klo
    double precision :: transverse
    if(.not.yinyang_enabled.or.lb_automatic) call mpistop('invalid Yin-Yang CT radial layout')
    do n=1,igridstail
      g=igrids(n)
      if(.not.allocated(s(g)%ws)) call mpistop('CT radial face state missing')
      do side=1,2
        if(side==1.and.node(pig1_,g)/=1) cycle
        if(side==2.and.node(pig1_,g)/=ng1(node(plevel_,g))) cycle
        anchor=ixMlo1;first=ixMlo1-1;last=ixGlo1;stride=-1
        if(side==2) then
          anchor=ixMhi1;first=ixMhi1+1;last=ixGhi1;stride=1
        end if
        ! Angular face halos at the anchor radius are already complete.
        do c=2,3
          jlo=ixGlo2;klo=ixGlo3
          if(c==2) jlo=jlo-1
          if(c==3) klo=klo-1
          do k=klo,ixGhi3
            do j=jlo,ixGhi2
              do i=first,last,stride
                s(g)%ws(i,j,k,c)=s(g)%ws(anchor,j,k,c)
              end do
            end do
          end do
        end do
        do k=ixGlo3,ixGhi3
          do j=ixGlo2,ixGhi2
            do i=first,last,stride
              transverse=s(g)%surfaceC(i,j,k,2)*s(g)%ws(i,j,k,2) &
                -s(g)%surfaceC(i,j-1,k,2)*s(g)%ws(i,j-1,k,2) &
                +s(g)%surfaceC(i,j,k,3)*s(g)%ws(i,j,k,3) &
                -s(g)%surfaceC(i,j,k-1,3)*s(g)%ws(i,j,k-1,3)
              if(side==1) then
                s(g)%ws(i-1,j,k,1)=(s(g)%surfaceC(i,j,k,1)*s(g)%ws(i,j,k,1)+transverse) &
                  /s(g)%surfaceC(i-1,j,k,1)
              else
                s(g)%ws(i,j,k,1)=(s(g)%surfaceC(i-1,j,k,1)*s(g)%ws(i-1,j,k,1)-transverse) &
                  /s(g)%surfaceC(i,j,k,1)
              end if
            end do
          end do
        end do
      end do
    end do
    }
    {^IFONED call mpistop('CT radial fill requires 3D') }
    {^IFTWOD call mpistop('CT radial fill requires 3D') }
  end subroutine

  ! Caller supplies mhd_ct_sync_center. Every ghost cell is visited exactly once:
  ! x slabs, y slabs at active x, z slabs at active x/y. Angular slabs run first.
  ! Active total energy is
  ! unchanged; ghosts retain their prescribed thermal energy when B is replaced.
  subroutine yy_sync_ct_state(s,sync_center)
    type(state),intent(inout) :: s(max_blocks)
    interface
      subroutine sync_center(ixO^L,s,keep_pressure)
        use mod_global_parameters
        integer,intent(in) :: ixO^L
        type(state),intent(inout) :: s
        logical,intent(in) :: keep_pressure
      end subroutine
    end interface
    integer :: n,g,d,side,i,anchor,ixO^L,lo(ndim),hi(ndim),amin(ndim),amax(ndim),gmin(ndim),gmax(ndim)
    if(ndim/=3.or..not.yinyang_enabled.or.lb_automatic) &
      call mpistop('CT state sync requires 3D Yin-Yang coordinates')
    {^D& amin(^D)=ixMlo^D;amax(^D)=ixMhi^D;gmin(^D)=ixGlo^D;gmax(^D)=ixGhi^D;\}
    do n=1,igridstail
      g=igrids(n)
      call sync_center(ixM^LL,s(g),.false.)
      do d=ndim,1,-1
        do side=1,2
          lo=gmin;hi=gmax
          lo(1:d-1)=amin(1:d-1);hi(1:d-1)=amax(1:d-1)
          if(side==1) then
            hi(d)=amin(d)-1
          else
            lo(d)=amax(d)+1
          end if
          {^IFTHREED
          if(d==1) then
            anchor=0
            if(side==1.and.node(pig1_,g)==1.and..not.special_inner_ct) anchor=ixMlo1
            if(side==2.and.node(pig1_,g)==ng1(node(plevel_,g))) anchor=ixMhi1
            if(anchor/=0) then
              ! Copy after active and angular ghost B/energy are consistent.
              do i=lo(1),hi(1)
                s(g)%w(i,ixGlo2:ixGhi2,ixGlo3:ixGhi3,:)=s(g)%w(anchor,ixGlo2:ixGhi2,ixGlo3:ixGhi3,:)
              end do
            end if
          end if
          }
          ixOmin^D=lo(^D);ixOmax^D=hi(^D);
          call sync_center(ixO^L,s(g),.true.)
        end do
      end do
    end do
  end subroutine

  ! Cache all owner angular-boundary edges, including radial endpoints.
  subroutine yy_runtime_ct_prepare(s)
    type(state),intent(in) :: s(max_blocks)
    {^IFTHREED
    integer,allocatable :: receivers(:,:)
    integer :: pass,n,g,c,i,j,k,idx(3),gi(3),bc(3),lo(3),hi(3),bn(3),nn(3),count,l
    logical :: local_needs_sync
    integer :: local_count(refine_max_level),global_count(refine_max_level)
    if(runtime_ct_epoch==particle_mesh_epoch) return
    if(allocated(runtime_ct)) deallocate(runtime_ct)
    local_count=0
    do n=1,igridstail
      l=node(plevel_,igrids(n))
      local_count(l)=local_count(l)+1
    end do
    call MPI_ALLREDUCE(local_count,global_count,refine_max_level,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierrmpi)
    allocate(runtime_ct(refine_max_level))
    do l=1,refine_max_level
      if(global_count(l)==0) cycle
      call yy_prepare_ct_patch(runtime_ct(l)%patch,s,l)
      call yy_prepare_ct_ghosts(runtime_ct(l)%ghosts,s,face_moments=yy_ct_face_moments,level=l)
      bn=[block_nx1,block_nx2,block_nx3];nn=bn*[ng1(l),ng2(l),ng3(l)]
      hi=[ixMhi1,ixMhi2,ixMhi3]
      do pass=1,2
        count=0
        do n=1,igridstail
          g=igrids(n)
          if(node(plevel_,g)/=l) cycle
          bc=[node(pig1_,g)-1,node(pig2_,g)-1,node(pig3_,g)-1]
          do c=1,3
            lo=[ixMlo1,ixMlo2,ixMlo3]-1;lo(c)=lo(c)+1
            do k=lo(3),hi(3)
              do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                  idx=[i,j,k];gi=bc*bn+idx-[ixMlo1,ixMlo2,ixMlo3]
                  if(.not.ct_seam_edge(gi,c,nn)) cycle
                  ! The user prescribes inner tangential EMFs in each patch basis.
                  if(special_inner_ct.and.c/=1.and.gi(1)==-1) cycle
                  if(.not.yy_ct_edge_is_owner(g,idx,c)) cycle
                  count=count+1
                  if(pass==2) receivers(:,count)=[g,idx,c]
                end do
              end do
            end do
          end do
        end do
        if(pass==1) allocate(receivers(5,count))
      end do
      call yy_prepare_ct_edges(runtime_ct(l)%edges,s,receivers,radial_extrapolate=.true.,level=l)
      deallocate(receivers)
    end do
    ! A fixed-grid sampler can pack owner edges directly. Check once per
    ! plan, so all ranks make the same decision about the first collective.
    sync_edges_before_exchange=.true.
    if(refine_max_level==1) then
      local_needs_sync=.false.
      associate(request=>runtime_ct(1)%edges%edge%samples%request)
        do n=1,size(request,2)
          idx=request(2:4,n)+[ixMlo1,ixMlo2,ixMlo3]
          if(.not.yy_ct_edge_is_owner(request(1,n),idx,request(5,n))) then
            local_needs_sync=.true.
            exit
          end if
        end do
      end associate
      call MPI_ALLREDUCE(local_needs_sync,sync_edges_before_exchange,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierrmpi)
    end if
    runtime_ct_epoch=particle_mesh_epoch
    }
  end subroutine

  subroutine yy_runtime_ct_edges(s,fluxes)
    type(state),intent(in) :: s(max_blocks)
    type(ct_stage_flux),intent(inout) :: fluxes(:)
    integer :: l
    call yy_runtime_ct_prepare(s)
    if(sync_edges_before_exchange) then
      do l=1,size(runtime_ct)
        if(runtime_ct(l)%patch%mesh_epoch<0) cycle
        call yy_sync_ct_edges(runtime_ct(l)%patch,fluxes)
      end do
    end if
    do l=1,size(runtime_ct)
      if(runtime_ct(l)%patch%mesh_epoch<0) cycle
      call yy_exchange_ct_edges(runtime_ct(l)%edges,fluxes,defer_write=.true.)
    end do
    ! A fine receiver may use a coarse donor at the same stage and vice versa.
    ! Commit no overset edge until every level has packed its native integrals.
    do l=1,size(runtime_ct)
      if(runtime_ct(l)%patch%mesh_epoch<0) cycle
      call yy_write_ct_edges(runtime_ct(l)%edges,fluxes)
      call yy_sync_ct_edges(runtime_ct(l)%patch,fluxes)
    end do
  end subroutine

  subroutine yy_runtime_ct_ghosts(s,sibling_faces_filled)
    type(state),intent(inout) :: s(max_blocks)
    logical,optional,intent(in) :: sibling_faces_filled
    integer :: l
    logical :: sync_interior
    sync_interior=.true.
    ! After fixed-grid sibling exchange, interpolation reads active donor faces;
    ! angular ghost closure does not depend on interior ghost copies.
    if(present(sibling_faces_filled)) sync_interior=refine_max_level>1.or..not.sibling_faces_filled
    if(.not.associated(phys_ct_sync_center)) call mpistop('CT variable sync procedure is missing')
    call yy_runtime_ct_prepare(s)
    do l=1,size(runtime_ct)
      if(runtime_ct(l)%patch%mesh_epoch<0) cycle
      if(sync_interior) call yy_sync_ct_faces(runtime_ct(l)%patch,s,.false.)
      call yy_fill_ct_ghosts(runtime_ct(l)%ghosts,s,.true.,level=l)
      call yy_sync_ct_faces(runtime_ct(l)%patch,s,.true.)
    end do
    call yy_ct_radial_faces(s)
    call yy_sync_ct_state(s,phys_ct_sync_center)
  end subroutine

  ! Supply radial cell donors before angular primitive interpolation. Final
  ! radial B/energy is recomputed only after all face and angular data are ready.
  subroutine yy_ct_radial_cells(s,time,qdt)
    use mod_usr_methods, only: usr_special_bc,usr_prepare_boundary
    double precision,intent(in) :: time,qdt
    type(state),target,intent(inout) :: s(max_blocks)
    {^IFTHREED
    integer :: n,g,side,anchor,lo,hi,i,ixO^L
    if(special_inner_ct.and.associated(usr_prepare_boundary)) call usr_prepare_boundary(time,qdt)
    do n=1,igridstail
      g=igrids(n)
      do side=1,2
        if(side==1.and.node(pig1_,g)/=1) cycle
        if(side==2.and.node(pig1_,g)/=ng1(node(plevel_,g))) cycle
        anchor=ixMlo1;lo=ixGlo1;hi=ixMlo1-1
        if(side==2) then
          anchor=ixMhi1;lo=ixMhi1+1;hi=ixGhi1
        end if
        do i=lo,hi
          s(g)%w(i,ixGlo2:ixGhi2,ixGlo3:ixGhi3,:)=s(g)%w(anchor,ixGlo2:ixGhi2,ixGlo3:ixGhi3,:)
        end do
        if(special_inner_ct.and.side==1) then
          block=>s(g)
          ^D&dxlevel(^D)=rnode(rpdx^D_,g);
          ixO^L=ixG^LL;ixOmax1=hi
          call usr_special_bc(qdt,time,ixG^LL,ixO^L,1,s(g)%w,s(g)%x)
        end if
      end do
    end do
    }
  end subroutine

  ! Invalidate a previous completion record BEFORE overwriting either snapshot.
  ! Called collectively on the world communicator, outside block/OpenMP loops.
  subroutine yy_checkpoint_begin(filename)
    character(len=*),intent(in) :: filename
    integer :: u,ios
    logical :: exists
    if(.not.yinyang_enabled) return
    if(mype==0) then
      inquire(file=trim(filename)//'.yy.complete',exist=exists)
      if(exists) then
        open(newunit=u,file=trim(filename)//'.yy.complete',status='old',iostat=ios)
        if(ios/=0) call mpistop('cannot invalidate yinyang checkpoint completion')
        close(u,status='delete',iostat=ios)
        if(ios/=0) call mpistop('cannot delete yinyang checkpoint completion')
      end if
    end if
    call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)
  end subroutine

  ! Version 2 completion records are published only after both data files and
  ! sidecars are closed. Size checks detect truncation, not same-size corruption.
  subroutine yy_checkpoint(filename,writing)
    character(len=*),intent(in) :: filename
    logical,intent(in) :: writing
    integer :: u,ios,header(8),expected(8),imin,imax,version,saved_ranks,ranks_min,ranks_max
    integer :: stored_version,version_min,version_max
    logical :: allow_repartition
    integer(kind=8) :: bytes,saved_bytes
    double precision :: stamp,tmin,tmax
    if(.not.yinyang_enabled) return
    {^IFTHREED
    expected=[2,yy_patch,yy_world_size,domain_nx1,domain_nx2,domain_nx3, &
      merge(1,0,stagger_grid),merge(1,0,B0field)]
    call MPI_ALLREDUCE(it,imin,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierrmpi)
    call MPI_ALLREDUCE(it,imax,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
    call MPI_ALLREDUCE(global_time,tmin,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierrmpi)
    call MPI_ALLREDUCE(global_time,tmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
    if(imin/=imax.or.tmin/=tmax) call mpistop('yinyang checkpoint patch times disagree')
    saved_ranks=yy_world_size
    stored_version=2
    allow_repartition=yy_ct_experimental.and.refine_max_level>1
    if(mype==0) then
      inquire(file=trim(filename),size=bytes,iostat=ios)
      if(ios/=0.or.bytes<=0) call mpistop('invalid yinyang checkpoint data size')
      if(writing) then
        open(newunit=u,file=trim(filename)//'.yy',status='replace',iostat=ios)
        if(ios/=0) call mpistop('cannot write yinyang checkpoint metadata')
        write(u,*,iostat=ios) expected
        if(ios/=0) call mpistop('cannot write yinyang checkpoint header')
        write(u,*,iostat=ios) it,global_time,bytes
        if(ios/=0) call mpistop('cannot write yinyang checkpoint time and size')
      else
        open(newunit=u,file=trim(filename)//'.yy',status='old',iostat=ios)
        if(ios/=0) call mpistop('missing yinyang checkpoint metadata')
        read(u,*,iostat=ios) version
        if(ios/=0) call mpistop('invalid yinyang checkpoint metadata')
        rewind(u)
        stored_version=version
        if(version==1) then
          read(u,*,iostat=ios) header(1:6)
          if(ios/=0) call mpistop('invalid legacy yinyang checkpoint metadata')
          if(any(header(2:6)/=expected(2:6))) call mpistop('legacy yinyang checkpoint mesh or MPI mismatch')
          read(u,*,iostat=ios) imin,stamp
          if(ios/=0) call mpistop('invalid legacy yinyang checkpoint time')
          if(yy_world_rank==0) write(*,'(a)') &
            'Yin-Yang: legacy v1 checkpoint has no completion or file-size guarantee'
        else if(version==2) then
          read(u,*,iostat=ios) header
          if(ios/=0) call mpistop('invalid yinyang checkpoint header')
          if(allow_repartition) then
            saved_ranks=header(3)
            if(saved_ranks<2.or.mod(saved_ranks,2)/=0) call mpistop('invalid yinyang checkpoint MPI size')
            expected(3)=saved_ranks
          end if
          if(any(header/=expected)) call mpistop('yinyang checkpoint mesh, patch, MPI or field mode mismatch')
          read(u,*,iostat=ios) imin,stamp,saved_bytes
          if(ios/=0) call mpistop('invalid yinyang checkpoint time and size')
          if(saved_bytes/=bytes) call mpistop('yinyang checkpoint data size mismatch')
        else
          call mpistop('unsupported yinyang checkpoint version')
        end if
        if(imin/=it.or.stamp/=global_time) call mpistop('yinyang checkpoint timestamp mismatch')
      end if
      close(u,iostat=ios)
      if(ios/=0) call mpistop('cannot close yinyang checkpoint metadata')
    end if
    if(.not.writing.and.allow_repartition) then
      call MPI_BCAST(stored_version,1,MPI_INTEGER,0,icomm,ierrmpi)
      call MPI_BCAST(saved_ranks,1,MPI_INTEGER,0,icomm,ierrmpi)
      call MPI_ALLREDUCE(stored_version,version_min,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierrmpi)
      call MPI_ALLREDUCE(stored_version,version_max,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
      call MPI_ALLREDUCE(saved_ranks,ranks_min,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierrmpi)
      call MPI_ALLREDUCE(saved_ranks,ranks_max,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
      if(ranks_min/=ranks_max.or.version_min/=version_max) &
        call mpistop('paired yinyang checkpoint MPI size or versions disagree')
    end if
    call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)
    if(mype==0) then
      if(writing) then
        open(newunit=u,file=trim(filename)//'.yy.complete',status='replace',iostat=ios)
        if(ios/=0) call mpistop('cannot publish yinyang checkpoint completion')
        write(u,*,iostat=ios) it,global_time,bytes
        if(ios/=0) call mpistop('cannot write yinyang checkpoint completion')
        close(u,iostat=ios)
        if(ios/=0) call mpistop('cannot close yinyang checkpoint completion')
      else if(version==2) then
        open(newunit=u,file=trim(filename)//'.yy.complete',status='old',iostat=ios)
        if(ios/=0) call mpistop('missing yinyang checkpoint completion')
        read(u,*,iostat=ios) imin,stamp,saved_bytes
        if(ios/=0) call mpistop('invalid yinyang checkpoint completion')
        close(u)
        if(imin/=it.or.stamp/=global_time.or.saved_bytes/=bytes) &
          call mpistop('stale yinyang checkpoint completion')
      end if
    end if
    call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)
    }
  end subroutine

  subroutine yy_report()
    double precision :: local(4),maximum(4),average(4)
    integer(kind=8) :: samples(2),totals(2),calls_min,calls_max
    if(.not.yinyang_enabled) return
    local=[build_seconds,pack_seconds,wait_seconds,unpack_seconds]
    call MPI_REDUCE(local,maximum,4,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierrmpi)
    call MPI_REDUCE(local,average,4,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierrmpi)
    samples=[int(raw_nr,8),int(nr,8)]
    call MPI_REDUCE(samples,totals,2,MPI_INTEGER8,MPI_SUM,0,MPI_COMM_WORLD,ierrmpi)
    call MPI_REDUCE(exchange_calls,calls_min,1,MPI_INTEGER8,MPI_MIN,0,MPI_COMM_WORLD,ierrmpi)
    call MPI_REDUCE(exchange_calls,calls_max,1,MPI_INTEGER8,MPI_MAX,0,MPI_COMM_WORLD,ierrmpi)
    if(yy_world_rank==0) then
      if(refine_max_level>1) write(*,'(a)') 'Yin-Yang level-one exchange statistics'
      average=average/yy_world_size
      write(*,'(a,4es14.6)') 'Yin-Yang max-rank build / pack / exchange / unpack seconds: ',maximum
      write(*,'(a,4es14.6)') 'Yin-Yang mean-rank build / pack / exchange / unpack seconds: ',average
      write(*,'(a,2i14)') 'Yin-Yang donor samples raw / unique per exchange: ',totals
      write(*,'(a,2i14)') 'Yin-Yang payload bytes raw / unique per exchange: ',totals*int(nw,8)*8_8
      write(*,'(a,2i14)') 'Yin-Yang exchange calls min / max: ',calls_min,calls_max
    end if
    if(refine_max_level>1) call yy_glm_report()
  end subroutine

  ! Deduplicate requests for each sending peer. Only schedule construction uses
  ! hashing; runtime slot lookups and their four-term summation order are unchanged.
  subroutine yy_unique_requests(outreq)
    integer,allocatable,intent(inout) :: outreq(:,:)
    integer,allocatable :: compact(:,:),remap(:),table(:),peer(:),unique_counts(:)
    integer :: capacity,p,j,q,h,found,nunique,t
    integer(kind=8) :: key
    raw_nr=nr
    if(nr>(huge(nr)-1)/2) call mpistop('yinyang request table exceeds integer capacity')
    capacity=max(1,2*nr+1)
    allocate(compact(4,nr),remap(nr),table(capacity),peer(nr),unique_counts(yy_world_size))
    table=0;unique_counts=0;nunique=0
    do p=1,yy_world_size
      do j=sd(p)+1,sd(p)+sc(p)
        key=int(p,8)
        do q=1,4
          key=modulo(key*1000003_8+int(outreq(q,j),8),int(capacity,8))
        end do
        h=int(key)+1
        do
          found=table(h)
          if(found==0) then
            nunique=nunique+1
            compact(:,nunique)=outreq(:,j);peer(nunique)=p
            table(h)=nunique;remap(j)=nunique
            unique_counts(p)=unique_counts(p)+1
            exit
          end if
          if(peer(found)==p) then
            if(all(compact(:,found)==outreq(:,j))) then
              remap(j)=found
              exit
            end if
          end if
          h=mod(h,capacity)+1
        end do
      end do
    end do
    do t=1,nt
      do q=1,4
        slot(q,t)=remap(slot(q,t))
      end do
    end do
    sc=unique_counts;sd(1)=0
    do p=2,yy_world_size
      sd(p)=sd(p-1)+sc(p-1)
    end do
    nr=nunique
    deallocate(outreq)
    allocate(outreq(4,nr))
    outreq=compact(:,1:nr)
  end subroutine

  pure subroutine yy_position(r,theta,phi,patch,v)
    double precision,intent(in) :: r,theta,phi
    integer,intent(in) :: patch
    double precision,intent(out) :: v(3)
    double precision :: b(3,3)
    call yy_basis(theta,phi,patch,b)
    v=r*b(:,1)
  end subroutine

  subroutine yy_validate()
    use mod_usr_methods, only: usr_special_bc,usr_set_electric_field
    integer :: iw
    if(yy_ct_face_moments.and..not.yy_ct_experimental)call mpistop('face moments require experimental YY CT')
    if(yy_ct_experimental.and..not.yinyang_enabled) call mpistop('YY CT requires yinyang coordinates')
    if(.not.yinyang_enabled) return
    {^IFTHREED
    if(yy_ct_experimental) then
      if(.not.stagger_grid.or.physics_type/='mhd'.or..not.ct_stage_split.or. &
        .not.associated(phys_ct_sync_center).or..not.associated(phys_get_ct_edges).or. &
        .not.associated(phys_apply_ct_edges)) call mpistop('YY CT needs staged classical UCT-HLL with cell-center synchronization')
      if(dimsplit.or.t_stepper/=threestep.or.t_integrator/=ssprk3.or.internalboundary) &
        call mpistop('YY CT needs unsplit SSPRK3; no internal boundaries')
      if(any(typeboundary(:,2)/=bc_cont)) call mpistop('YY CT requires a continuous outer boundary')
      special_inner_ct=any(typeboundary(:,1)==bc_special)
      do iw=1,nw
        if(special_inner_ct.and..not.any(iw_mag==iw)) then
          if(typeboundary(iw,1)/=bc_special) call mpistop('YY CT inner fluid variables must all be special')
        else
          if(typeboundary(iw,1)/=bc_cont) call mpistop('YY CT requires continuous inner magnetic boundaries')
        end if
      end do
      if(special_inner_ct.and.(.not.associated(usr_special_bc).or..not.associated(usr_set_electric_field))) &
        call mpistop('YY CT special inner boundary requires fluid and electric field callbacks')
    else if(stagger_grid) then
      call mpistop('Yin-Yang CT requires explicit yy_ct_experimental=T')
    end if
    if(lb_automatic.or.local_timestep) &
      call mpistop('Yin-Yang requires global time steps and no cost-weighted balancing')
    if(refine_max_level>1) then
      if(yy_ct_experimental.and.any(int([domain_nx1,domain_nx2,domain_nx3],8) &
        *2_8**(refine_max_level-1)>huge(1))) call mpistop('Yin-Yang CT AMR cell indices exceed capacity')
      if(physics_type/='mhd'.or.internalboundary) &
        call mpistop('Yin-Yang AMR requires MHD without internal boundaries')
      if(any(int([domain_nx1/block_nx1,domain_nx2/block_nx2,domain_nx3/block_nx3],8) &
        *2_8**(refine_max_level-1)>huge(1))) call mpistop('Yin-Yang AMR native block indices exceed capacity')
      if(2*nghostcells>min(block_nx1,block_nx2,block_nx3)) &
        call mpistop('Yin-Yang AMR requires block widths of at least twice the ghost width')
      if(level_io>0.or.level_io_min/=1.or.level_io_max/=nlevelshi) &
        call mpistop('Yin-Yang AMR does not support forced output refinement levels')
    end if
    if(B0field.or.use_particles.or.use_multigrid) &
      call mpistop('yinyang v1: B0, particles and multigrid are unsupported')
    if(any(periodB).or.any(aperiodB)) call mpistop('yinyang has no periodic boundary')
    if(xprobmin1<=0.d0.or.xprobmin2<=0.d0.or.xprobmax2>=dpi) &
      call mpistop('yinyang requires a positive-radius shell excluding polar axes')
    if(stretched_dim(2).or.stretched_dim(3)) call mpistop('yinyang angular mesh must be uniform')
    if(physics_type/='mhd'.and.physics_type/='hd'.and.physics_type/='rho') &
      call mpistop('yinyang supports rho, hd and ideal GLM mhd only')
    if(nw/=nwgc.or.iwstart/=1.or.nwaux>0.or.nwextra>0) &
      call mpistop('yinyang requires complete primary-state exchange')
    if(yy_world_rank==0) write(*,'(a)') 'Yin-Yang v1: interpolation coupling (not flux conservative)'
    }
  end subroutine

  ! Called by both patches after the ordinary AMR mesh update. Layout epochs
  ! are patch-local and may diverge even though the next exchange is collective.
  subroutine yy_invalidate_after_regrid()
    integer :: changed,any_changed
    changed=merge(1,0,epoch/=particle_mesh_epoch)
    call MPI_ALLREDUCE(changed,any_changed,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
    if(any_changed/=0) then
      epoch=-1_8
      runtime_ct_epoch=-1_8
      {^IFTHREED ct_directory_epoch=-1_8 }
      glm_mode_known=.false.
      call yy_glm_invalidate()
    end if
  end subroutine

  ! All ranks use the level-one schedule while both patches are fully coarse.
  subroutine yy_refresh_glm_mode(mixed)
    logical,intent(out) :: mixed
    integer :: highest
    mixed=.false.
    if(.not.yinyang_enabled.or.refine_max_level==1) return
    if(.not.glm_mode_known) then
      call MPI_ALLREDUCE(levmax,highest,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
      glm_mode_mixed=highest>1
      glm_mode_known=.true.
      ! The mesh epoch is consumed by whichever schedule is active.
      if(glm_mode_mixed) epoch=particle_mesh_epoch
    end if
    mixed=glm_mode_mixed
  end subroutine

  ! Intersect this level-one block with the angular domain.
  subroutine yy_angular_interior(g,lo,hi)
    integer,intent(in) :: g
    integer,intent(out) :: lo(2),hi(2)
    {^IFTHREED
    integer :: offset(2)
    offset=([node(pig2_,g),node(pig3_,g)]-1)*[block_nx2,block_nx3]
    lo=[ixMlo2,ixMlo3]-offset
    hi=lo+[domain_nx2,domain_nx3]-1
    }
  end subroutine

  subroutine yy_build(s)
    type(state),intent(in) :: s(max_blocks)
    {^IFTHREED
    integer :: nblocks,nlocal,n,g,a,b,j,k,i,q,t,p,off,dj,dk,bi,bj,bk,ci,cj,ck,pass,gi,lo(2),hi(2)
    integer,allocatable :: counts(:),displs(:),local(:,:),allblocks(:,:),owners(:,:,:,:),ids(:,:,:,:)
    integer,allocatable :: outreq(:,:),cursor(:),sc4(:),rc4(:),sd4(:),rd4(:)
    double precision :: v(3),theta,phi,aa,bb,wa,wb,dteta,dphi,mat(3,3)
    if(allocated(rec)) then
      deallocate(rec,slot,request,sc,rc,sd,rd,handles,statuses,weights,rbasis,sbasis,sendbuf,recvbuf,scratch,shead,snext,rhead,rnext)
    end if
    allocate(counts(yy_world_size),displs(yy_world_size))
    nlocal=igridstail
    if(any(node(plevel_,igrids(1:nlocal))/=1)) call mpistop('Yin-Yang level-one exchange needs a fully coarse mesh')
    allocate(local(5,nlocal))
    do n=1,nlocal
      g=igrids(n)
      local(:,n)=[node(pig1_,g),node(pig2_,g),node(pig3_,g),g,yy_world_rank]
    end do
    call MPI_ALLGATHER(5*nlocal,1,MPI_INTEGER,counts,1,MPI_INTEGER,MPI_COMM_WORLD,ierrmpi)
    displs(1)=0
    do p=2,yy_world_size
      displs(p)=displs(p-1)+counts(p-1)
    end do
    nblocks=sum(counts)/5
    allocate(allblocks(5,nblocks))
    call MPI_ALLGATHERV(local,5*nlocal,MPI_INTEGER,allblocks,counts,displs,MPI_INTEGER,MPI_COMM_WORLD,ierrmpi)
    allocate(owners(ng1(1),ng2(1),ng3(1),0:1),ids(ng1(1),ng2(1),ng3(1),0:1))
    owners=-1; ids=-1
    do n=1,nblocks
      p=allblocks(5,n)/(yy_world_size/2)
      bi=allblocks(1,n);bj=allblocks(2,n);bk=allblocks(3,n)
      owners(bi,bj,bk,p)=allblocks(5,n)
      ids(bi,bj,bk,p)=allblocks(4,n)
    end do
    if(any(owners<0)) call mpistop('yinyang missing level-one block')
    allocate(sc(yy_world_size),rc(yy_world_size),sd(yy_world_size),rd(yy_world_size),cursor(yy_world_size))
    sc=0; nt=0
    dteta=(xprobmax2-xprobmin2)/domain_nx2
    dphi=(xprobmax3-xprobmin3)/domain_nx3
    do pass=1,2
      t=0
      do n=1,nlocal
        g=local(4,n)
        call yy_angular_interior(g,lo,hi)
        if(all(lo<=[ixGlo2,ixGlo3]).and.all(hi>=[ixGhi2,ixGhi3])) cycle
        do k=ixGlo3,ixGhi3
          do j=ixGlo2,ixGhi2
            if(j>=lo(1).and.j<=hi(1).and.k>=lo(2).and.k<=hi(2)) cycle
            theta=s(g)%x(ixMlo1,j,k,2); phi=s(g)%x(ixMlo1,j,k,3)
            v=[-sin(theta)*cos(phi),cos(theta),sin(theta)*sin(phi)]
            theta=atan2(sqrt(v(1)**2+v(2)**2),v(3));phi=atan2(v(2),v(1))
            aa=(theta-xprobmin2)/dteta-0.5d0;bb=(phi-xprobmin3)/dphi-0.5d0
            dj=floor(aa);dk=floor(bb);wa=aa-dj;wb=bb-dk
            if(dj<0.or.dj+1>=domain_nx2.or.dk<0.or.dk+1>=domain_nx3) &
              call mpistop('yinyang ghost lacks interior donors: increase angular overlap')
            do i=ixGlo1,ixGhi1
              gi=(node(pig1_,g)-1)*block_nx1+i-ixMlo1
              bi=max(1,min(ng1(1),gi/block_nx1+1));ci=gi-(bi-1)*block_nx1+ixMlo1
              t=t+1
              if(pass==2) then
                rec(:,t)=[g,i,j,k]
                call yy_basis(s(g)%x(i,j,k,2),s(g)%x(i,j,k,3),yy_patch,mat)
                rbasis(:,:,t)=transpose(mat)
              end if
              q=0
              do b=0,1
                do a=0,1
                  q=q+1
                  bj=(dj+a)/block_nx2+1;bk=(dk+b)/block_nx3+1
                  cj=mod(dj+a,block_nx2)+ixMlo2;ck=mod(dk+b,block_nx3)+ixMlo3
                  p=owners(bi,bj,bk,1-yy_patch)+1
                  if(pass==1) then
                    sc(p)=sc(p)+1
                  else
                    off=sd(p)+cursor(p)+1;cursor(p)=cursor(p)+1
                    slot(q,t)=off
                    weights(q,t)=merge(1.d0-wa,wa,a==0)*merge(1.d0-wb,wb,b==0)
                    outreq(:,off)=[ids(bi,bj,bk,1-yy_patch),ci,cj,ck]
                  end if
                end do
              end do
            end do
          end do
        end do
      end do
      if(pass==1) then
        nt=t;nr=sum(sc)
        allocate(rec(4,nt),slot(4,nt),weights(4,nt),rbasis(3,3,nt),outreq(4,nr))
        sd(1)=0
        do p=2,yy_world_size
          sd(p)=sd(p-1)+sc(p-1)
        end do
        cursor=0
      end if
    end do
    call yy_unique_requests(outreq)
    call MPI_ALLTOALL(sc,1,MPI_INTEGER,rc,1,MPI_INTEGER,MPI_COMM_WORLD,ierrmpi)
    rd(1)=0
    do p=2,yy_world_size
      rd(p)=rd(p-1)+rc(p-1)
    end do
    ns=sum(rc)
    allocate(request(4,ns),sbasis(3,3,ns),sendbuf(nw,ns),recvbuf(nw,nr))
    allocate(sc4(yy_world_size),rc4(yy_world_size),sd4(yy_world_size),rd4(yy_world_size))
    sc4=4*sc;rc4=4*rc;sd4=4*sd;rd4=4*rd
    call MPI_ALLTOALLV(outreq,sc4,sd4,MPI_INTEGER,request,rc4,rd4,MPI_INTEGER,MPI_COMM_WORLD,ierrmpi)
    do t=1,ns
      g=request(1,t);i=request(2,t);j=request(3,t);k=request(4,t)
      call yy_basis(s(g)%x(i,j,k,2),s(g)%x(i,j,k,3),yy_patch,sbasis(:,:,t))
    end do
    allocate(scratch(ixG^T,nw),handles(2*yy_world_size),statuses(MPI_STATUS_SIZE,2*yy_world_size))
    allocate(shead(max_blocks),rhead(max_blocks),snext(ns),rnext(nt))
    shead=0;rhead=0
    do t=ns,1,-1
      g=request(1,t);snext(t)=shead(g);shead(g)=t
    end do
    do t=nt,1,-1
      g=rec(1,t);rnext(t)=rhead(g);rhead(g)=t
    end do
    if(yy_world_rank==0.and.(refine_max_level==1.or.exchange_calls==0_8)) &
      write(*,'(a,i0,a,i0)') 'Yin-Yang schedule: rank 0 receivers=',nt,', samples=',ns
    epoch=particle_mesh_epoch
    }
  end subroutine

  subroutine yy_exchange(s,nwstart,nwbc)
    type(state),target,intent(inout) :: s(max_blocks)
    integer,intent(in) :: nwstart,nwbc
    {^IFTHREED
    integer :: n,g,t,i,j,k,q,iv,iw,p,nh,side,axis,ixO^L,lo(3),hi(3),angular_lo(2),angular_hi(2)
    double precision :: v(3),val(nw),t0
    if(nwstart/=1.or.nwbc/=nw) call mpistop('yinyang partial-state exchange unsupported')
    if(epoch==-1_8.or.(refine_max_level==1.and.epoch/=particle_mesh_epoch)) then
      t0=MPI_WTIME()
      call yy_build(s)
      build_seconds=build_seconds+MPI_WTIME()-t0
    end if
    exchange_calls=exchange_calls+1_8
    t0=MPI_WTIME()
    do n=1,igridstail
      g=igrids(n)
      if(shead(g)==0) cycle
      block=>s(g)
      scratch=s(g)%w
      ixO^L=ixG^LL;
      call yy_angular_interior(g,angular_lo,angular_hi)
      ixOmin2=max(ixGlo2,angular_lo(1));ixOmax2=min(ixGhi2,angular_hi(1))
      ixOmin3=max(ixGlo3,angular_lo(2));ixOmax3=min(ixGhi3,angular_hi(2))
      call phys_to_primitive(ixG^LL,ixO^L,scratch,s(g)%x)
      t=shead(g)
      do while(t/=0)
        i=request(2,t);j=request(3,t);k=request(4,t)
        sendbuf(:,t)=scratch(i,j,k,:)
        do iv=1,nvector
          iw=iw_vector(iv);v=scratch(i,j,k,iw+1:iw+3)
          sendbuf(iw+1:iw+3,t)=matmul(sbasis(:,:,t),v)
        end do
        t=snext(t)
      end do
    end do
    pack_seconds=pack_seconds+MPI_WTIME()-t0
    nh=0;t0=MPI_WTIME()
    do p=1,yy_world_size
      if(sc(p)>0) then
        nh=nh+1
        call MPI_IRECV(recvbuf(1,sd(p)+1),nw*sc(p),MPI_DOUBLE_PRECISION,p-1,27131,MPI_COMM_WORLD,handles(nh),ierrmpi)
      end if
      if(rc(p)>0) then
        nh=nh+1
        call MPI_ISEND(sendbuf(1,rd(p)+1),nw*rc(p),MPI_DOUBLE_PRECISION,p-1,27131,MPI_COMM_WORLD,handles(nh),ierrmpi)
      end if
    end do
    if(nh>0) call MPI_WAITALL(nh,handles,statuses,ierrmpi)
    wait_seconds=wait_seconds+MPI_WTIME()-t0
    t0=MPI_WTIME()
    do n=1,igridstail
      g=igrids(n)
      if(rhead(g)==0) cycle
      block=>s(g)
      t=rhead(g)
      do while(t/=0)
        i=rec(2,t);j=rec(3,t);k=rec(4,t)
        val=0.d0
        do q=1,4
          val=val+weights(q,t)*recvbuf(:,slot(q,t))
        end do
        do iv=1,nvector
          iw=iw_vector(iv);v=val(iw+1:iw+3)
          val(iw+1:iw+3)=matmul(rbasis(:,:,t),v)
        end do
        scratch(i,j,k,:)=val
        t=rnext(t)
      end do
      ! Convert rectangular ghost slabs in bulk; every read cell was filled above.
      call yy_angular_interior(g,angular_lo,angular_hi)
      do axis=2,3
        do side=1,2
          lo=[ixGlo1,ixGlo2,ixGlo3];hi=[ixGhi1,ixGhi2,ixGhi3]
          if(side==1) then
            hi(axis)=angular_lo(axis-1)-1
          else
            lo(axis)=angular_hi(axis-1)+1
          end if
          if(axis==3) then
            lo(2)=max(ixGlo2,angular_lo(1))
            hi(2)=min(ixGhi2,angular_hi(1))
          end if
          if(any(lo>hi)) cycle
          ixOmin1=lo(1);ixOmin2=lo(2);ixOmin3=lo(3)
          ixOmax1=hi(1);ixOmax2=hi(2);ixOmax3=hi(3)
          call phys_to_conserved(ixG^LL,ixO^L,scratch,s(g)%x)
          s(g)%w(ixO^S,:)=scratch(ixO^S,:)
        end do
      end do
    end do
    unpack_seconds=unpack_seconds+MPI_WTIME()-t0
    }
  end subroutine
end module mod_yinyang
