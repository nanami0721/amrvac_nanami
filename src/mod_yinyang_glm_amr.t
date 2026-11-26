!> Sparse overset geometry and exchanges for adaptive cell-centered GLM.
module mod_yinyang_glm_amr
  use, intrinsic :: iso_fortran_env, only: int64
  use mod_global_parameters
  use mod_comm_lib, only: mpistop
  use mod_forest, only: level_head, tree_node, igrid_to_node, refine, coarsen
  use mod_physics, only: phys_to_primitive, phys_to_conserved
  implicit none
  private
  public :: yy_glm_balance, yy_glm_exchange, yy_glm_basis, yy_glm_invalidate, yy_glm_report
  {^IFTHREED
  public :: leaf_directory, index_directory, covering_leaf, domain_cells, radial_edge
  public :: build_directory, find_leaf, key_level, key_index, key_rank, key_grid, check_level_pair
  integer, parameter :: ndata=7, key_level=1, key_index=2, key_rank=5, key_grid=6, key_action=7
  integer, parameter :: unchanged=0, split_leaf=1, merge_siblings=-1, deferred_leaf=2
  integer, parameter :: change_coarse=1, defer_split=2
  integer, parameter :: primitive_sample=0, conserved_sample=-1
  type leaf_directory
    ! Rows are level, native i/j/k block indices, world rank, grid ID and edit action.
    integer, allocatable :: leaf(:,:), bucket(:)
    integer :: deepest
  end type
  abstract interface
    subroutine check_seam_levels(dir,want)
      import leaf_directory
      type(leaf_directory),intent(in) :: dir
      integer,intent(inout) :: want(:)
    end subroutine
  end interface

  integer, parameter :: query_size=8
  type mixed_plan
    integer,allocatable :: receivers(:,:), ends(:,:), slot(:)
    integer(int64),allocatable :: request(:,:)
    integer,allocatable :: send_count(:), recv_count(:), send_off(:), recv_off(:)
    integer,allocatable :: first(:), next(:), handles(:), statuses(:,:), source_slot(:), stencil_start(:,:)
    double precision,allocatable :: source_basis(:,:,:,:,:), target_basis(:,:,:), displacement(:,:), span_inv(:,:,:)
    double precision,allocatable :: restriction_center(:,:)
    double precision,allocatable :: weights(:,:), volume(:), send(:,:), recv(:,:)
  end type
  type(mixed_plan),save :: interior_plan,corner_plan
  double precision,allocatable,save :: primitive(:,:,:,:)
  integer(int64),save :: calls(2)=0_int64, samples(2)=0_int64, sites(4)=0_int64
  double precision,save :: seconds(4)=0.d0

  type cell_stencil
    integer(int64) :: site(3,4)
    double precision :: weight(4)
  end type
  }

contains

  pure subroutine yy_glm_basis(theta,phi,patch,b)
    double precision,intent(in) :: theta,phi
    integer,intent(in) :: patch
    double precision,intent(out) :: b(3,3)
    double precision :: tmp(3)
    b(:,1)=[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]
    b(:,2)=[cos(theta)*cos(phi),cos(theta)*sin(phi),-sin(theta)]
    b(:,3)=[-sin(phi),cos(phi),0.d0]
    if(patch==1) then
      b(1,:)=-b(1,:)
      tmp=b(2,:); b(2,:)=b(3,:); b(3,:)=tmp
    end if
  end subroutine

  subroutine yy_glm_balance(changed,deferred{^IFTHREED ,check_ct})
    logical,intent(out) :: changed
    logical,intent(inout) :: deferred(max_blocks,0:npe-1)
    {^IFTHREED
    procedure(check_seam_levels),optional :: check_ct
    }
    changed=.false.
    {^IFTHREED
    call balance_proposal(changed,deferred,check_ct)
    }
  end subroutine

  {^IFTHREED
  pure function domain_cells(level) result(nn)
    integer,intent(in) :: level
    integer(int64) :: nn(3)
    nn=int([domain_nx1,domain_nx2,domain_nx3],int64)*2_int64**(level-1)
  end function

  pure function block_cells() result(bn)
    integer :: bn(3)
    bn=[block_nx1,block_nx2,block_nx3]
  end function

  ! Directory keys use the native level and block tuple. Sparse storage grows
  ! with actual leaves, independent of the largest permitted level.
  integer function key_hash(level,idx,patch,capacity) result(h)
    integer,intent(in) :: level,idx(3),patch,capacity
    integer :: axis
    integer(kind=8) :: acc,limit
    limit=int(capacity,8)
    acc=int(2*level+patch,8)
    do axis=1,3
      acc=modulo(acc*1000003_8+int(idx(axis),8),limit)
    end do
    h=int(acc)+1
  end function

  subroutine build_directory(dir,deferred)
    type(leaf_directory),intent(out) :: dir
    logical,intent(in),optional :: deferred(max_blocks,0:npe-1)
    integer,allocatable :: local(:,:),counts(:),offsets(:)
    type(tree_node),pointer :: leaf,parent,first
    integer :: n,pass,level,l,nlocal,action,idx(3),rank
    do pass=1,2
      nlocal=0
      do level=levmin,levmax
        leaf=>level_head(level)%node
        do while(associated(leaf))
          if(leaf%ipe==mype) then
            l=level
            idx=[leaf%ig1,leaf%ig2,leaf%ig3]
            action=unchanged
            if(present(deferred)) then
              if(coarsen(leaf%igrid,mype)) then
                parent=>leaf%parent%node
                first=>parent%child(1,1,1)%node
                if(associated(first,leaf)) then
                  l=level-1
                  idx=[parent%ig1,parent%ig2,parent%ig3]
                  action=merge_siblings
                else
                  leaf=>leaf%next%node
                  cycle
                end if
              else if(deferred(leaf%igrid,mype)) then
                action=deferred_leaf
              else if(refine(leaf%igrid,mype)) then
                action=split_leaf
              end if
            end if
            if(action==split_leaf) then
              do n=0,7
                nlocal=nlocal+1
                if(pass==2) local(:,nlocal)=[level+1,2*idx-1+[iand(n,1),iand(ishft(n,-1),1),ishft(n,-2)], &
                  yy_world_rank,leaf%igrid,action]
              end do
            else
              nlocal=nlocal+1
              if(pass==2) local(:,nlocal)=[l,idx,yy_world_rank,leaf%igrid,action]
            end if
          end if
          leaf=>leaf%next%node
        end do
      end do
      if(pass==1) allocate(local(ndata,nlocal))
    end do

    allocate(counts(yy_world_size),offsets(yy_world_size))
    call MPI_ALLGATHER(ndata*nlocal,1,MPI_INTEGER,counts,1,MPI_INTEGER,MPI_COMM_WORLD,ierrmpi)
    offsets(1)=0
    do rank=2,yy_world_size
      offsets(rank)=offsets(rank-1)+counts(rank-1)
    end do
    n=sum(counts)/ndata
    allocate(dir%leaf(ndata,n))
    call MPI_ALLGATHERV(local,ndata*nlocal,MPI_INTEGER,dir%leaf,counts,offsets,MPI_INTEGER,MPI_COMM_WORLD,ierrmpi)
    call index_directory(dir)
  end subroutine

  subroutine index_directory(dir)
    type(leaf_directory),intent(inout) :: dir
    integer :: n,idx(3),h,patch
    n=size(dir%leaf,2)
    if(int(n,int64)*2>=huge(n)) call mpistop('Yin-Yang AMR leaf directory exceeds index capacity')
    allocate(dir%bucket(2*n+1))
    dir%bucket=0
    dir%deepest=maxval(dir%leaf(key_level,:))
    do n=1,size(dir%leaf,2)
      patch=dir%leaf(key_rank,n)/(yy_world_size/2)
      idx=dir%leaf(key_index:key_index+2,n)
      h=key_hash(dir%leaf(key_level,n),idx,patch,size(dir%bucket))
      do while(dir%bucket(h)/=0)
        h=mod(h,size(dir%bucket))+1
      end do
      dir%bucket(h)=n
    end do
  end subroutine

  integer function find_leaf(dir,level,idx,patch) result(found)
    type(leaf_directory),intent(in) :: dir
    integer,intent(in) :: level,idx(3),patch
    integer :: h,n
    h=key_hash(level,idx,patch,size(dir%bucket))
    found=0
    do
      n=dir%bucket(h)
      if(n==0) return
      if(patch==dir%leaf(key_rank,n)/(yy_world_size/2).and.dir%leaf(key_level,n)==level) then
        if(all(dir%leaf(key_index:key_index+2,n)==idx)) then
          found=n
          return
        end if
      end if
      h=mod(h,size(dir%bucket))+1
    end do
  end function

  ! The virtual cell lies in one native block even at level+1, because native
  ! block widths are even. For gaps larger than one, its low child is a witness
  ! of a too-fine donor and the proposed layout is rebuilt after promotion.
  integer function covering_leaf(dir,level,virtual,patch) result(found)
    type(leaf_directory),intent(in) :: dir
    integer,intent(in) :: level,patch
    integer(int64),intent(in) :: virtual(3)
    integer :: l,ig(3),bn(3)
    bn=block_cells()
    do l=min(level,dir%deepest),1,-1
      ig=int(virtual/(int(bn,int64)*2_int64**(level-l)))+1
      found=find_leaf(dir,l,ig,patch)
      if(found/=0) return
    end do
    do l=level+1,dir%deepest
      ig=int((2_int64**(l-level)*virtual)/bn)+1
      found=find_leaf(dir,l,ig,patch)
      if(found/=0) return
    end do
    call mpistop('Yin-Yang AMR virtual cell has no owning leaf')
  end function

  subroutine angular_stencil(level,global_j,global_k,stencil)
    integer,intent(in) :: level
    integer(int64),intent(in) :: global_j,global_k
    type(cell_stencil),intent(out) :: stencil
    integer :: a,b,q
    integer(int64) :: nn(3),j0,k0
    double precision :: theta,phi,t,p,ac,bc,dteta,dphi,v(3)
    nn=domain_cells(level)
    dteta=(xprobmax2-xprobmin2)/nn(2)
    dphi=(xprobmax3-xprobmin3)/nn(3)
    theta=xprobmin2+(global_j+0.5d0)*dteta
    phi=xprobmin3+(global_k+0.5d0)*dphi
    v=[-sin(theta)*cos(phi),cos(theta),sin(theta)*sin(phi)]
    t=(atan2(sqrt(v(1)**2+v(2)**2),v(3))-xprobmin2)/dteta-0.5d0
    p=(atan2(v(2),v(1))-xprobmin3)/dphi-0.5d0
    j0=floor(t,kind=int64); k0=floor(p,kind=int64)
    ac=t-j0; bc=p-k0
    if(j0<0.or.j0+1>=nn(2).or.k0<0.or.k0+1>=nn(3)) &
      call mpistop('Yin-Yang AMR ghost lacks angular interior donors: increase overlap')
    q=0
    do b=0,1
      do a=0,1
        q=q+1
        stencil%site(:,q)=[0_int64,j0+a,k0+b]
        stencil%weight(q)=merge(1.d0-ac,ac,a==0)*merge(1.d0-bc,bc,b==0)
      end do
    end do
  end subroutine

  subroutine balance_proposal(changed,deferred,check_ct)
    logical,intent(out) :: changed
    logical,intent(inout) :: deferred(max_blocks,0:npe-1)
    procedure(check_seam_levels),optional :: check_ct
    type(leaf_directory) :: dir
    type(cell_stencil) :: stencil
    integer,allocatable :: want(:),wanted(:)
    integer :: n,l,bn(3),j,k,q,other,local_rank,grid,ic,jc,kc
    integer(int64) :: nn(3),off(3),i
    type(tree_node),pointer :: first,child
    call build_directory(dir,deferred)
    changed=.false.
    if(yy_ct_experimental.and.dir%deepest-minval(dir%leaf(key_level,:))<=1) return
    allocate(want(size(dir%leaf,2)),wanted(size(dir%leaf,2)))
    want=0
    bn=block_cells()
    do n=1,size(dir%leaf,2)
      if(dir%leaf(key_rank,n)/=yy_world_rank) cycle
      l=dir%leaf(key_level,n)
      nn=domain_cells(l)
      off=(int(dir%leaf(key_index:key_index+2,n),int64)-1)*bn
      if(all(off(2:3)>=nghostcells).and.all(off(2:3)+bn(2:3)+nghostcells<=nn(2:3))) cycle
      do k=-nghostcells,bn(3)+nghostcells-1
        do j=-nghostcells,bn(2)+nghostcells-1
          if(off(2)+j>=0.and.off(2)+j<nn(2).and.off(3)+k>=0.and.off(3)+k<nn(3)) cycle
          call angular_stencil(l,off(2)+j,off(3)+k,stencil)
          do i=max(0_int64,off(1)-nghostcells),min(nn(1)-1,off(1)+bn(1)+nghostcells-1)
            do q=1,4
              stencil%site(1,q)=i
              other=covering_leaf(dir,l,stencil%site(:,q),1-yy_patch)
              call check_level_pair(dir,n,other,want)
            end do
          end do
        end do
      end do
    end do
    if(present(check_ct)) call check_ct(dir,want)
    call MPI_ALLREDUCE(want,wanted,size(want),MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierrmpi)
    changed=any(wanted/=0)
    do n=1,size(wanted)
      if(wanted(n)==0) cycle
      if(dir%leaf(key_rank,n)/(yy_world_size/2)/=yy_patch) cycle
      local_rank=mod(dir%leaf(key_rank,n),npe)
      grid=dir%leaf(key_grid,n)
      if(wanted(n)==defer_split) then
        refine(grid,local_rank)=.false.
        coarsen(grid,local_rank)=.false.
        deferred(grid,local_rank)=.true.
        cycle
      end if
      select case(dir%leaf(key_action,n))
      case(merge_siblings)
        first=>igrid_to_node(grid,local_rank)%node%parent%node
        do kc=1,2
          do jc=1,2
            do ic=1,2
              child=>first%child(ic,jc,kc)%node
              coarsen(child%igrid,child%ipe)=.false.
            end do
          end do
        end do
      case(unchanged)
        refine(grid,local_rank)=.true.
        coarsen(grid,local_rank)=.false.
      end select
    end do
  end subroutine
  subroutine check_level_pair(dir,n,other,want)
    type(leaf_directory),intent(in) :: dir
    integer,intent(in) :: n,other
    integer,intent(inout) :: want(:)
    integer :: l,l2,coarse,fine
    l=dir%leaf(key_level,n)
    l2=dir%leaf(key_level,other)
    if(abs(l-l2)<=1) return
    coarse=merge(n,other,l<l2)
    fine=merge(other,n,l<l2)
    select case(dir%leaf(key_action,coarse))
    case(unchanged,merge_siblings)
      want(coarse)=max(want(coarse),change_coarse)
    case(split_leaf,deferred_leaf)
      ! New geometry can expose a fine donor that the old stencil
      ! did not reach. Defer the split until its partner can catch up.
      select case(dir%leaf(key_action,fine))
      case(split_leaf)
        want(fine)=defer_split
      case(merge_siblings)
        want(fine)=max(want(fine),change_coarse)
      case default
        if(dir%leaf(key_action,coarse)/=split_leaf) &
          call mpistop('Yin-Yang AMR cannot balance existing seam leaves')
        want(coarse)=defer_split
      end select
    end select
  end subroutine

  }

  subroutine yy_glm_invalidate()
    {^IFTHREED
    type(mixed_plan) :: empty
    interior_plan=empty
    corner_plan=empty
    }
  end subroutine

  subroutine yy_glm_exchange(s,radial_only)
    type(state),target,intent(inout) :: s(max_blocks)
    logical,intent(in) :: radial_only
    {^IFTHREED
    type(leaf_directory) :: directory
    integer :: phase
    double precision :: before
    if(.not.allocated(interior_plan%receivers)) then
      before=MPI_WTIME()
      call build_directory(directory)
      call build_plan(interior_plan,directory,.false.)
      call cache_geometry(interior_plan,s)
      call build_plan(corner_plan,directory,.true.)
      call cache_geometry(corner_plan,s)
      if(.not.allocated(primitive)) allocate(primitive(ixG^T,nw))
      seconds(1)=seconds(1)+MPI_WTIME()-before
    end if
    phase=merge(2,1,radial_only)
    calls(phase)=calls(phase)+1_int64
    if(radial_only) then
      samples(phase)=samples(phase)+size(corner_plan%request,2,kind=int64)
      call exchange_plan(s,corner_plan,.true.)
    else
      samples(phase)=samples(phase)+size(interior_plan%request,2,kind=int64)
      call exchange_plan(s,interior_plan,.false.)
    end if
    }
  end subroutine

  subroutine yy_glm_report()
    {^IFTHREED
    double precision :: maximum(4)
    integer(int64) :: min_calls(2),max_calls(2),total_samples(2),total_sites(4)
    call MPI_REDUCE(calls,min_calls,2,MPI_INTEGER8,MPI_MIN,0,MPI_COMM_WORLD,ierrmpi)
    call MPI_REDUCE(calls,max_calls,2,MPI_INTEGER8,MPI_MAX,0,MPI_COMM_WORLD,ierrmpi)
    call MPI_REDUCE(samples,total_samples,2,MPI_INTEGER8,MPI_SUM,0,MPI_COMM_WORLD,ierrmpi)
    call MPI_REDUCE(sites,total_sites,4,MPI_INTEGER8,MPI_SUM,0,MPI_COMM_WORLD,ierrmpi)
    call MPI_REDUCE(seconds,maximum,4,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierrmpi)
    if(yy_world_rank==0) then
      write(*,'(a,4i14)') 'Yin-Yang mixed interior/corner exchange calls min then max: ',min_calls,max_calls
      write(*,'(a,2i14)') 'Yin-Yang mixed interior/corner samples sent: ',total_samples
      write(*,'(a,2i14)') 'Yin-Yang mixed interior/corner payload bytes: ',total_samples*nw*8_int64
      write(*,'(a,4i14)') 'Yin-Yang mixed layout sites coarse/same/fine and local corners: ',total_sites
      write(*,'(a,4es14.6)') 'Yin-Yang mixed max-rank build / pack / MPI / unpack seconds: ',maximum
    end if
    }
  end subroutine

  {^IFTHREED
  ! The radial faces match native spherical AMR, including both stretch modes.
  ! Integer edge positions let the planner identify descendants without
  ! creating a finest-level array or storing replicated block geometries.
  pure recursive function radial_edge(level,edge) result(r)
    integer,intent(in) :: level
    integer(int64),intent(in) :: edge
    integer(int64) :: left, right, nn(3)
    double precision :: r
    select case(stretch_type(1))
    case(stretch_uni)
      r=xprobmin1+dxfirst_1mq(level,1)*(1.d0-qstretch(level,1)**edge)
    case(stretch_symm)
      nn=domain_cells(level)
      left=int(nstretchedblocks(level,1),int64)*block_nx1/2
      right=nn(1)-left
      if(edge<left) then
        r=xprobmin1+xstretch1-dxfirst_1mq(level,1)*(1.d0-qstretch(level,1)**(left-edge))
      else if(edge>right) then
        r=xprobmin1+xstretch1+(right-left)*dxmid(level,1) &
          +dxfirst_1mq(level,1)*(1.d0-qstretch(level,1)**(edge-right))
      else
        r=xprobmin1+xstretch1+(edge-left)*dxmid(level,1)
      end if
    case default
      r=xprobmin1+edge*dx(1,level)
    end select
  end function

  pure function virtual_center(level,cell) result(center)
    integer,intent(in) :: level
    integer(int64),intent(in) :: cell(3)
    double precision :: center(3), left, right
    left=radial_edge(level,cell(1))
    right=radial_edge(level,cell(1)+1)
    center=[0.5d0*(left+right), &
      xprobmin2+(cell(2)+0.5d0)*dx(2,level),xprobmin3+(cell(3)+0.5d0)*dx(3,level)]
  end function

  pure function child_volume(level,cell) result(volume)
    integer,intent(in) :: level
    integer(int64),intent(in) :: cell(3)
    double precision :: left,right,dr,r,theta,volume
    left=radial_edge(level,cell(1))
    right=radial_edge(level,cell(1)+1)
    r=0.5d0*(left+right)
    dr=right-left
    theta=xprobmin2+(cell(2)+0.5d0)*dx(2,level)
    volume=dr*(r*r+dr*dr/12.d0)*sin(theta)
    ! The factors 2*sin(dtheta/2)*dphi cancel in the eight-cell mean.
  end function

  subroutine build_plan(plan,dir,radial_only)
    type(mixed_plan),intent(out) :: plan
    type(leaf_directory),intent(in) :: dir
    logical,intent(in) :: radial_only
    type(cell_stencil) :: stencil
    integer,allocatable :: cursor(:),sc(:),rc(:),sd(:),rd(:)
    integer(int64),allocatable :: outgoing(:,:)
    integer :: pass,n,g,l,ld,j,k,i,bn(3)
    integer(int64) :: nn(3),off(3),nnative(3),local(3)
    integer :: nt,nterm,t,q,leaves(4),peer,offset,ix_lo(3),nrestrict
    integer(int64) :: cell(3),child(3)
    integer :: u,v,w,mode,ns,nr,grid
    double precision :: volumes(2,2,2),total
    logical :: outside,corner_is_same
    bn=block_cells()
    ix_lo=[ixMlo1,ixMlo2,ixMlo3]
    allocate(sc(yy_world_size),rc(yy_world_size),sd(yy_world_size),rd(yy_world_size),cursor(yy_world_size))
    sc=0
    do pass=1,2
      nt=0
      nterm=0
      nrestrict=0
      do n=1,igridstail
        g=igrids(n)
        l=node(plevel_,g)
        nn=domain_cells(l)
        off=(int([node(pig1_,g),node(pig2_,g),node(pig3_,g)],int64)-1)*bn
        if(all(off(2:3)>=nghostcells).and.all(off(2:3)+bn(2:3)+nghostcells<=nn(2:3))) cycle
        if(radial_only.and.off(1)>=nghostcells.and.off(1)+bn(1)+nghostcells<=nn(1)) cycle
        do k=ixGlo3,ixGhi3
          do j=ixGlo2,ixGhi2
            local=off+[0,j-ixMlo2,k-ixMlo3]
            if(all(local(2:3)>=0).and.all(local(2:3)<nn(2:3))) cycle
            call angular_stencil(l,local(2),local(3),stencil)
            do i=ixGlo1,ixGhi1
              local(1)=off(1)+i-ixMlo1
              outside=local(1)<0.or.local(1)>=nn(1)
              if(outside.neqv.radial_only) cycle
              do q=1,4
                stencil%site(1,q)=max(0_int64,min(nn(1)-1,local(1)))
                leaves(q)=covering_leaf(dir,l,stencil%site(:,q),1-yy_patch)
              end do
              if(any(abs(dir%leaf(key_level,leaves)-l)>1)) &
                call mpistop('Yin-Yang AMR ghost donor exceeds one level')
              corner_is_same=all(dir%leaf(key_level,leaves)==l)
              if(outside.and..not.corner_is_same) then
                if(pass==2) sites(4)=sites(4)+1_int64
                cycle
              end if
              stencil%site(1,:)=local(1)
              nt=nt+1
              if(pass==2) then
                plan%receivers(:,nt)=[g,i,j,k]
                plan%weights(:,nt)=stencil%weight
              end if
              do q=1,4
                ld=dir%leaf(key_level,leaves(q))
                if(pass==2) sites(ld-l+2)=sites(ld-l+2)+1_int64
                peer=dir%leaf(key_rank,leaves(q))+1
                grid=dir%leaf(key_grid,leaves(q))
                nnative=(int(dir%leaf(key_index:key_index+2,leaves(q)),int64)-1)*bn
                cell=stencil%site(:,q)
                if(ld<=l) then
                  mode=primitive_sample
                  if(ld<l) mode=l
                  child=cell/2**(l-ld)-nnative+ix_lo
                  call append_query([int(grid,int64),child,int(mode,int64),merge(cell,0*cell,ld<l)],1.d0)
                else
                  nrestrict=nrestrict+1
                  if(pass==2) plan%restriction_center(:,nrestrict)=virtual_center(l,cell)
                  do w=0,1
                    do v=0,1
                      do u=0,1
                        child=2*cell+[u,v,w]
                        volumes(u+1,v+1,w+1)=child_volume(ld,child)
                      end do
                    end do
                  end do
                  total=sum(volumes)
                  if(total<=0.d0) call mpistop('Yin-Yang AMR restriction has nonpositive volume')
                  do w=0,1
                    do v=0,1
                      do u=0,1
                        child=2*cell+[u,v,w]-nnative+ix_lo
                        call append_query([int(grid,int64),child,int(conserved_sample,int64),0_int64,0_int64,0_int64],volumes(u+1,v+1,w+1)/total)
                      end do
                    end do
                  end do
                end if
                if(pass==2) plan%ends(q,nt)=nterm
              end do
            end do
          end do
        end do
      end do
      if(pass==1) then
        allocate(plan%receivers(4,nt),plan%ends(4,nt),plan%weights(4,nt))
        allocate(plan%volume(nterm),plan%slot(nterm),outgoing(query_size,nterm))
        allocate(plan%restriction_center(3,nrestrict))
        call displace(sc,sd)
        cursor=0
      end if
    end do
    call unique_queries(outgoing,plan%slot,sc,sd,nr)
    call MPI_ALLTOALL(sc,1,MPI_INTEGER,rc,1,MPI_INTEGER,MPI_COMM_WORLD,ierrmpi)
    call displace(rc,rd)
    ns=sum(rc)
    allocate(plan%request(query_size,ns))
    call MPI_ALLTOALLV(outgoing,query_size*sc,query_size*sd,MPI_INTEGER8,plan%request, &
      query_size*rc,query_size*rd,MPI_INTEGER8,MPI_COMM_WORLD,ierrmpi)
    plan%recv_count=sc; plan%recv_off=sd
    plan%send_count=rc; plan%send_off=rd
    allocate(plan%send(nw,ns),plan%recv(nw,nr),plan%first(max_blocks),plan%next(ns))
    allocate(plan%handles(2*yy_world_size),plan%statuses(MPI_STATUS_SIZE,2*yy_world_size))
    plan%first=0
    do t=ns,1,-1
      g=int(plan%request(1,t))
      plan%next(t)=plan%first(g)
      plan%first(g)=t
    end do

  contains

    subroutine append_query(query,weight)
      integer(int64),intent(in) :: query(query_size)
      double precision,intent(in) :: weight
      if(int(nterm,int64)*query_size>=huge(nterm)-query_size) call mpistop('Yin-Yang AMR request count exceeds MPI capacity')
      nterm=nterm+1
      if(pass==1) then
        sc(peer)=sc(peer)+1
      else
        offset=sd(peer)+cursor(peer)+1
        cursor(peer)=cursor(peer)+1
        outgoing(:,offset)=query
        plan%slot(nterm)=offset
        plan%volume(nterm)=weight
      end if
    end subroutine

  end subroutine

  subroutine displace(counts,offsets)
    integer,intent(in) :: counts(:)
    integer,intent(out) :: offsets(:)
    integer :: p
    offsets(1)=0
    do p=2,size(counts)
      offsets(p)=offsets(p-1)+counts(p-1)
    end do
  end subroutine

  subroutine unique_queries(queries,slots,counts,offsets,nunique)
    integer(int64),allocatable,intent(inout) :: queries(:,:)
    integer,intent(inout) :: slots(:),counts(:),offsets(:)
    integer,intent(out) :: nunique
    integer,allocatable :: bucket(:),remap(:),owner(:)
    integer(int64),allocatable :: compact(:,:)
    integer :: n,p,k,q,h,found,limit,first,last
    integer(kind=8) :: hash
    limit=2*size(queries,2)+1
    allocate(bucket(limit),remap(size(queries,2)),owner(size(queries,2)),compact(query_size,size(queries,2)))
    nunique=0
    bucket=0
    do p=1,size(counts)
      first=offsets(p)+1
      last=offsets(p)+counts(p)
      counts(p)=0
      do n=first,last
        hash=int(p,int64)
        do k=1,query_size
          hash=modulo(hash*1000003_8+int(queries(k,n),8),int(limit,8))
        end do
        h=int(hash)+1
        do
          found=bucket(h)
          if(found==0) then
            nunique=nunique+1
            compact(:,nunique)=queries(:,n)
            bucket(h)=nunique
            owner(nunique)=p
            counts(p)=counts(p)+1
            remap(n)=nunique
            exit
          end if
          if(owner(found)==p.and.all(compact(:,found)==queries(:,n))) then
            remap(n)=found
            exit
          end if
          h=mod(h,limit)+1
        end do
      end do
    end do
    call displace(counts,offsets)
    slots=remap(slots)
    deallocate(queries)
    allocate(queries(query_size,nunique))
    queries=compact(:,1:nunique)
  end subroutine

  subroutine cache_geometry(plan,s)
    type(mixed_plan),intent(inout) :: plan
    type(state),intent(in) :: s(max_blocks)
    integer :: n,g,j,k,i,t,axis,idx(3),point(3),low(3),high(3),ndonor,slot,level
    double precision :: pos(3),coord(3)
    low=[ixMlo1,ixMlo2,ixMlo3]
    high=[ixMhi1,ixMhi2,ixMhi3]
    ndonor=count(plan%first/=0)
    allocate(plan%source_slot(max_blocks))
    allocate(plan%source_basis(3,3,ixMlo2:ixMhi2,ixMlo3:ixMhi3,ndonor))
    allocate(plan%target_basis(3,3,size(plan%receivers,2)))
    allocate(plan%displacement(3,size(plan%request,2)),plan%stencil_start(3,size(plan%request,2)))
    allocate(plan%span_inv(2,3,size(plan%request,2)))
    slot=0
    do n=1,igridstail
      g=igrids(n)
      if(plan%first(g)==0) cycle
      slot=slot+1
      plan%source_slot(g)=slot
      do k=ixMlo3,ixMhi3
        do j=ixMlo2,ixMhi2
          call yy_glm_basis(s(g)%x(ixMlo1,j,k,2),s(g)%x(ixMlo1,j,k,3),yy_patch,plan%source_basis(:,:,j,k,slot))
        end do
      end do
    end do
    do t=1,size(plan%receivers,2)
      g=plan%receivers(1,t)
      i=plan%receivers(2,t);j=plan%receivers(3,t);k=plan%receivers(4,t)
      call yy_glm_basis(s(g)%x(i,j,k,2),s(g)%x(i,j,k,3),yy_patch,plan%target_basis(:,:,t))
    end do
    do t=1,size(plan%request,2)
      level=int(plan%request(5,t))
      if(level<=0) cycle
      g=int(plan%request(1,t))
      idx=int(plan%request(2:4,t))
      pos=virtual_center(level,plan%request(6:8,t))
      plan%displacement(:,t)=pos-s(g)%x(idx(1),idx(2),idx(3),:)
      do axis=1,3
        if(high(axis)-low(axis)<2) then
          plan%span_inv(:,axis,t)=0.d0
          cycle
        end if
        point=idx
        point(axis)=max(low(axis),min(high(axis)-2,idx(axis)-1))
        plan%stencil_start(axis,t)=point(axis)
        do i=1,3
          coord(i)=s(g)%x(point(1),point(2),point(3),axis)
          point(axis)=point(axis)+1
        end do
        plan%span_inv(:,axis,t)=1.d0/(coord(2:3)-coord(1:2))
      end do
    end do
  end subroutine

  subroutine rotate_cartesian(values,frame)
    double precision,intent(inout) :: values(nw)
    double precision,intent(in) :: frame(3,3)
    integer :: iv,iw
    do iv=1,nvector
      iw=iw_vector(iv)
      values(iw+1:iw+3)=matmul(frame,values(iw+1:iw+3))
    end do
  end subroutine

  subroutine reconstruct(plan,g,t,idx,v)
    type(mixed_plan),intent(in) :: plan
    integer,intent(in) :: g,t,idx(3)
    double precision,intent(inout) :: v(nw)
    double precision :: sample(nw,3),value(nw),a(nw),b(nw),slopes(nw),center(nw)
    integer :: axis,k,point(3),slot
    center=v
    slot=plan%source_slot(g)
    do axis=1,3
      if(plan%span_inv(1,axis,t)==0.d0) cycle
      point=idx
      point(axis)=plan%stencil_start(axis,t)
      do k=1,3
        value=primitive(point(1),point(2),point(3),:)
        call rotate_cartesian(value,plan%source_basis(:,:,point(2),point(3),slot))
        sample(:,k)=value
        point(axis)=point(axis)+1
      end do
      a=(sample(:,2)-sample(:,1))*plan%span_inv(1,axis,t)
      b=(sample(:,3)-sample(:,2))*plan%span_inv(2,axis,t)
      where((a>0.d0.and.b>0.d0).or.(a<0.d0.and.b<0.d0))
        slopes=sign(min(abs(a),abs(b)),a)
      elsewhere
        slopes=0.d0
      end where
      v=v+slopes*plan%displacement(axis,t)
    end do
    if(v(iw_rho)<=0.d0) v(iw_rho)=center(iw_rho)
    if(v(iw_e)<=0.d0) v(iw_e)=center(iw_e)
  end subroutine

  subroutine exchange_plan(s,plan,radial_only)
    type(state),target,intent(inout) :: s(max_blocks)
    type(mixed_plan),intent(inout) :: plan
    logical,intent(in) :: radial_only
    integer :: n,g,t,idx(3),mode,p,nh,tag,ixO^L,i,j,k,iv,iw,q,offset,last,restriction
    double precision :: value(nw),buf(nw),one_cell(1,1,1,nw),one_point(1,1,1,ndim),before
    before=MPI_WTIME()
    do n=1,igridstail
      g=igrids(n)
      if(plan%first(g)==0) cycle
      block=>s(g)
      primitive=s(g)%w
      ixO^L=ixM^LL;
      if(radial_only) then
        ixOmin1=ixGlo1
        ixOmax1=ixGhi1
      end if
      call phys_to_primitive(ixG^LL,ixO^L,primitive,s(g)%x)
      t=plan%first(g)
      do while(t/=0)
        idx=int(plan%request(2:4,t))
        mode=int(plan%request(5,t))
        if(mode==conserved_sample) then
          value=s(g)%w(idx(1),idx(2),idx(3),:)
        else
          value=primitive(idx(1),idx(2),idx(3),:)
        end if
        call rotate_cartesian(value,plan%source_basis(:,:,idx(2),idx(3),plan%source_slot(g)))
        if(mode>0) call reconstruct(plan,g,t,idx,value)
        plan%send(:,t)=value
        t=plan%next(t)
      end do
    end do

    seconds(2)=seconds(2)+MPI_WTIME()-before
    before=MPI_WTIME()
    nh=0
    tag=merge(27133,27132,radial_only)
    do p=1,yy_world_size
      if(plan%recv_count(p)>0) then
        nh=nh+1
        call MPI_IRECV(plan%recv(1,plan%recv_off(p)+1),nw*plan%recv_count(p),MPI_DOUBLE_PRECISION, &
          p-1,tag,MPI_COMM_WORLD,plan%handles(nh),ierrmpi)
      end if
      if(plan%send_count(p)>0) then
        nh=nh+1
        call MPI_ISEND(plan%send(1,plan%send_off(p)+1),nw*plan%send_count(p),MPI_DOUBLE_PRECISION, &
          p-1,tag,MPI_COMM_WORLD,plan%handles(nh),ierrmpi)
      end if
    end do
    if(nh>0) call MPI_WAITALL(nh,plan%handles,plan%statuses,ierrmpi)
    seconds(3)=seconds(3)+MPI_WTIME()-before
    before=MPI_WTIME()

    offset=1
    restriction=0
    do t=1,size(plan%receivers,2)
      g=plan%receivers(1,t)
      i=plan%receivers(2,t); j=plan%receivers(3,t); k=plan%receivers(4,t)
      block=>s(g)
      value=0.d0
      do q=1,4
        last=plan%ends(q,t)
        if(offset==last) then
          buf=plan%recv(:,plan%slot(offset))
        else
          buf=0.d0
          do n=offset,last
            buf=buf+plan%volume(n)*plan%recv(:,plan%slot(n))
          end do
          one_cell(1,1,1,:)=buf
          restriction=restriction+1
          one_point(1,1,1,:)=plan%restriction_center(:,restriction)
          call phys_to_primitive(1,1,1,1,1,1,1,1,1,1,1,1,one_cell,one_point)
          buf=one_cell(1,1,1,:)
        end if
        value=value+plan%weights(q,t)*buf
        offset=last+1
      end do
      do iv=1,nvector
        iw=iw_vector(iv)
        value(iw+1:iw+3)=matmul(transpose(plan%target_basis(:,:,t)),value(iw+1:iw+3))
      end do
      s(g)%w(i,j,k,:)=value
      call phys_to_conserved(ixG^LL,i,j,k,i,j,k,s(g)%w,s(g)%x)
    end do
    seconds(4)=seconds(4)+MPI_WTIME()-before
  end subroutine
  }

end module mod_yinyang_glm_amr
