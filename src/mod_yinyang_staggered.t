! Geometry-only staggered sampling. No solver state is updated.
module mod_yinyang_staggered
  use mpi
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  integer, parameter, public :: yy_face_location=1, yy_edge_location=2
  public :: yy_staggered_table, yy_staggered_build, yy_staggered_exchange
  public :: yy_staggered_clear, yy_staggered_basis, yy_staggered_offset
  integer,parameter,public :: yy_staggered_max_sources=128
  public :: yy_staggered_sources, yy_staggered_coordinates, yy_staggered_bracket
  abstract interface
    subroutine yy_staggered_sources(family,c,idx,rank,request,weight,level,n)
      integer,intent(in) :: family,c,idx(3)
      integer,intent(out) :: rank(:),request(:,:),level(:),n
      double precision,intent(out) :: weight(:)
    end subroutine
  end interface
  type yy_staggered_table
    integer :: comm=MPI_COMM_NULL, stencil=24
    integer, allocatable :: request(:,:), slot(:,:), sc(:),rc(:),sd(:),rd(:)
    logical, allocatable :: covered(:)
    integer,allocatable :: term_end(:,:),term_slot(:),donor_level_range(:,:)
    double precision,allocatable :: term_weight(:)
    double precision, allocatable :: donor_position(:,:), donor_basis(:,:)
    double precision, allocatable :: coefficient(:,:,:), send(:),recv(:)
  end type
contains
  subroutine fail(comm,code)
    integer,intent(in) :: comm,code
    integer :: ierr,rank
    character(len=96) :: message
    select case(code)
    case(10);message='invalid patch or sampling family'
    case(11,12);message='invalid grid sizes or array shapes'
    case(13);message='invalid radial faces or angular extents'
    case(14,15);message='invalid block directory entry'
    case(16,17);message='duplicate or missing level-one block'
    case(18);message='singular or ill-conditioned Cartesian reconstruction'
    case(19,20);message='donor request routed to the wrong owner or patch'
    case(21);message='exchange arrays disagree with table'
    case(23);message='duplicate local block ID on one rank'
    case(24);message='empty global block directory'
    case(25);message='nonfinite coordinate input'
    case default;message='invalid staggered table'
    end select
    call MPI_COMM_RANK(comm,rank,ierr)
    write(*,'(a,i0,2a)') 'Yin-Yang staggered table rank ',rank,': ',trim(message)
    call MPI_ABORT(comm,code,ierr)
  end subroutine

  pure subroutine yy_staggered_offset(family,c,offset)
    integer,intent(in) :: family,c
    double precision,intent(out) :: offset(3)
    offset=0.5d0
    if(family==yy_face_location) then
      offset(c)=1.d0
    else
      offset=1.d0;offset(c)=0.5d0
    end if
  end subroutine

  pure subroutine yy_staggered_basis(theta,phi,patch,b)
    double precision,intent(in) :: theta,phi
    integer,intent(in) :: patch
    double precision,intent(out) :: b(3,3)
    double precision :: tmp(3)
    b(:,1)=[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]
    b(:,2)=[cos(theta)*cos(phi),cos(theta)*sin(phi),-sin(theta)]
    b(:,3)=[-sin(phi),cos(phi),0.d0]
    if(patch==1) then
      b(1,:)=-b(1,:);tmp=b(2,:);b(2,:)=b(3,:);b(3,:)=tmp
    end if
  end subroutine

  subroutine yy_staggered_clear(t)
    type(yy_staggered_table),intent(inout) :: t
    type(yy_staggered_table) :: empty
    t=empty
  end subroutine

  ! Directory rows: zero-based block coordinates, local block ID, world rank,
  ! patch (0/1). Every rank supplies the same complete level-one directory.
  ! Receiver points are (r,theta,phi) in its own patch. Donors are scalar normal
  ! face components or tangent edge components, NOT dt*length edge integrals.
  ! Request rows: local block ID, three zero-based local cell indices, component.
  ! An index -1 denotes the physical lower face/edge, not an unfilled ghost donor.
  subroutine yy_staggered_build(t,comm,patch,family,bn,nn,radial,angular,directory,points,radial_scalar,edge_radial_extrapolate,face_widths,face_directions,sources)
    type(yy_staggered_table),intent(inout) :: t
    integer,intent(in) :: comm,patch,family,bn(3),nn(3),directory(:,:)
    double precision,intent(in) :: radial(0:),angular(2,2),points(:,:)
    logical,optional,intent(in) :: radial_scalar,edge_radial_extrapolate
    double precision,optional,intent(in) :: face_widths(:,:)
    integer,optional,intent(in) :: face_directions(:)
    procedure(yy_staggered_sources),optional :: sources
    integer,parameter :: max_sources=yy_staggered_max_sources
    integer :: np,rank,ierr,n,p,c,q,a,b,d,h,pass,target,base(3),idx(3),bc(3),loc(3),components
    integer :: qrank(max_sources,24),qreq(5,max_sources,24),qlevel(max_sources,24),qterms(24)
    integer :: owner,nt,ns,nr,off,ib,nb(3),term,total_terms
    double precision :: qweight(max_sources,24)
    integer,allocatable :: owners(:,:,:,:),ids(:,:,:,:),cursor(:),out(:,:),idmap(:)
    double precision :: xyz(3),sph(3),step(2),offset(3),f(3),coord(3),frame(3,3)
    double precision :: w(24),basis(3,24),normal(3,3),inv(3,3),det
    logical :: valid,scalar,extrapolate,moments
    double precision :: design(24,11),target_row(11),fit(24),width(3),scale
    integer :: fit_status
    call yy_staggered_clear(t)
    extrapolate=.false.
    if(present(edge_radial_extrapolate)) extrapolate=edge_radial_extrapolate
    if(extrapolate.and.family/=yy_edge_location) call fail(comm,10)
    scalar=.false.
    if(present(radial_scalar)) scalar=radial_scalar
    moments=present(face_widths).and.present(face_directions)
    if(present(face_widths).neqv.present(face_directions))call fail(comm,12)
    if(moments)then
      if(family/=yy_face_location.or.scalar)call fail(comm,10)
      if(size(face_widths,1)/=3.or.size(face_widths,2)/=size(points,2).or. &
        size(face_directions)/=size(points,2))call fail(comm,12)
      if(any(face_directions<1).or.any(face_directions>3).or.any(face_widths<=0.d0))call fail(comm,12)
      if(.not.all(ieee_is_finite(face_widths)))call fail(comm,25)
    end if
    components=3
    if(scalar) then
      if(family/=yy_face_location) call fail(comm,10)
      components=1;t%stencil=8
    end if
    call MPI_COMM_SIZE(comm,np,ierr);call MPI_COMM_RANK(comm,rank,ierr)
    if(patch<0.or.patch>1.or.(family/=yy_face_location.and.family/=yy_edge_location)) &
      call fail(comm,10)
    if(any(bn<=0).or.any(nn<2)) call fail(comm,11)
    if(any(mod(nn,bn)/=0).or.size(radial)/=nn(1)+1.or.size(points,1)/=3.or.size(directory,1)/=6) &
      call fail(comm,12)
    if(.not.all(ieee_is_finite(radial)).or..not.all(ieee_is_finite(angular)).or. &
      .not.all(ieee_is_finite(points))) call fail(comm,25)
    if(any(radial(1:)<=radial(:nn(1)-1)).or.radial(0)<=0.d0.or.any(angular(:,2)<=angular(:,1))) &
      call fail(comm,13)
    nb=nn/bn;t%comm=comm;nt=size(points,2)
    if(.not.present(sources)) then
      if(size(directory,2)==0) call fail(comm,24)
      if(any(directory(4,:)<1)) call fail(comm,15)
      allocate(owners(0:nb(1)-1,0:nb(2)-1,0:nb(3)-1,0:1),ids(0:nb(1)-1,0:nb(2)-1,0:nb(3)-1,0:1))
      allocate(idmap(maxval(directory(4,:))))
      idmap=0;owners=-1;ids=-1
      do n=1,size(directory,2)
        bc=directory(1:3,n);p=directory(6,n)
        if(any(bc<0).or.any(bc>=nb).or.p<0.or.p>1) call fail(comm,14)
        if(directory(5,n)<0.or.directory(5,n)>=np.or.directory(4,n)<1) call fail(comm,15)
        if(owners(bc(1),bc(2),bc(3),p)/=-1) call fail(comm,16)
        if(directory(5,n)==rank) then
          if(idmap(directory(4,n))/=0) call fail(comm,23)
          idmap(directory(4,n))=n
        end if
        owners(bc(1),bc(2),bc(3),p)=directory(5,n);ids(bc(1),bc(2),bc(3),p)=directory(4,n)
      end do
      if(any(owners<0)) call fail(comm,17)
    end if
    if(present(sources)) then
      allocate(t%term_end(0:t%stencil,nt),t%donor_level_range(2,nt))
      t%term_end=0
      t%donor_level_range=0
    end if
    allocate(t%sc(np),t%rc(np),t%sd(np),t%rd(np),cursor(np),t%covered(nt),t%coefficient(3,t%stencil,nt))
    if(.not.present(sources)) allocate(t%slot(t%stencil,nt),source=0)
    t%sc=0;t%coefficient=0.d0;t%covered=.false.
    step=(angular(:,2)-angular(:,1))/nn(2:3)
    do pass=1,2
      total_terms=0
      do target=1,nt
        call yy_staggered_coordinates(points(:,target),patch,sph,xyz)
        valid=.true.;q=0
        if(moments)then
          scale=max(face_widths(1,target),points(1,target)*maxval(step))
          call face_row(patch,face_directions(target),points(:,target),face_widths(:,target),xyz,scale,target_row)
        end if
        do c=1,components
          call yy_staggered_offset(family,c,offset)
          call yy_staggered_bracket(sph,offset,nn,radial,angular,step,extrapolate.and.c==1,moments,base,f,valid)
          if(.not.valid) exit
          do h=0,1
            do b=0,1
              do a=0,1
                q=q+1;idx=base+[a,b,h];bc=max(0,idx/bn);loc=idx-bc*bn
                if(present(sources)) then
                  call sources(family,c,idx,qrank(:,q),qreq(:,:,q),qweight(:,q),qlevel(:,q),qterms(q))
                  if(qterms(q)==0) valid=.false.
                else
                  owner=owners(bc(1),bc(2),bc(3),1-patch)
                  if(owner<0) valid=.false.
                  qrank(1,q)=owner+1;qreq(:,1,q)=[ids(bc(1),bc(2),bc(3),1-patch),loc,c]
                  qterms(q)=1
                end if
                w(q)=merge(1.d0-f(1),f(1),a==0)*merge(1.d0-f(2),f(2),b==0)*merge(1.d0-f(3),f(3),h==0)
                if(present(sources).and..not.moments) then
                  if(abs(w(q))<=32.d0*epsilon(1.d0)) qterms(q)=0
                end if
                coord=[radial_sample(radial,idx(1),offset(1)),angular(:,1)+(idx(2:3)+offset(2:3))*step]
                call yy_staggered_basis(coord(2),coord(3),1-patch,frame)
                basis(:,q)=frame(:,c)
                if(moments)then
                  width=[1.d0,step]
                  if(c/=1)width(1)=radial(idx(1)+1)-radial(idx(1))
                  call face_row(1-patch,c,coord,width,xyz,scale,design(q,:))
                end if
              end do
            end do
          end do
        end do
        if(.not.valid) cycle
        t%covered(target)=.true.
        if(pass==1) then
          do q=1,t%stencil
            do term=1,qterms(q)
              t%sc(qrank(term,q))=t%sc(qrank(term,q))+1
            end do
          end do
        else
          ! Weighted Cartesian constant fit. Normalization makes any uniform
          ! Cartesian vector exact despite non-collocated component samples.
          if(moments)then
            ! Keep locality of the trilinear stencil without zero-weight rank loss.
            ! A small positive floor supplies off-plane derivative information.
            w=sqrt(max(w,1.d-2))
            do q=1,t%stencil
              design(q,:)=w(q)*design(q,:)
            end do
            call face_fit(design,target_row,fit,fit_status)
            fit=fit*w
            if(fit_status/=0)call fail(comm,18)
          else if(.not.scalar) then
          normal=0.d0
          do q=1,t%stencil
            do d=1,3
              normal(:,d)=normal(:,d)+w(q)*basis(:,q)*basis(d,q)
            end do
          end do
          inv(:,1)=cross(normal(:,2),normal(:,3))
          inv(:,2)=cross(normal(:,3),normal(:,1))
          inv(:,3)=cross(normal(:,1),normal(:,2))
          det=dot_product(normal(:,1),inv(:,1))
          if(det<1.d-8) call fail(comm,18)
          inv=transpose(inv)/det
          end if
          if(present(sources)) then
            t%term_end(0,target)=total_terms
            t%donor_level_range(:,target)=[huge(1),0]
          end if
          do q=1,t%stencil
            do term=1,qterms(q)
              p=qrank(term,q);off=t%sd(p)+cursor(p)+1;cursor(p)=cursor(p)+1
              out(:,off)=qreq(:,term,q)
              if(present(sources)) then
                total_terms=total_terms+1
                t%term_slot(total_terms)=off
                t%term_weight(total_terms)=qweight(term,q)
                t%donor_level_range(1,target)=min(t%donor_level_range(1,target),qlevel(term,q))
                t%donor_level_range(2,target)=max(t%donor_level_range(2,target),qlevel(term,q))
              else
                t%slot(q,target)=off
              end if
            end do
            if(present(sources)) t%term_end(q,target)=total_terms
            if(moments)then
              call yy_staggered_basis(points(2,target),points(3,target),patch,frame)
              t%coefficient(:,q,target)=fit(q)*frame(:,face_directions(target))
            else if(scalar) then
              ! Br is a scalar under patch rotation at a fixed physical point.
              ! Its own face lattice reaches radial walls without extrapolation.
              t%coefficient(:,q,target)=w(q)*xyz/points(1,target)
            else
              t%coefficient(:,q,target)=w(q)*matmul(inv,basis(:,q))
            end if
          end do
        end if
      end do
      if(pass==1) then
        t%sd(1)=0
        do p=2,np
          t%sd(p)=t%sd(p-1)+t%sc(p-1)
        end do
        allocate(out(5,sum(t%sc)));cursor=0
        if(present(sources)) allocate(t%term_slot(sum(t%sc)),t%term_weight(sum(t%sc)))
      end if
    end do
    call compact_requests(t,out)
    call MPI_ALLTOALL(t%sc,1,MPI_INTEGER,t%rc,1,MPI_INTEGER,comm,ierr)
    t%rd(1)=0
    do p=2,np
      t%rd(p)=t%rd(p-1)+t%rc(p-1)
    end do
    ns=sum(t%rc);nr=sum(t%sc)
    allocate(t%request(5,ns),t%send(ns),t%recv(nr))
    call MPI_ALLTOALLV(out,5*t%sc,5*t%sd,MPI_INTEGER,t%request,5*t%rc,5*t%rd,MPI_INTEGER,comm,ierr)
    if(present(sources)) return
    allocate(t%donor_position(3,ns),t%donor_basis(3,ns))
    do n=1,ns
      ib=t%request(1,n);owner=idmap(ib)
      if(owner==0) call fail(comm,19)
      if(directory(6,owner)/=patch) call fail(comm,20)
      c=t%request(5,n);idx=directory(1:3,owner)*bn+t%request(2:4,n)
      call yy_staggered_offset(family,c,offset)
      coord=[radial_sample(radial,idx(1),offset(1)),angular(:,1)+(idx(2:3)+offset(2:3))*step]
      call yy_staggered_basis(coord(2),coord(3),patch,frame)
      t%donor_position(:,n)=coord(1)*frame(:,1);t%donor_basis(:,n)=frame(:,c)
    end do
  end subroutine

  pure subroutine yy_staggered_coordinates(point,patch,sph,xyz)
    double precision,intent(in) :: point(:)
    integer,intent(in) :: patch
    double precision,intent(out) :: sph(3),xyz(3)
    double precision :: frame(3,3),v(3)
    call yy_staggered_basis(point(2),point(3),patch,frame)
    xyz=point(1)*frame(:,1)
    v=xyz
    if(patch==0) v=[-xyz(1),xyz(3),xyz(2)]
    sph=[point(1),atan2(sqrt(v(1)**2+v(2)**2),v(3)),atan2(v(2),v(1))]
  end subroutine

  pure subroutine yy_staggered_bracket(sph,offset,nn,radial,angular,step,extrapolate,moments,base,f,valid)
    double precision,intent(in) :: sph(3),offset(3),radial(0:),angular(2,2),step(2)
    integer,intent(in) :: nn(3)
    logical,intent(in) :: extrapolate,moments
    integer,intent(out) :: base(3)
    double precision,intent(out) :: f(3)
    logical,intent(out) :: valid
    double precision :: left,right,x,tolerance,u
    integer :: lo,hi,mid,d
    valid=.false.
    ! Bracket the shared nonuniform radial mesh without extrapolation.
    lo=0;hi=nn(1)-1
    if(offset(1)==1.d0) lo=-1
    left=radial_sample(radial,lo,offset(1));right=radial_sample(radial,hi,offset(1))
    x=sph(1)
    tolerance=1.d-12*max(1.d0,abs(left),abs(right))
    ! Optional one-sided linear closure only for radially centered E_r.
    ! It uses the first/last two physical samples, never ghost donors.
    if(extrapolate.or.moments) then
      if(x<radial(0)-tolerance.or.x>radial(nn(1))+tolerance) then
        return
      end if
      x=max(radial(0),min(radial(nn(1)),x))
    else
    if(x<left-tolerance.or.x>right+tolerance) then
        return
      end if
      x=max(left,min(right,x))
    end if
    do while(hi-lo>1)
      mid=(lo+hi)/2
      if(radial_sample(radial,mid,offset(1))<=x) then
        lo=mid
      else
        hi=mid
      end if
    end do
    ! A positive-weight fit uses both sides even at a lattice node.
    ! Canonicalize ties so roundoff cannot select a different stencil.
    if(moments)then
      if(abs(x-radial_sample(radial,lo+1,offset(1)))<=tolerance)then
        x=radial_sample(radial,lo+1,offset(1));lo=min(lo+1,nn(1)-2)
      else if(abs(x-radial_sample(radial,lo,offset(1)))<=tolerance)then
        x=radial_sample(radial,lo,offset(1))
      end if
    end if
    base(1)=lo;left=radial_sample(radial,lo,offset(1));right=radial_sample(radial,lo+1,offset(1))
    f(1)=(x-left)/(right-left)
    do d=2,3
      u=(sph(d)-angular(d-1,1))/step(d-1)-offset(d)
      if(moments.and.abs(u-dble(nint(u)))<1.d-12)u=dble(nint(u))
      lo=0;if(offset(d)==1.d0) lo=-1
      if(u<lo-1.d-12.or.u>nn(d)-1+1.d-12) then
        return
      end if
      u=max(dble(lo),min(dble(nn(d)-1),u))
      base(d)=min(floor(u),nn(d)-2);f(d)=u-base(d)
    end do
    valid=.true.
  end subroutine

  pure double precision function radial_sample(radial,i,fraction) result(r)
    double precision,intent(in) :: radial(0:),fraction
    integer,intent(in) :: i
    if(fraction==1.d0) then
      r=radial(i+1)
    else
      r=0.5d0*(radial(i)+radial(i+1))
    end if
  end function

  ! A native face/edge appears in many virtual stencils. Compact per-peer
  ! request lists in first-use order while preserving every reconstruction sum.
  subroutine compact_requests(t,out)
    type(yy_staggered_table),intent(inout) :: t
    integer,intent(inout) :: out(:,:)
    integer,allocatable :: bucket(:),mapping(:)
    integer :: capacity,peer,start,stop,nunique,n,q,h,found
    integer(kind=8) :: key
    if(maxval(t%sc)>(huge(capacity)-1)/2) call fail(t%comm,12)
    capacity=2*maxval(t%sc)+1
    allocate(bucket(capacity),mapping(size(out,2)))
    nunique=0
    do peer=1,size(t%sc)
      start=t%sd(peer)+1
      stop=t%sd(peer)+t%sc(peer)
      t%sd(peer)=nunique
      if(stop<start) cycle
      bucket=0
      do n=start,stop
        key=0_8
        do q=1,size(out,1)
          key=modulo(key*1000003_8+int(out(q,n),8),int(capacity,8))
        end do
        h=int(key)+1
        do
          found=bucket(h)
          if(found==0) then
            nunique=nunique+1
            out(:,nunique)=out(:,n)
            bucket(h)=nunique
            found=nunique
            exit
          end if
          if(all(out(:,found)==out(:,n))) exit
          h=modulo(h,capacity)+1
        end do
        mapping(n)=found
      end do
      t%sc(peer)=nunique-t%sd(peer)
    end do
    if(allocated(t%term_slot)) then
      do n=1,size(t%term_slot)
        t%term_slot(n)=mapping(t%term_slot(n))
      end do
    else
      do n=1,size(t%covered)
        if(.not.t%covered(n)) cycle
        do q=1,t%stencil
          t%slot(q,n)=mapping(t%slot(q,n))
        end do
      end do
    end if
  end subroutine

  ! Call once per rank outside block loops. Input samples follow request order;
  ! output is Cartesian. Uncovered receivers are zero and MUST be masked.
  subroutine yy_staggered_exchange(t,samples,cartesian)
    type(yy_staggered_table),intent(inout) :: t
    double precision,intent(in) :: samples(:)
    double precision,intent(out) :: cartesian(:,:)
    integer :: ierr,n,q,term
    double precision :: value
    if(size(samples)/=size(t%send).or.size(cartesian,1)/=3.or.size(cartesian,2)/=size(t%covered)) &
      call fail(t%comm,21)
    t%send=samples
    call MPI_ALLTOALLV(t%send,t%rc,t%rd,MPI_DOUBLE_PRECISION,t%recv,t%sc,t%sd,MPI_DOUBLE_PRECISION,t%comm,ierr)
    cartesian=0.d0
    if(allocated(t%term_end)) then
      do n=1,size(t%covered)
        if(.not.t%covered(n)) cycle
        do q=1,t%stencil
          value=0.d0
          do term=t%term_end(q-1,n)+1,t%term_end(q,n)
            value=value+t%term_weight(term)*t%recv(t%term_slot(term))
          end do
          cartesian(:,n)=cartesian(:,n)+t%coefficient(:,q,n)*value
        end do
      end do
    else
      do n=1,size(t%covered)
        if(.not.t%covered(n)) cycle
        do q=1,t%stencil
          cartesian(:,n)=cartesian(:,n)+t%coefficient(:,q,n)*t%recv(t%slot(q,n))
        end do
      end do
    end if
  end subroutine

  ! Face-average moments of B(x)=b+G*(x-origin)/scale with trace(G)=0.
  ! Four-point Gauss integration in each surface coordinate; normalized by the
  ! analytic native spherical face area. Used only during table construction.
  pure subroutine face_row(patch,c,center,width,origin,scale,row)
    integer,intent(in) :: patch,c
    double precision,intent(in) :: center(:),width(:),origin(:),scale
    double precision,intent(out) :: row(:)
    double precision,parameter :: nodes(4)=[-.861136311594052575d0,-.339981043584856265d0, &
      .339981043584856265d0,.861136311594052575d0]
    double precision,parameter :: weights(4)=[.347854845137453857d0,.652145154862546143d0, &
      .652145154862546143d0,.347854845137453857d0]
    double precision :: x(3),frame(3,3),offset(3),normal(3),mom(3,3),jac,area,w
    integer :: a,b,d,e,j,k,n
    d=mod(c,3)+1;e=mod(c+1,3)+1;row=0.d0;mom=0.d0
    do b=1,4
      do a=1,4
        x=center;x(d)=x(d)+.5d0*width(d)*nodes(a);x(e)=x(e)+.5d0*width(e)*nodes(b)
        call yy_staggered_basis(x(2),x(3),patch,frame)
        offset=(x(1)*frame(:,1)-origin)/scale;normal=frame(:,c)
        select case(c)
        case(1);jac=x(1)**2*sin(x(2))
        case(2);jac=x(1)*sin(x(2))
        case(3);jac=x(1)
        end select
        w=.25d0*width(d)*width(e)*weights(a)*weights(b)*jac
        row(1:3)=row(1:3)+w*normal
        do k=1,3
          mom(:,k)=mom(:,k)+w*normal*offset(k)
        end do
      end do
    end do
    select case(c)
    case(1);area=center(1)**2*2.d0*sin(center(2))*sin(.5d0*width(2))*width(3)
    case(2);area=center(1)*width(1)*sin(center(2))*width(3)
    case(3);area=center(1)*width(1)*width(2)
    end select
    n=3
    do k=1,3
      do j=1,3
        if(j==3.and.k==3)cycle
        n=n+1;row(n)=mom(j,k)
        if(j==k)row(n)=row(n)-mom(3,3)
      end do
    end do
    row=row/area
  end subroutine

  ! Reorthogonalized QR avoids squaring the fit condition number. All 24 face
  ! samples have positive weights, including when a target lies on a lattice.
  pure subroutine face_fit(design,target,coefficient,status)
    double precision,intent(in) :: design(24,11),target(11)
    double precision,intent(out) :: coefficient(24)
    integer,intent(out) :: status
    double precision :: q(24,11),r(11,11),z(11),projection
    integer :: i,j,pass
    q=design;r=0.d0;z=0.d0;coefficient=0.d0;status=0
    do j=1,11
      do pass=1,2
        do i=1,j-1
          projection=dot_product(q(:,i),q(:,j));r(i,j)=r(i,j)+projection
          q(:,j)=q(:,j)-projection*q(:,i)
        end do
      end do
      r(j,j)=sqrt(sum(q(:,j)**2))
      if(r(j,j)<1.d-10)then
        status=1;return
      end if
      q(:,j)=q(:,j)/r(j,j)
    end do
    do j=1,11
      z(j)=(target(j)-dot_product(r(1:j-1,j),z(1:j-1)))/r(j,j)
    end do
    coefficient=matmul(q,z)
    if(maxval(abs(matmul(transpose(design),coefficient)-target))>1.d-10)status=2
  end subroutine

  pure function cross(a,b) result(c)
    double precision,intent(in) :: a(3),b(3)
    double precision :: c(3)
    c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
  end function
end module
