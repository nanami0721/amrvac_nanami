program sparse_index
  use, intrinsic :: iso_fortran_env, only: int64
  use mod_global_parameters, only: domain_nx1,domain_nx2,domain_nx3,block_nx1,block_nx2,block_nx3,yy_world_size
  use mod_yinyang_glm_amr, only: leaf_directory,index_directory,covering_leaf,domain_cells
  implicit none
  type(leaf_directory) :: dir
  integer(int64) :: center(3), top(3), coarse(3)
  integer :: bn(3)

  domain_nx1=16
  domain_nx2=8192
  domain_nx3=8192
  block_nx1=8
  block_nx2=8
  block_nx3=8
  yy_world_size=2
  bn=[block_nx1,block_nx2,block_nx3]
  top=domain_cells(20)
  if(top(2)/=4294967296_int64.or.top(3)/=4294967296_int64) error stop 'level-20 domain overflow'
  center=[42_int64,3500000004_int64,3900000004_int64]
  coarse=center/2
  allocate(dir%leaf(7,3))
  dir%leaf(:,1)=[20,int(center/bn)+1,1,5,0]
  dir%leaf(:,2)=[19,int(coarse/bn)+1,0,6,0]
  dir%leaf(:,3)=[20,int((top-1)/bn)+1,1,7,0]
  call index_directory(dir)
  if(size(dir%bucket)/=7) error stop 'directory is not proportional to actual leaves'
  if(covering_leaf(dir,20,center,1)/=1) error stop 'high-index same-level lookup'
  if(covering_leaf(dir,19,coarse,1)/=1) error stop 'high-index fine lookup'
  if(covering_leaf(dir,20,center,0)/=2) error stop 'high-index coarse lookup'
  if(covering_leaf(dir,20,top-1,1)/=3) error stop 'upper endpoint lookup'
  print '(a)', 'Yin-Yang level-20 sparse lookup passed'
end program
