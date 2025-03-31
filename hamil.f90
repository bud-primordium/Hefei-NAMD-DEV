! 哈密顿量模块
! 本模块用于构建和处理分子动力学模拟中的哈密顿量
! 主要功能包括：
! 1. 构建时间依赖的Kohn-Sham哈密顿量
! 2. 处理非绝热耦合矩阵
! 3. 管理波函数和布居数数据
! 4. 提供多种时间演化算法
! 5. 支持多电子态计算

module hamil
  use prec      !! 精度定义模块
  use fileio    !! 文件IO模块
  use utils     !! 工具函数模块
  use couplings !! 耦合计算模块
  use constants !! 常量定义模块
  implicit none

  private :: calcDij_i, calcDij_r

  interface calcDij
    procedure :: calcDij_i, calcDij_r
  end interface

  ! TDKS类型定义
  ! 包含波函数、哈密顿量、布居数等关键变量
  type TDKS
    integer :: ndim
    ! _[c,p,n] means current, previous, next
    complex(kind=q), allocatable, dimension(:) :: psi_c, psi_p, psi_n  !! 当前、前一个、下一个时间步的波函数
    complex(kind=q), allocatable, dimension(:,:) :: psi_a             !! 所有时间步的波函数
    ! the result of hamiltonian acting on a vector
    ! complex(kind=q), allocatable, dimension(:) :: hpsi
    !! population
    real(kind=q), allocatable, dimension(:,:) :: pop_a, mppop_a       !! 布居数
    !!!!real(kind=q), allocatable, dimension(:) :: pop_c, pop_n
    ! real(kind=q) :: norm_c !! unused

    complex(kind=q), allocatable, dimension(:,:) :: ham_c             !! 当前哈密顿量

    ! KS eigenvalues & Non-adiabatic couplings
    real(kind=q), allocatable, dimension(:,:) :: eigKs    !! 这些在长版本中已弃用
    real(kind=q), allocatable, dimension(:,:,:) :: NAcoup

    !! surface hopping related
    real(kind=q), allocatable, dimension(:,:) :: sh_pops, sh_mppops   !! 表面跳跃布居数
    real(kind=q), allocatable, dimension(:,:,:) :: sh_prob            !! 跃迁概率 P_ij (j,i,N-1)

    !! decoherence induced surface hopping
    real(kind=q), allocatable, dimension(:,:) :: dish_pops, dish_mppops !! 退相干诱导的表面跳跃布居数（不再使用？）
    real(kind=q), allocatable, dimension(:,:) :: recom_pops            !! 重组布居数
    real(kind=q), allocatable, dimension(:)   :: dish_decmoment       !! 退相干动量 t_i(t)
    ! whether the memory has been allocated
    logical :: LALLO = .FALSE.                                        !! 内存是否已分配

  end type

contains

  ! initTDKS
  ! 输入：
  ! - ks: TDKS类型变量
  ! 流程：
  ! 1. 分配内存
  ! 2. 初始化波函数
  ! 3. 设置初始布居数
  subroutine initTDKS(ks)
    implicit none

    type(TDKS), intent(inout)  :: ks

    integer :: N

    !! 内存分配
    ks%ndim = inp%NBASIS
    N = inp%NBASIS

    if (.NOT. ks%LALLO) then
      allocate(ks%psi_c(N))
      allocate(ks%psi_p(N))
      allocate(ks%psi_n(N))
      ! allocate(ks%hpsi(N))

      allocate(ks%ham_c(N,N))

      ! allocate(ks%eigKs(N, inp%NSW))
      ! allocate(ks%NAcoup(N, N, inp%NSW))

      select case (inp%ALGO)
      case ('FSSH')
        allocate(ks%psi_a(N, inp%NAMDTIME))
        allocate(ks%pop_a(N, inp%NAMDTIME))
        allocate(ks%sh_pops(N, inp%NAMDTIME))
        allocate(ks%sh_prob(N, N, inp%NAMDTIME-1))
        if (inp%LSPACE) then
          allocate(ks%mppop_a(inp%NBADNS, inp%NAMDTIME))
          allocate(ks%sh_mppops(inp%NBADNS, inp%NAMDTIME))
        end if
      case ('DISH')
        allocate(ks%dish_decmoment(N))
        allocate(ks%dish_pops(1,1))
      case default
      end select

      ! 现在复制 olap%eig&Dij => ks%eig&Dij
      ! ks%eigKs = olap%Eig
      ! ks%NAcoup = olap%Dij

      ks%LALLO = .TRUE.
    end if

    !! 初始化，将载流子放在初始能带上
    ks%psi_c = con%cero
    ks%psi_c(inp%INIBAND) = con%uno

    select case (inp%ALGO)
    case ('FSSH')
      ks%psi_a = con%cero
      ks%psi_a(inp%INIBAND, 1) = con%uno
      ks%pop_a = 0.0_q 
      ks%pop_a(inp%INIBAND, 1) = 1.0_q
      ks%sh_pops = 0.0_q
      ks%sh_prob = 0.0_q
    case ('DISH')
    end select

  end subroutine

  ! make_hamil
  ! 输入：
  ! - tion: 当前时间步
  ! - tele: 当前电子步
  ! - olap: 重叠矩阵
  ! - ks: TDKS类型变量
  ! 流程：
  ! 1. 构建非绝热耦合矩阵
  ! 2. 添加能量本征值
  subroutine make_hamil(tion, tele, olap, ks)
    implicit none

    integer, intent(in) :: tion, tele
    type(overlap), intent(in) :: olap
    type(TDKS), intent(inout) :: ks

    integer :: RTIME,XTIME !! 左右时间步
    integer :: i
    RTIME = tion
    XTIME = RTIME + 1

    ks%ham_c = (0.0_q, 0.0_q)
    if (tele <= (inp%NELM / 2)) then
      ks%ham_c(:,:) = interpolate((tele + inp%NELM/2.0_q - 0.5_q) / inp%NELM, &
                                          olap%Dij(:,:,RTIME-1), olap%Dij(:,:,RTIME))
    else 
      ks%ham_c(:,:) = interpolate((tele - inp%NELM/2.0_q - 0.5_q) / inp%NELM, &
                                          olap%Dij(:,:,RTIME), olap%Dij(:,:,XTIME))
    end if

    !! 乘以 -i * hbar
    ks%ham_c = -con%I * con%hbar * ks%ham_c 
    
    !! 能量本征值部分
    do i=1, ks%ndim
      ks%ham_c(i,i) = interpolate((tele - 0.5_q) / inp%NELM, &
                                  olap%Eig(i,RTIME), olap%Eig(i,XTIME))
    end do

  end subroutine

  ! make_hamil2
  ! 输入：
  ! - tion: 当前时间步
  ! - tele: 当前电子步
  ! - olap: 重叠矩阵
  ! - ks: TDKS类型变量
  ! 流程：
  ! 1. 构建改进的非绝热耦合矩阵
  ! 2. 添加能量本征值
  subroutine make_hamil2(tion, tele, olap, ks)
    implicit none

    integer, intent(in) :: tion, tele
    type(overlap), intent(in) :: olap
    type(TDKS), intent(inout) :: ks

    integer :: RTIME,XTIME !! 左右时间步
    integer :: i
    RTIME = tion
    XTIME = RTIME + 1

    ks%ham_c = (0.0_q, 0.0_q)
    if (tele <= (inp%NELM / 2)) then
      ks%ham_c(:,:) = interpolate((tele + inp%NELM/2.0_q) / inp%NELM, &
                                          olap%Dij(:,:,RTIME-1), olap%Dij(:,:,RTIME))
    else 
      ks%ham_c(:,:) = interpolate((tele - inp%NELM/2.0_q) / inp%NELM, &
                                          olap%Dij(:,:,RTIME), olap%Dij(:,:,XTIME))
    end if

    !! 乘以 -i * hbar
    ks%ham_c = -con%I * con%hbar * ks%ham_c 
    
    !! 能量本征值部分
    do i=1, ks%ndim
      ks%ham_c(i,i) = interpolate(REAL(tele, kind=q) / inp%NELM, &
                                  olap%Eig(i,RTIME), olap%Eig(i,XTIME))
    end do

  end subroutine

  ! make_hamil_wrong
  ! 输入：
  ! - tion: 当前时间步
  ! - tele: 当前电子步
  ! - olap: 重叠矩阵
  ! - ks: TDKS类型变量
  ! 流程：
  ! 1. 构建错误的哈密顿量（用于测试）
  subroutine make_hamil_wrong(tion, tele, olap, ks)
    implicit none

    integer, intent(in) :: tion, tele
    type(overlap), intent(in) :: olap
    type(TDKS), intent(inout) :: ks

    integer :: RTIME,XTIME
    integer :: i
    RTIME = tion
    XTIME = RTIME + 1

    ks%ham_c = (0.0_q, 0.0_q)
    ks%ham_c(:,:) = interpolate(REAL(tele, kind=q) / inp%NELM, &
                                olap%Dij(:,:,RTIME), olap%Dij(:,:,XTIME))

    !! 乘以 -i * hbar
    ks%ham_c = -con%I * con%hbar * ks%ham_c 
    
    !! 能量本征值部分
    do i=1, ks%ndim
      ks%ham_c(i,i) = interpolate(REAL(tele, kind=q) / inp%NELM, &
                                  (olap%Eig(i,RTIME)+olap%Eig(i,XTIME))/2.0_q, &
                                  (olap%Eig(i,XTIME)+olap%Eig(i,XTIME+1))/2.0_q)
    end do

  end subroutine

  ! calcDij_r
  ! 输入：
  ! - dij1, dij2, dij3: 三个时间步的非绝热耦合
  ! - tele: 当前电子步
  ! 输出：
  ! - calcDij_r: 插值后的非绝热耦合
  elemental function calcDij_r(dij1, dij2, dij3, tele)
    implicit none
    real(kind=q), intent(in) :: dij1, dij2, dij3
    real(kind=q), intent(in) :: tele
    real(kind=q) :: calcDij_r

    if (tele <= (inp%NELM / 2)) then
      calcDij_r = interpolate((tele + inp%NELM/2.0_q) / inp%NELM, &
                                          dij1, dij2)
    else 
      calcDij_r = interpolate((tele - inp%NELM/2.0_q) / inp%NELM, &
                                          dij2, dij3)
    end if
  end function

  ! calcDij_i
  ! 输入：
  ! - dij1, dij2, dij3: 三个时间步的非绝热耦合
  ! - tele: 当前电子步
  ! 输出：
  ! - calcDij_i: 插值后的非绝热耦合
  elemental function calcDij_i(dij1, dij2, dij3, tele)
    implicit none
    real(kind=q), intent(in) :: dij1, dij2, dij3
    integer, intent(in) :: tele
    real(kind=q) :: calcDij_i

    calcDij_i = calcDij_r(dij1, dij2, dij3, REAL(tele, kind=q))
  end function

  end module
