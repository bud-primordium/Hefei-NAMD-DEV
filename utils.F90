! 工具函数模块
! 本模块提供了各种通用工具函数
! 主要功能包括：
! 1. 随机数生成和初始化
! 2. 插值计算
! 3. 复数运算
! 4. 矩阵操作
! 5. 文件输出控制

module utils
  use prec      ! 提供精度控制
  use constants ! 提供常量定义
#ifdef ENABLEMPI
  use mpi       ! 提供MPI并行计算支持
#endif
  implicit none

  private :: interpolate_c, interpolate_r

  interface interpolate
    procedure :: interpolate_c, interpolate_r
  end interface

  type, private :: mathLib
  contains
    ! private
    ! procedure, private :: interpolate_r, interpolate_c
    ! generic, public :: interpolate => interpolate_r, interpolate_c
  end type
  type(mathLib), protected :: maths

contains

  !! 从系统时钟初始化随机数种子
  !! 代码来源：http://fortranwiki.org/fortran/show/random_seed
  subroutine init_random_seed(salt)
    implicit none
    integer, intent(in) :: salt
    integer :: i, n, clock, ierr
    integer, dimension(:), allocatable :: seed
    real :: r

    call random_seed(size = n)
    allocate(seed(n))

    call system_clock(count=clock)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    seed = seed + (modulo(clock, 65536)+300) * salt
    call random_seed(put = seed)

    deallocate(seed)

#ifdef ENABLEMPI
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
    call random_number(r)
    write(*,'(A,I6,A,F6.3)') '[I] Salt: ', salt, ', First random number: ', r
#ifdef ENABLEMPI
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif

  end subroutine

  !! 实数线性插值
  elemental function interpolate_r(x, a0, a1)
    real(q), intent(in) :: x, a0, a1
    real(q)             :: interpolate_r

    interpolate_r = (1.0_q-x) * a0 + x * a1
  end function
    
  !! 复数线性插值
  elemental function interpolate_c(x, a0, a1)
    real(q), intent(in)    :: x
    complex(q), intent(in) :: a0, a1
    complex(q)             :: interpolate_c

    interpolate_c = (con%uno-x) * a0 + x * a1
  end function

  !! 计算复数的模
  elemental function cmplx_norm(x)
    complex(q), intent(in) :: x
    real(q)                :: cmplx_norm

    cmplx_norm = REAL(CONJG(x) * x, kind=q)
  end function

  !! 从向量创建对角矩阵
  pure function zDIAG(N, V)
    implicit none
    integer, intent(in)   :: N
    complex(q), intent(in), dimension(N) :: V
    complex(q), dimension(N,N) :: zDIAG

    logical, dimension(N,N) :: MASK
    integer :: i

    zDIAG = 0.0_q
    MASK = reshape([(MOD(i,N) == i/N, i=0,N*N-1)], shape=[N,N])
    zDIAG = unpack(V, MASK, zDIAG)
  end function

  !! 支持MPI的消息屏幕输出
  subroutine printToScreen(msg, fh)
    character(len=*), intent(in) :: msg
    integer, intent(in) :: fh
    integer :: ierr
#ifdef ENABLEMPI
    call MPI_FILE_WRITE_SHARED(fh, msg, len(msg), MPI_CHAR, MPI_STATUS_IGNORE, ierr)
    ! call MPI_FILE_WRITE_ORDERED(fh, msg, len(msg), MPI_CHAR, MPI_STATUS_IGNORE, ierr)
#else
    write(fh,*) msg
#endif
  end subroutine

end module utils