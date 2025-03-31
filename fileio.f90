! 文件IO模块
! 本模块用于处理分子动力学模拟中的文件输入输出操作
! 主要功能包括：
! 1. 读取和处理输入参数
! 2. 管理模拟配置信息
! 3. 提供文件操作接口
! 4. 支持并行计算

module fileio
  use omp_lib    ! 提供OpenMP并行计算支持
  use prec       ! 提供精度控制

  implicit none

  ! 定义NAMD信息类型，用于存储模拟参数和配置
  type, private :: namdInfo
    logical, private :: isCreated = .false.
    integer :: NPROG !! TRAJPAR - 轨迹并行数
    integer :: IPROG !! MPI组内进程索引，从0开始
    integer :: NPAR  ! 并行计算参数
    integer :: COLOR ! MPI进程组颜色
    integer :: COMMUNICATOR ! MPI通信器

    integer :: BMIN   ! 最低能带索引
    integer :: BMAX   ! 最高能带索引
    integer :: NBADNS ! 能带数量
    integer :: HORB   ! 最高能量轨道，通常是N(e) 或 1(h)
    integer :: LORB   ! 最低能量轨道，通常是1(e) 或 N(h)
    integer, allocatable, dimension(:,:) :: BASIS !! 活性空间基态宏状态
    integer :: NBASIS ! 绝热态基组数量
    integer :: INIBAND ! 激发电子/空穴的初始绝热态
    integer :: NAMDTINI ! NAMD初始时间步
    integer :: NAMDTIME ! NAMD时间步数
    integer, allocatable, dimension(:), private :: NAMDTINI_A ! NAMD初始时间步数组
    integer, allocatable, dimension(:), private :: INIBAND_A  ! 初始绝热态数组

    integer :: NSW      ! MD步数
    integer :: NTRAJ    ! 表面跳跃轨迹数
    integer :: NELM     ! 电子波传播步数
    integer :: NSAMPLE  ! 采样步数

    integer :: NACBASIS ! 非绝热耦合基组数
    integer :: NACELE   ! 非绝热耦合电子数
    real(kind=q) :: POTIM ! MD时间步长
    real(kind=q) :: TEMP  ! MD温度

    logical :: LHOLE    ! 空穴或电子表面跳跃
    logical :: LSHP     ! 是否进行表面跳跃, 当下值为.TRUE.
    character(len=256) :: ALGO !! 表面跳跃算法
    integer :: ALGO_INT ! 算法整数标识
    logical :: LBINOUT  ! 是否使用二进制输出
    logical :: LCPTXT   ! 是否使用文本格式耦合
    logical :: LSPACE   ! 是否使用活性空间
    ! 运行目录
    character(len=256) :: RUNDIR  ! 运行目录
    character(len=256) :: TBINIT  ! 初始条件文件
    character(len=256) :: DIINIT  !! 纯退相干时间矩阵输入文件
    ! DISH参数
    real(kind=q), allocatable, dimension(:,:) :: DEPHMATR !! 任意两个绝热态之间的纯退相干率，来自DEPHTIME文件，单位为1/fs，对称矩阵
    integer :: NPARDISH ! DISH并行参数
    integer :: DEBUGLEVEL ! 调试级别
  contains
    procedure :: getInstance
    procedure :: getUserInp, printUserInp
    procedure :: setIni, setHLORB, setMPI
    procedure, private :: checkUserInp
  end type

  type(namdInfo), public :: inp

contains

  ! 获取NAMD信息实例
  ! 输入：
  ! - this: 输入输出，NAMD信息对象
  ! 流程：
  ! 1. 检查是否已创建
  ! 2. 设置创建标志
  subroutine getInstance(this)
    class(namdInfo), intent(inout) :: this
    if (this%isCreated) then
      return
    else
      this%isCreated = .TRUE.
    end if
  end subroutine

  ! 设置初始条件
  ! 输入：
  ! - this: 输入输出，NAMD信息对象
  ! - value: 输入，初始条件索引
  ! 流程：
  ! 1. 设置初始能带
  ! 2. 设置初始时间步
  subroutine setIni(this, value)
    class(namdInfo), intent(inout) :: this
    integer, intent(in) :: value
    this%INIBAND = this%INIBAND_A(value)
    this%NAMDTINI = this%NAMDTINI_A(value)
  end subroutine

  ! 设置最高和最低能量轨道
  ! 输入：
  ! - this: 输入输出，NAMD信息对象
  ! - horb: 输入，最高能量轨道
  ! - lorb: 输入，最低能量轨道
  ! 流程：
  ! 1. 设置最高能量轨道
  ! 2. 设置最低能量轨道
  subroutine setHLORB(this, horb, lorb)
    class(namdInfo), intent(inout) :: this
    integer, intent(in) :: horb, lorb
    this%HORB = horb
    this%LORB = lorb
  end subroutine

  ! 设置MPI参数
  ! 输入：
  ! - this: 输入输出，NAMD信息对象
  ! - iprog: 输入，进程索引
  ! - nprog: 输入，总进程数
  ! - color: 输入，进程组颜色
  ! - communicator: 输入，通信器
  ! 流程：
  ! 1. 设置进程参数
  ! 2. 设置通信参数
  subroutine setMPI(this, iprog, nprog, color, communicator)
    implicit none
    class(namdInfo), intent(inout) :: this
    integer, intent(in) :: iprog, nprog, color, communicator

    this%IPROG = iprog
    this%NPROG = nprog
    this%COLOR = color
    this%COMMUNICATOR = communicator
  end subroutine

  ! 读取用户输入参数
  ! 输入：
  ! - this: 输入输出，NAMD信息对象
  ! - nprog: 输入，总进程数
  ! 流程：
  ! 1. 定义本地变量
  ! 2. 读取输入文件
  ! 3. 设置调试级别
  ! 4. 读取初始条件
  ! 5. 读取活性空间
  ! 6. 读取退相干时间
  ! 7. 分配参数
  subroutine getUserInp(this, nprog)
    implicit none

    class(namdInfo), intent(inout) :: this
    integer, intent(in) :: nprog

    ! 本地变量，与"inp"中的变量同名
    integer :: npar      = 1

    integer :: bmin
    integer :: bmax

    integer :: nsw
    integer :: ntraj     = 1000
    integer :: nelm      = 10
    integer :: nsample   

    integer :: nacbasis  = 100
    integer :: nacele    = 1

    integer :: namdtime  = 0
    real(kind=q) :: potim= 1.0_q
    real(kind=q) :: temp = 300_q

    ! 空穴或电子表面跳跃
    logical :: lhole          = .FALSE.
    logical :: lshp           = .TRUE.  !! 是否进行表面跳跃
    character(len=16) :: algo = 'DISH'
    integer :: algo_int       = 0
    logical :: lcptxt         = .TRUE.
    logical :: lspace         = .FALSE.
    logical :: lbinout        = .FALSE.
    ! 运行目录
    character(len=256) :: rundir = 'run'
    character(len=256) :: tbinit = 'INICON'
    character(len=256) :: diinit = 'DEPHTIME'
    character(len=256) :: spinit = 'ACSPACE'

    character(len=256) :: debuglevel = 'I'

    namelist /NAMDPARA/ bmin, bmax,                          &
                          nsample, ntraj, nsw, nelm,           &
                          temp, namdtime, potim,               &
                          lhole, lshp, algo, algo_int, lcptxt, &
                          lspace, nacbasis, nacele,            &
                          npar, lbinout,                   &
                          rundir, tbinit, diinit, spinit,      &
                          debuglevel

    integer :: ierr, i, j, nthread = 1
    logical :: lext

    ! 打开输入文件
    open(file="inp", unit=8, status='unknown', action='read', iostat=ierr)
    if ( ierr /= 0 ) then
      write(*,*) "[E] IOError: I/O error with input file: 'inp'"
      stop
    end if

    ! 读取输入文件
    read(unit=8, nml=NAMDPARA)
    close(unit=8)

    ! 设置OpenMP线程数
    ! nthread = MIN(npardish, omp_get_max_threads())
    ! call omp_set_num_threads(nthread)

    ! 设置调试级别
    select case (debuglevel)
    case ('D')
      this%DEBUGLEVEL = 0
    case ('I')
      this%DEBUGLEVEL = 1
    case ('W')
      this%DEBUGLEVEL = 2
    case ('E')
      this%DEBUGLEVEL = 3
    end select

    ! 读取初始条件
    allocate(this%INIBAND_A(nsample), this%NAMDTINI_A(nsample))
    inquire(file=tbinit, exist=lext)
    if (.NOT. lext) then
      write(*,*) "[E] IOError: File containing initial conditions does NOT exist!"
      stop
    else
      open(unit=9, file=tbinit, action='read')
      do i=1, nsample
        read(unit=9,fmt=*) this%NAMDTINI_A(i), this%INIBAND_A(i)
      end do
      close(9)
      ! 假设BMIN足够大，已弃用
      if (this%INIBAND_A(1) >= bmin) then
        this%INIBAND_A = this%INIBAND_A - bmin + 1
      end if
    end if

    ! 读取活性空间
    if (lspace) then
      allocate(this%BASIS(nacele, nacbasis))
      inquire(file=spinit, exist=lext)
      if (.NOT. lext) then
        write(*,*) "[E] IOError: File containing Active Space does NOT exist!"
        stop
      else
        open(unit=11, file=spinit, action='read')
        do i = 1, nacbasis
          read(unit=11, fmt=*) (this%BASIS(j,i), j=1,nacele)
        end do
        close(11)
        this%NBASIS = nacbasis
      end if
    else
      allocate(this%BASIS(bmax - bmin + 1, 1))
      this%BASIS = reshape([(i, i=bmin, bmax)], shape=[1,bmax-bmin+1])
      this%NBASIS = bmax - bmin + 1
    end if

    ! 读取退相干时间
    if (algo == 'DISH') then
      allocate(this%DEPHMATR(this%NBASIS, this%NBASIS))   !! 读取纯退相干时间矩阵，对角元素为零
      inquire(file=diinit, exist=lext)
      if (.NOT. lext) then
        write(*,*) "[E] IOError: File containing initial conditions of DISH does NOT exist!"
        stop
      else
        open(unit=10, file=diinit, action='read')
        read(unit=10, fmt=*) ((this%DEPHMATR(i,j), j=1, this%NBASIS), i=1, this%NBASIS)
        do i = 1, this%NBASIS
          do j = 1, this%NBASIS
            if (i /= j) then
              this%DEPHMATR(j,i) = 1.0_q / this%DEPHMATR(j,i)
            end if
          end do
        end do
        close(10)
      end if
    end if

    ! 分配参数
    this%NPAR     = MIN(npar, nprog, nsample)

    this%BMIN     = bmin
    this%BMAX     = bmax
    this%NBADNS   = bmax - bmin + 1

    this%NSW      = nsw
    this%NTRAJ    = ntraj
    this%NELM     = nelm
    this%NSAMPLE  = nsample

    this%NAMDTIME = namdtime 
    this%POTIM    = potim
    this%TEMP     = temp

    this%LHOLE    = lhole
    this%LSHP     = lshp
    this%LCPTXT   = lcptxt
    this%ALGO     = algo
    this%ALGO_INT = algo_int
    this%LBINOUT  = lbinout

    this%RUNDIR   = trim(rundir)
    this%TBINIT   = trim(tbinit)
    this%DIINIT   = trim(diinit)

    this%LSPACE   = lspace
    this%NACBASIS = nacbasis
    this%NACELE   = nacele

    call this%checkUserInp()
  end subroutine

  ! 检查用户输入参数
  ! 输入：
  ! - this: 输入，NAMD信息对象
  ! 流程：
  ! 1. 检查能带范围
  ! 2. 检查FSSH算法参数
  subroutine checkUserInp(this)
    implicit none
    class(namdInfo), intent(in) :: this

    integer :: i 

    ! 进行参数检查...
    ! 以下检查将在未来版本中实现
    if (this%BMIN <= 0 .OR. this%BMAX <= 0 .OR. this%BMIN >= this%BMAX) then
      write(*,*) "[E] Please specify the correct BMIN/BMAX"
      stop
    end if

    ! 检查FSSH算法参数
    if (this%ALGO == 'FSSH') then
      do i=1, this%NSAMPLE
        if (this%NAMDTINI_A(i) + this%NAMDTIME - 1 > this%NSW) then
          write(*,*) "[E] NAMDTIME too long..."
          stop
        else if (this%NAMDTINI_A(i) == 1) then
          write(*,*) "[E] NAMDTINI should > 1..."
          stop
        end if
      end do
    end if
  end subroutine

  ! 打印用户输入参数
  ! 输入：
  ! - this: 输入，NAMD信息对象
  ! 流程：
  ! 1. 格式化输出参数
  ! 2. 打印参数值
  subroutine printUserInp(this)
    implicit none
    class(namdInfo), intent(in) :: this
    character(len=4096) :: buf

    write(buf,'(A60,A)') "------------------------------------------------------------", new_line('')
    write(buf,'(A,A30,A3,I10,A)')   trim(buf), 'NPAR',     ' = ', this%NPAR,             new_line('')
    write(buf,'(A,A30,A3,I10,A)')   trim(buf), 'TRAJPAR',  ' = ', this%NPROG,            new_line('')
    write(buf,'(A,A30,A3,I10,A)')   trim(buf), 'MPIGROUP', ' = ', this%COLOR,            new_line('')
    write(buf,'(A,A)')              trim(buf),                                           new_line('')
    write(buf,'(A,A30,A3,I10,A)')   trim(buf), 'BMIN',     ' = ', this%BMIN,             new_line('')
    write(buf,'(A,A30,A3,I10,A)')   trim(buf), 'BMAX',     ' = ', this%BMAX,             new_line('')
    write(buf,'(A,A30,A3,I10,A)')   trim(buf), 'INIBAND',  ' = ', this%INIBAND,          new_line('')
    write(buf,'(A,A)')              trim(buf),                                           new_line('')
    write(buf,'(A,A30,A3,F10.1,A)') trim(buf), 'POTIM',    ' = ', this%POTIM,            new_line('')
    write(buf,'(A,A30,A3,F10.1,A)') trim(buf), 'TEMP',     ' = ', this%TEMP,             new_line('')
    write(buf,'(A,A30,A3,I10,A)')   trim(buf), 'NAMDTINI', ' = ', this%NAMDTINI,         new_line('')
    write(buf,'(A,A30,A3,I10,A)')   trim(buf), 'NAMDTIME', ' = ', this%NAMDTIME,         new_line('')
    write(buf,'(A,A)')              trim(buf) ,                                          new_line('')
    write(buf,'(A,A30,A3,I10,A)')   trim(buf), 'NSW',      ' = ', this%NSW,              new_line('')
    write(buf,'(A,A30,A3,I10,A)')   trim(buf), 'NTRAJ',    ' = ', this%NTRAJ,            new_line('')
    write(buf,'(A,A30,A3,I10,A)')   trim(buf), 'NELM',     ' = ', this%NELM,             new_line('')
    write(buf,'(A,A)')              trim(buf),                                           new_line('')
    write(buf,'(A,A30,A3,L10,A)')   trim(buf), 'LHOLE',    ' = ', this%LHOLE,            new_line('')
    write(buf,'(A,A30,A3,L10,A)')   trim(buf), 'LSHP',     ' = ', this%LSHP,             new_line('')
    write(buf,'(A,A30,A3,A10,A)')   trim(buf), 'ALGO',     ' = ', TRIM(this%ALGO),       new_line('')
    write(buf,'(A,A30,A3,I10,A)')   trim(buf), 'ALGO_INT', ' = ', this%ALGO_INT,         new_line('')
    write(buf,'(A,A30,A3,L10,A)')   trim(buf), 'LBINOUT',  ' = ', this%LBINOUT,          new_line('')
    write(buf,'(A,A30,A3,L10,A)')   trim(buf), 'LCPTXT',   ' = ', this%LCPTXT,           new_line('')
    write(buf,'(A,A)')              trim(buf),                                           new_line('')
    write(buf,'(A,A30,A3,L10,A)')   trim(buf), 'LSPACE',   ' = ', this%LSPACE,           new_line('')
    if (this%LSPACE) then
      write(buf,'(A,A30,A3,I10,A)') trim(buf), 'NACBASIS', ' = ', this%NACBASIS,         new_line('')
      write(buf,'(A,A30,A3,I10,A)') trim(buf), 'NACELE',   ' = ', this%NACELE,           new_line('')
    end if
    write(buf,'(A,A60)') trim(buf), "------------------------------------------------------------"
    write(*,'(A)') trim(buf)
  end subroutine

end module
