! 非绝热分子动力学主程序
! 本程序实现了基于FSSH和DISH的非绝热分子动力学模拟
! 主要功能包括：
! 1. 初始化MPI并行环境
! 2. 读取和处理输入参数
! 3. 计算非绝热耦合
! 4. 执行时间演化
! 5. 实现表面跳跃
! 6. 输出模拟结果
! 7. 支持多能带计算

Program main
  use prec      !! 精度定义模块
  use fileio    !! 文件IO模块
  use utils     !! 工具函数模块
  use parallel  !! 并行计算模块
  use couplings !! 耦合计算模块
  use hamil     !! 哈密顿量模块
  use fssh      !! 最少跳跃表面跳跃模块
  use dish      !! 退相干诱导表面跳跃模块
#ifdef ENABLEMPI
  use mpi       !! MPI并行计算模块
#endif

  implicit none

  type(TDKS) :: ks           !! TDKS类型变量
  type(overlap) :: olap, olap_sp  !! 重叠矩阵

  integer :: ns, cr, cm, t1, t2, ttot1, ttot2
  integer :: nprog = 1, iprog = 0, ierr
  integer :: communicator = 0, color = 0
  integer :: lower, upper, ncount
  real(kind=qs) :: tottime
  character(len=256) :: buf

#ifdef ENABLEMPI
  !! 初始化MPI
  call MPI_INIT(ierr)
  if (ierr /= 0) then
    write(*,*) "[E] MPI initialization failed. Aborting..."
    call MPI_FINALIZE(ierr)
    stop
  end if
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprog, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, iprog, ierr)
#endif

  !! 初始化系统时钟
  call system_clock(count_rate=cr)
  call system_clock(count_max=cm)
  
  if (iprog == 0) call printWelcome()
  !! 初始化随机数种子
  call init_random_seed(salt=iprog)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! 第一步：获取用户输入
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call inp%getInstance()
  call inp%getUserInp(nprog)

#ifdef ENABLEMPI
  !! 设置MPI通信器
  color = MODULO(iprog, inp%NPAR)
  call MPI_COMM_SPLIT(MPI_COMM_WORLD, color, iprog, communicator, ierr)
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  call MPI_COMM_SIZE(communicator, nprog, ierr)
  call MPI_COMM_RANK(communicator, iprog, ierr)
#endif
  call inp%setMPI(iprog, nprog, color, communicator)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! 第二步：获取耦合信息
  !! 第一次运行时，以下子程序将从WAVECAR计算非绝热耦合NAC并写入二进制文件COUPCAR中
  !! 后续运行时，直接从该文件读取NAC
  !! 对于一般NAMD运行，该文件过大，解决方案是只将需要的信息写入纯文本文件
  !! 如果存在此类文件（在inp中设置LCPTXT = .TRUE.），则可跳过大型二进制文件
  !! 直接读取纯文本文件，这将在'initTDKS'子程序中完成
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call TDCoupIJ(olap, olap_sp)
  call system_clock(ttot1)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! 分配任务
  call divideTasks(color, inp%NPAR, inp%NSAMPLE, lower, upper, ncount)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ns=lower, upper
    call inp%setIni(ns)
    if (iprog == 0) call inp%printUserInp()

    !! 初始化KS矩阵
    call system_clock(t1)
    call initTDKS(ks)
    call system_clock(t2)
    if (iprog == 0) then
      write(buf,'(A3,I3,A5,I4,A,T48,F11.3)') "MPI", color, " TINI", inp%NAMDTINI, ": CPU Time in initTDKS [s]:", &
                                              MODULO(t2-t1, cm) / REAL(cr)
      write(*,*) trim(buf)
    end if

    select case(inp%ALGO)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case ('FSSH')
      if (iprog == 0) then
        !! 时间演化
        t1 = t2
        call runSE(ks, olap)
        call system_clock(t2)
        write(buf,'(A3,I3,A5,I4,A,T48,F11.3)') "MPI", color, " TINI", inp%NAMDTINI, ": CPU Time in runSE [s]:", &
                                                MODULO(t2-t1, cm) / REAL(cr)
        write(*,*) trim(buf)

        t1 = t2
        call printSE(ks, olap)
        call system_clock(t2)
        write(buf,'(A3,I3,A5,I4,A,T48,F11.3)') "MPI", color, " TINI", inp%NAMDTINI, ": CPU Time in printSE [s]:", &
                                                MODULO(t2-t1, cm) / REAL(cr)
        write(*,*) trim(buf)
        !! 运行表面跳跃
        if (inp%LSHP) then
          t1 = t2
          call runSH(ks, olap)
          call system_clock(t2)
          write(buf,'(A3,I3,A5,I4,A,T48,F11.3)') "MPI", color, " TINI", inp%NAMDTINI, ": CPU Time in runSH [s]:", &
                                                  MODULO(t2-t1, cm) / REAL(cr)
          write(*,*) trim(buf)

          t1 = t2
          call printSH(ks, olap)
          call system_clock(t2)
          write(buf,'(A3,I3,A5,I4,A,T48,F11.3)') "MPI", color, " TINI", inp%NAMDTINI, ": CPU Time in printSH [s]:", &
                                                  MODULO(t2-t1, cm) / REAL(cr)
          write(*,*) trim(buf)
        end if
        if (inp%LSPACE) then
          t1 = t2
          call printMPFSSH(ks)
          call system_clock(t2)
          write(buf,'(A3,I3,A5,I4,A,T48,F11.3)') "MPI", color, " TINI", inp%NAMDTINI, ": CPU Time in printMPFSSH [s]:", &
                                                  MODULO(t2-t1, cm) / REAL(cr)
          write(*,*) trim(buf)
        end if
      end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case ('DISH')
      t1 = t2
      call runDISH(ks, olap)
      call system_clock(t2)
      if (iprog == 0) then
        write(buf,'(A3,I3,A5,I4,A,T48,F11.3)') "MPI", color, " TINI", inp%NAMDTINI, ": CPU Time in runDISH [s]:", &
                                                MODULO(t2-t1, cm) / REAL(cr)
        write(*,*) trim(buf)
      end if
    end select
  end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call system_clock(ttot2)
  tottime = MODULO(ttot2-ttot1, cm) / REAL(cr)
#ifdef ENABLEMPI
  call MPI_REDUCE(tottime+0.0, tottime, 1, MPI_REAL, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
#endif
  
  if (iprog == 0 .AND. color == 0) then
    write(*,'(A)') "------------------------------------------------------------"
    write(buf,'(A,T48,F11.3)') "All Time Elapsed [s]:", tottime
    write(*,*) trim(buf)
  end if

#ifdef ENABLEMPI
  call MPI_COMM_FREE(communicator, ierr)
  call MPI_FINALIZE(ierr)
#endif

contains

  ! printWelcome
  ! 流程：
  ! 1. 打印欢迎信息
  ! 2. 显示程序版本和作者信息
  ! 3. 显示支持的输入参数和文件
  subroutine printWelcome()
    implicit none
    character(len=8)  :: date
    character(len=10) :: time
    character(len=5)  :: zone
    !! 打印ASCII大标题
    write(*,'(A)') "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,'(A)') "| __          ________ _      _____ ____  __          __ ______   _______ ____   |"
    write(*,'(A)') "| \ \        / /  ____| |    / ____/ __ \|  \        /  |  ____| |__   __/ __ \  |"
    write(*,'(A)') "|  \ \  /\  / /| |__  | |   | |   | |  | | \ \      / / | |__       | | | |  | | |"
    write(*,'(A)') "|   \ \/  \/ / |  __| | |   | |   | |  | | |\ \    / /| |  __|      | | | |  | | |"
    write(*,'(A)') "|    \  /\  /  | |____| |___| |___| |__| | | \ \  / / | | |____     | | | |__| | |"
    write(*,'(A)') "|     \/_ \/ _ |______|______\_____\____/| |  \ \/ /  | |______|  _ |_|  \____/  |"
    write(*,'(A)') "|      | |  | |  ____|    | \ | |   /\   | |   \  /   | |  __ \  | |             |"
    write(*,'(A)') "|      | |__| | |__ ______|  \| |  /  \  | |    \/    | | |  | | | |             |"
    write(*,'(A)') "|      |  __  |  __|______| . ` | / /\ \ | |          | | |  | | | |             |"
    write(*,'(A)') "|      | |  | | |         | |\  |/ ____ \| |          | | |__| | |_|             |"
    write(*,'(A)') "|      |_|  |_|_|         |_| \_/_/    \_\_|          |_|_____/  (_)             |"
    write(*,'(A)') "|                                                                                |"
    write(*,'(A)') "| Version: Jan. 2023.                                                            |"
    write(*,'(A)') "| Authors:                                                                       |"
    write(*,'(A)') "| Big Thanks to:                                                                 |"
    write(*,'(A)') "|                                                                                |"
    write(*,'(A)') "| Supported input paramters:                                                     |"
    write(*,'(A)') "|     BMIN, BMAX,                                                                |"
    write(*,'(A)') "|     NSAMPLE, NTRAJ, NSW, NELM,                                                 |"
    write(*,'(A)') "|     TEMP, NAMDTIME, POTIM,                                                     |"
    write(*,'(A)') "|     LHOLE, LSHP, ALGO, ALGO_INT, LCPTXT,                                       |"
    write(*,'(A)') "|     LSPACE, NACBASIS, NACELE,                                                  |"
    write(*,'(A)') "|     NPAR, LBINOUT,                                                             |"
    write(*,'(A)') "|     RUNDIR, TBINIT, DIINIT, SPINIT,                                            |"
    write(*,'(A)') "|     DEBUGLEVEL                                                                 |"
    write(*,'(A)') "| Supported input files:                                                         |"
    write(*,'(A)') "|     inp, COUPCAR, EIGTXT, NATXT, INICON, DEPHTIME, ACSPACE                     |"
    write(*,'(A)') "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,'(A)') "                                                                                  "
    !! 打印程序启动时间
    call date_and_time(date, time, zone)
    write(*,'(A21,X,A4,".",A2,".",A2,X,A2,":",A2,":",A6,X,"UTC",A5,/)') 'The program starts at',          &
                                                                        date(1:4), date(5:6), date(7:8),  &
                                                                        time(1:2), time(3:4), time(5:10), &
                                                                        zone
  end subroutine

end Program
