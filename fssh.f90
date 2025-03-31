! 最少开关跳跃(FSSH)模块
! 本模块实现了Tully的FSSH方法，用于模拟非绝热分子动力学
! 主要功能包括：
! 1. 计算电子态之间的跃迁概率
! 2. 实现表面跳跃动力学
! 3. 处理波函数的时间演化
! 4. 输出波函数和布居数信息
! 5. 支持多能带FSSH计算

module fssh
  use prec      !! 精度定义模块
  use constants !! 常量定义模块
  use fileio    !! 文件IO模块
  use utils     !! 工具函数模块
  use couplings !! 耦合计算模块
  use hamil     !! 哈密顿量模块
  use TimeProp  !! 时间演化模块
  implicit none

  private :: whichToHop, calcProb, calcProbwithDBC

  contains

  ! whichToHop
  ! 输入：
  ! - indion: 当前离子步
  ! - ks: TDKS类型变量，包含波函数信息
  ! - cstat: 当前电子态
  ! 输出：
  ! - which: 选择跃迁的目标态
  ! 流程：
  ! 1. 生成随机数
  ! 2. 根据跃迁概率选择目标态
  subroutine whichToHop(indion, ks, cstat, which)
    implicit none

    integer, intent(in) :: indion, cstat
    integer, intent(out) :: which
    type(TDKS), intent(in) :: ks

    integer :: i
    real(kind=q) :: lower, upper, r

    which = 0
    call random_number(r)  !! 生成[0,1]之间的随机数

    do i=1, ks%ndim
      if (i == 1) then
        lower = 0.0_q
        ! upper = ks%sh_prop(i, indion+1)  !! 旧版本的概率计算方式
        upper = ks%sh_prob(i, cstat, indion)  !! 新版本的概率计算方式
      else
        lower = upper
        ! upper = upper + ks%sh_prop(i, indion+1)  !! 旧版本的累积概率
        upper = upper + ks%sh_prob(i, cstat, indion)  !! 新版本的累积概率
      end if
      if (lower <= r .AND. r < upper) then
        which = i
        exit
      end if
    end do

  end subroutine

  ! calcProb
  ! 输入：
  ! - tion: 当前时间步
  ! - indion: 当前离子步
  ! - tele: 当前电子步
  ! - ks: TDKS类型变量
  ! - olap: 重叠矩阵
  ! 流程：
  ! 1. 计算跃迁概率矩阵
  ! 2. 考虑非绝热耦合效应
  subroutine calcProb(tion, indion, tele, ks, olap)
    implicit none

    type(TDKS), intent(inout) :: ks
    type(overlap), intent(in) :: olap
    integer, intent(in) :: tion, indion, tele

    integer :: i
    real(kind=q) :: Akk_p, Akk_c      !! |c_k|^2，波函数系数的模方
    real(kind=q), dimension(ks%ndim) :: Bkm_p, Bkm_c  !! 非绝热耦合项
    real(kind=q) :: edt               !! 时间步长

    edt = inp%POTIM / inp%NELM

    ! Akk = 0 => P_ij = NaN  !! 当波函数系数为0时，概率计算会出现NaN
    do i = 1, ks%ndim
      Akk_p = REAL(CONJG(ks%psi_p(i)) * ks%psi_p(i), kind=q)  !! 前一步波函数系数的模方
      Akk_c = REAL(CONJG(ks%psi_c(i)) * ks%psi_c(i), kind=q)  !! 当前步波函数系数的模方

      !Changed the matrix structure for NAC. 
      !Now NAC(i,j) is stored as  olap%Dij(j,i,time)
      !! 非绝热耦合矩阵结构改变：NAC(i,j)现在存储在olap%Dij(j,i,time)中
      Bkm_c = REAL(CONJG(ks%psi_c(i)) * ks%psi_c(:) * &
                   calcDij(olap%Dij(:,i,tion-1),olap%Dij(:,i,tion),olap%Dij(:,i,tion+1),tele), kind=q)
      Bkm_p = REAL(CONJG(ks%psi_p(i)) * ks%psi_p(:) * &
                   calcDij(olap%Dij(:,i,tion-1),olap%Dij(:,i,tion),olap%Dij(:,i,tion+1),tele-1), kind=q)

      !! 计算最少开关跳跃方法的跃迁概率
      ks%sh_prob(:, i, indion) = ks%sh_prob(:, i, indion) + &
                                 (Bkm_p / Akk_p + Bkm_c / Akk_c) * edt
    end do

  end subroutine
    
  ! calcProbwithDBC
  ! 输入：
  ! - tion: 当前时间步
  ! - indion: 当前离子步
  ! - ks: TDKS类型变量
  ! - olap: 重叠矩阵
  ! 流程：
  ! 1. 应用细致平衡条件
  ! 2. 考虑能量差对跃迁概率的影响
  subroutine calcProbwithDBC(tion, indion, ks, olap)
    implicit none

    type(TDKS), intent(inout) :: ks
    type(overlap), intent(in) :: olap
    integer, intent(in) :: tion, indion

    integer :: i, j
    real(kind=q) :: dE(ks%ndim), kbT

    kbT = inp%TEMP * BOLKEV

    ! Detail Balance Condition
    ! Beware there is a tricky part: 
    ! P_ij * EXP(-dE/kT) != SUM(P_ij/N * EXP(-dE/N/kT)) = P_ij * EXP(-dE/NkT)
    !! 细致平衡条件
    !! 注意：这里有一个微妙的问题：
    !! P_ij * EXP(-dE/kT) != SUM(P_ij/N * EXP(-dE/N/kT)) = P_ij * EXP(-dE/NkT)
    do i = 1, ks%ndim        
      !! 计算相邻时间步的平均能量差
      dE = ((olap%Eig(:,tion) + olap%Eig(:,tion+1)) - &
            (olap%Eig(i,tion) + olap%Eig(i,tion+1))) / 2.0_q
      if (inp%LHOLE) then
        dE = -dE  !! 空穴情况下的能量差取反
      end if
      !! 应用玻尔兹曼因子修正跃迁概率
      where ( dE > 0.0_q )
        ks%sh_prob(:,i,indion) = ks%sh_prob(:,i,indion) * exp(-dE / kbT)
      end where

      !! 确保概率非负
      forall (j=1:ks%ndim, ks%sh_prob(j,i,indion) < 0) ks%sh_prob(j,i,indion) = 0.0_q
    end do
  end subroutine


  ! runSE
  ! 输入：
  ! - ks: TDKS类型变量
  ! - olap: 重叠矩阵
  ! 流程：
  ! 1. 初始化电子态
  ! 2. 时间演化
  ! 3. 计算跃迁概率
  ! 4. 应用细致平衡条件
  ! 说明：
  ! - 使用不同的时间演化算法
  ! - 包括Trotter、Euler和EulerMod方法
  subroutine runSE(ks, olap)
    implicit none
    type(TDKS), intent(inout) :: ks
    type(overlap), intent(in) :: olap

    integer :: tion, indion, tele
    integer :: istat, cstat
    real(kind=q) :: norm, edt
    logical :: init

    ! init
    istat = inp%INIBAND
    cstat = istat
    edt = inp%POTIM / inp%NELM
    init = .TRUE.

    ! Enter the Loop
    do indion=1, inp%NAMDTIME - 1

      tion = indion + inp%NAMDTINI - 1

      do tele = 1, inp%NELM
        select case (inp%ALGO_INT)
        case (0)
          !! Trotter分解方法
          ks%psi_p = ks%psi_c
          call make_hamil(tion, tele, olap, ks)
          call Trotter(ks, edt)
        case (11) 
          !! 欧拉方法（注意：此方法会累积误差，可能需要重新归一化波函数）
          ks%psi_p = ks%psi_c
          call make_hamil(tion, tele, olap, ks)
          call Euler(ks, edt)
        case (10)
          !! 改进的欧拉方法
          if (init) then
            ks%psi_p = ks%psi_c
            call make_hamil(tion, tele, olap, ks)
            call Euler(ks, edt)
            init = .FALSE.
          else
            call make_hamil2(tion, tele-1, olap, ks)
            call EulerMod(ks, edt)
          end if
        case (2)
          !! 对角化方法
          ks%psi_p = ks%psi_c
          call make_hamil(tion, tele, olap, ks)
          call DiagonizeMat(ks, edt)
          ks%psi_c = HamPsi(ks%ham_c, ks%psi_c, 'N')
        end select
        if (inp%LSHP) call calcProb(tion, indion, tele, ks, olap)
      end do
      !! 检查波函数归一化
      norm = REAL(SUM(CONJG(ks%psi_c) * ks%psi_c), kind=q) 
      if ( norm <= 0.99_q .OR. norm >= 1.01_q)  then
        write(*,*) "[E] Error in Electronic Propagation"
        stop
      end if

      if (inp%LSHP) call calcProbwithDBC(tion, indion, ks, olap)
      !! 保存波函数和布居数
      ks%pop_a(:, indion+1) = REAL(CONJG(ks%psi_c) * ks%psi_c, kind=q)
      ks%psi_a(:, indion+1) = ks%psi_c 

    end do
  end subroutine

  ! runSH
  ! 输入：
  ! - ks: TDKS类型变量
  ! - olap: 重叠矩阵
  ! 流程：
  ! 1. 初始化电子态
  ! 2. 根据跃迁概率进行表面跳跃
  ! 3. 统计各态布居数
  subroutine runSH(ks, olap)
    implicit none

    type(TDKS), intent(inout) :: ks
    type(overlap), intent(in) :: olap
    integer :: i, tion, indion
    integer :: istat, cstat, which

    istat = inp%INIBAND

    do i=1, inp%NTRAJ
      ! 第一步时，当前态等于初始态
      cstat = istat
      do indion=1, inp%NAMDTIME - 1
        tion = indion + inp%NAMDTINI - 1

        call whichToHop(indion, ks, cstat, which)

        if (which > 0) then
          cstat = which  !! 发生跃迁，更新当前态
        end if
        !! 统计各态的布居数
        ks%sh_pops(cstat, indion+1) = ks%sh_pops(cstat, indion+1) + 1
      end do
    end do

    !! 归一化布居数
    ks%sh_pops = ks%sh_pops / inp%NTRAJ
    ks%sh_pops(istat, 1) = 1.0_q

  end subroutine

  ! printSE
  ! 输入：
  ! - ks: TDKS类型变量
  ! - olap: 重叠矩阵
  ! 流程：
  ! 1. 输出波函数和布居数信息
  subroutine printSE(ks, olap)
    implicit none
    type(TDKS), intent(in) :: ks
    type(overlap), intent(in) :: olap

    integer :: i, ierr
    integer :: tion, indion
    character(len=48) :: buf
    character(len=48) :: out_fmt, out_fmt_cmplx

    write(buf, *) inp%NAMDTINI
    write (out_fmt, '( "(f13.2,f11.6, ", I5, "(f11.6))" )' )  ks%ndim
    write (out_fmt_cmplx, '( "(f13.2,f11.6, ", I5, "(SS,f11.6,SP,f9.6", A3, "))" )' )  ks%ndim, '"i"'
    
    !! 打开输出文件
    open(unit=25, file='PSICT.' // trim(adjustl(buf)), status='unknown', action='write', iostat=ierr)
    open(unit=26, file='POPRT.' // trim(adjustl(buf)), status='unknown', action='write', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "[E] IOError: PSICT file I/O error!"
      stop
    end if

    !! 输出波函数和布居数信息
    do indion=1, inp%NAMDTIME
      tion = indion + inp%NAMDTINI - 1
      write(unit=25, fmt=out_fmt_cmplx) indion * inp%POTIM, &
                                  SUM(olap%Eig(:,tion) * ks%pop_a(:,indion)), &
                                  (ks%psi_a(i,indion), i=1, ks%ndim)
      write(unit=26, fmt=out_fmt) indion * inp%POTIM, &
                                  SUM(olap%Eig(:,tion) * ks%pop_a(:,indion)), &
                                  (ks%pop_a(i,indion), i=1, ks%ndim)
    end do

    close(25)
    close(26)

  end subroutine

  ! printSH
  ! 输入：
  ! - ks: TDKS类型变量
  ! - olap: 重叠矩阵
  ! 流程：
  ! 1. 输出表面跳跃结果
  ! 说明：
  ! - 输出到SHPROP文件
  ! - 包含时间、能量和表面跳跃布居数
  subroutine printSH(ks, olap)
    implicit none
    type(TDKS), intent(in) :: ks
    type(overlap), intent(in) :: olap

    integer :: i, ierr
    integer :: tion, indion
    character(len=48) :: buf
    character(len=48) :: out_fmt

    write(buf, *) inp%NAMDTINI
    write (out_fmt, '( "(f13.2,f11.6, ", I5, "(f11.6))" )' )  ks%ndim
    
    !! 打开输出文件
    open(unit=24, file='SHPROP.' // trim(adjustl(buf)), status='unknown', action='write', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "[E] IOError: SHPROP file I/O error!"
      stop
    end if

    !! 输出表面跳跃结果
    do indion=1, inp%NAMDTIME
      tion = indion + inp%NAMDTINI - 1
      write(unit=24, fmt=out_fmt) indion * inp%POTIM, &
                                  SUM(olap%Eig(:,tion) * ks%sh_pops(:,indion)), &
                                  (ks%sh_pops(i,indion), i=1, ks%ndim)
    end do

    close(24)

  end subroutine

  ! printMPFSSH
  ! 输入：
  ! - ks: TDKS类型变量
  ! 流程：
  ! 1. 计算多能带布居数
  ! 2. 输出多能带表面跳跃结果
  ! 说明：
  ! - 输出到MPPOPRT和MPSHPROP文件
  ! - 考虑多能带效应
  subroutine printMPFSSH(ks)
    implicit none
    type(TDKS), intent(inout) :: ks

    integer :: i, ierr
    integer :: bi, tion, indion
    integer, dimension(inp%NACELE) :: bands

    character(len=48) :: buf
    character(len=48) :: out_fmt

    write(buf, *) inp%NAMDTINI
    write (out_fmt, '( "(f13.2,f11.6, ", I5, "(f11.6))" )' )  inp%NBADNS
    
    !! 打开输出文件
    open(unit=51, file='MPPOPRT.' // trim(adjustl(buf)), status='unknown', action='write', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "[E] IOError: MPPOPRT file I/O error!"
      stop
    end if

    !! 计算多能带布居数
    ks%mppop_a = 0.0_q
    do i=1, ks%ndim
      bands = inp%BASIS(:,i) - inp%BMIN + 1
      do bi=1, inp%NACELE
        ks%mppop_a(bands(bi),:) = ks%mppop_a(bands(bi),:) + &
                                  ks%pop_a(i,:)
      end do
    end do

    !! 输出多能带布居数
    do indion=1, inp%NAMDTIME
      tion = indion + inp%NAMDTINI - 1
      write(unit=51, fmt=out_fmt) indion * inp%POTIM, & 
                                  ! SUM(ks%eigKs(:,tion) * ks%pop_a(:,indion)) / inp%NACELE, &
                                  (ks%mppop_a(i,indion), i=1, inp%NBADNS)
    end do

    close(51)

    !! 表面跳跃结果
    if (inp%LSHP) then
      !! 打开输出文件
      open(unit=52, file='MPSHPROP.' // trim(adjustl(buf)), status='unknown', action='write', iostat=ierr)
      if (ierr /= 0) then
        write(*,*) "[E] IOError: MPSHPROP file I/O error!"
        stop
      end if
  
      !! 计算多能带表面跳跃布居数
      ks%sh_mppops = 0.0_q
      do i=1, ks%ndim
        bands = inp%BASIS(:,i) - inp%BMIN + 1
        do bi=1, inp%NACELE
          ks%sh_mppops(bands(bi),:) = ks%sh_mppops(bands(bi),:) + &
                                      ks%sh_pops(i,:)
        end do
      end do
  
      !! 输出多能带表面跳跃结果
      do indion=1, inp%NAMDTIME
        tion = indion + inp%NAMDTINI - 1
        write(unit=52, fmt=out_fmt) indion * inp%POTIM, & 
                                    ! SUM(ks%eigKs(:,tion) * ks%sh_pops(:,indion)) / inp%NACELE, &
                                    (ks%sh_mppops(i,indion), i=1, inp%NBADNS)
      end do
  
      close(52)
    end if

  end subroutine

end module
