! 非绝热耦合模块
! 本模块用于处理分子动力学模拟中的非绝热耦合计算
! 主要功能包括：
! 1. 读取和处理能量本征值和非绝热耦合矩阵
! 2. 支持长时程模拟和短时程模拟
! 3. 处理单粒子和多粒子态的非绝热耦合
! 4. 提供文件IO接口用于数据存储和读取

module couplings
  use prec      ! 提供精度控制和常量定义
  use constants ! 提供常量定义，如MAXSIZE等
  use fileio    ! 提供文件操作接口
  implicit none

  private :: readNaEig, initSpace  ! 将这两个子程序设为私有，只能在模块内部使用

  ! 定义重叠矩阵类型，用于存储非绝热耦合和能量信息
  type overlap
    logical :: LLONG        ! 是否进行长时程模拟
    integer :: LBD          ! 下边界
    integer :: UBD          ! 上边界
    integer :: NBANDS       !! NBASIS - 能带数量
    integer :: TSTEPS       !! NSW - 时间步数
    real(kind=q) :: dt      ! 时间步长
    real(kind=q), allocatable, dimension(:,:,:) :: Dij !! 反厄米矩阵，用于存储非绝热耦合
    ! complex(kind=q), allocatable, dimension(:,:,:) :: DijR  ! 复数形式的实部（已注释）
    ! complex(kind=q), allocatable, dimension(:,:,:) :: DijI  ! 复数形式的虚部（已注释）
    real(kind=q), allocatable, dimension(:,:) :: Eig  ! 能量本征值
  end type

contains

  ! 初始化非绝热耦合矩阵
  ! 输入：
  ! - olap: 输出，多粒子重叠矩阵
  ! - olap_sp: 输出，单粒子重叠矩阵
  ! 流程：
  ! 1. 根据模拟类型设置存储模式
  ! 2. 初始化基本参数
  ! 3. 读取耦合信息
  subroutine TDCoupIJ(olap, olap_sp)
    implicit none
    type(overlap), intent(out) :: olap
    type(overlap), intent(out) :: olap_sp !! 单粒子重叠矩阵

    ! 根据模拟类型设置存储模式
    ! AIMD通常NSW < MAXSIZE，这里存储所有耦合系数
    ! DPMD通常NSW > MAXSIZE，这里存储M+2个系数
    olap%LLONG    = (inp%NSW > MAXSIZE)
    olap_sp%LLONG = (inp%NSW > MAXSIZE)
    if (inp%ALGO == 'FSSH' .AND. olap%LLONG) then
      write(*,*) "[E] The NSW is too long for a FSSH calculation."
      stop
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 初始化部分：设置基本参数
    olap%NBANDS = inp%NBASIS
    olap%TSTEPS = inp%NSW
    olap%dt = inp%POTIM

    olap_sp%NBANDS = inp%BMAX - inp%BMIN + 1
    olap_sp%TSTEPS = inp%NSW
    olap_sp%dt = inp%POTIM
    ! Trotter因子化积分器与复数非绝热耦合不兼容
    if (.NOT. olap%LLONG) then
      ! 短时程模拟：一次性分配所有内存
      allocate(olap%Dij(olap%NBANDS, olap%NBANDS, olap%TSTEPS))
      allocate(olap%Eig(olap%NBANDS, olap%TSTEPS))
      allocate(olap_sp%Dij(olap_sp%NBANDS, olap_sp%NBANDS, olap_sp%TSTEPS))
      allocate(olap_sp%Eig(olap_sp%NBANDS, olap_sp%TSTEPS))
      ! 从EIGTXT和NATXT读取耦合信息
      if (inp%LCPTXT) then
        call readNaEig(olap_sp)
        if (inp%LSPACE) then
          call initSpace(olap, olap_sp)
        else
          olap = olap_sp
        end if
      else
        write(*,*) "[E] This version does not support coupling from COUPCAR, please use NATXT"
        stop
      end if
    else !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! 长时程模拟：分批处理数据
      allocate(olap%Dij(olap%NBANDS, olap%NBANDS, MAXSIZE+2))
      allocate(olap%Eig(olap%NBANDS, MAXSIZE+2))
      allocate(olap_sp%Dij(olap_sp%NBANDS, olap_sp%NBANDS, MAXSIZE+2))
      allocate(olap_sp%Eig(olap_sp%NBANDS, MAXSIZE+2))
      ! 从EIGTXT和NATXT读取耦合信息
      if (inp%LCPTXT) then
        call readLongNaEig(olap, olap_sp)
      else
        write(*,*) "[E] This version does not support coupling from COUPCAR, please use NATXT"
        stop
      end if
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! deallocate(olap_sp%Dij, olap_sp%Eig)
  end subroutine

  ! 读取非绝热耦合和能量本征值
  ! 输入：
  ! - olap: 输入输出，重叠矩阵
  ! 流程：
  ! 1. 打开输入文件
  ! 2. 读取能量本征值和非绝热耦合矩阵
  ! 3. 设置最高占据轨道和最低未占据轨道
  subroutine readNaEig(olap)
    implicit none

    type(overlap), intent(inout) :: olap
    integer :: i, j, k, N, ierr

    ! 打开输入文件
    open(unit=22, file='EIGTXT', status='unknown', action='read', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "[E] IOError: EIGTXT does NOT exist!"
      stop
    end if
    open(unit=23, file='NATXT', status='unknown', action='read', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "[E] IOError: NATXT does NOT exist!"
      stop
    end if

    ! 读取能量本征值和非绝热耦合矩阵
    N = inp%NSW 
    do j=1, N 
      read(unit=22, fmt=*) (olap%Eig(i,j), i=1, olap%NBANDS)
    end do
    do k=1, N
      read(unit=23, fmt=*) ((olap%Dij(j, i, k), j=1, olap%NBANDS), &
                                                   i=1, olap%NBANDS)
    end do

    ! 设置最高占据轨道和最低未占据轨道
    if (inp%LHOLE) then
      call inp%setHLORB(1, olap%NBANDS)
    else
      call inp%setHLORB(olap%NBANDS, 1)
    end if

    ! 将非绝热耦合除以时间步长的两倍
    olap%Dij = olap%Dij / (2*inp%POTIM)

    close(unit=22)
    close(unit=23)
  end subroutine

  ! 读取长时程模拟的非绝热耦合和能量本征值
  ! 输入：
  ! - olap: 输入输出，多粒子重叠矩阵
  ! - olap_sp: 输入输出，单粒子重叠矩阵
  ! 流程：
  ! 1. 分批读取数据
  ! 2. 处理单粒子态
  ! 3. 计算多粒子态
  ! 4. 保存到临时文件
  subroutine readLongNaEig(olap, olap_sp)
    implicit none

    type(overlap), intent(inout) :: olap, olap_sp
    integer :: i, ii, j, k, N, ierr

    !!!!!!!!!!!!!!!!!!!
    ! 存储结构说明：
    ! 数据分三段存储，每段大小为M+2
    ! SAVE:    1 ~  M+2      ! 第一段数据
    !        M+3 ~ 2M+4      ! 第二段数据
    !       2M+5 ~ 3M+6      ! 第三段数据
    ! 这种分段存储方式用于处理大规模数据，避免内存溢出

    ! 打开输入输出文件
    open(unit=22, file='EIGTXT', status='unknown', action='read', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "[E] IOError: EIGTXT does NOT exist!"
      stop
    end if
    open(unit=23, file='NATXT', status='unknown', action='read', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "[E] IOError: NATXT does NOT exist!"
      stop
    end if
    ! 打开临时文件用于存储处理后的数据
    open(unit=42, file='EIGTMP', form='unformatted', status='unknown', &
         access='stream', action='write', iostat=ierr)
    open(unit=43, file='NATMP', form='unformatted', status='unknown',  &
         access='stream', action='write', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "[E] IOError: Can not write to EIGTMP file."
      stop
    end if

    ! 初始化
    N = olap_sp%NBANDS

    ! 设置最高占据轨道和最低未占据轨道
    if (inp%LHOLE) then
      call inp%setHLORB(1, N)
    else
      call inp%setHLORB(N, 1)
    end if

    ! 分批读取和写入数据
    do ii=1, inp%NSW, MAXSIZE+2
      read(unit=22, fmt=*, iostat=ierr) ((olap_sp%Eig(i,k), i=1, N), k=1, MAXSIZE+2)
      read(unit=23, fmt=*, iostat=ierr) (((olap_sp%Dij(j, i, k), j=1, N), i=1, N), k=1, MAXSIZE+2)

      olap_sp%Dij = olap_sp%Dij / (2*inp%POTIM)

      !! 初始化空间
      if (inp%LSPACE) then
        call initSpace(olap, olap_sp)
      else
        olap = olap_sp
      end if

      ! 写入临时文件
      write(unit=42) olap%Eig
      write(unit=43) olap%Dij
    end do

    close(22)
    close(23)
    close(42)
    close(43)
  end subroutine

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 计算多电子态的能量和非绝热耦合
  ! H_ii = SUM_{ik \in ES} Eig_ik - SUM_{ik \in GS} Eig_ik
  ! D_ij = SUM_k d_{ik.jk} PROD_{k'!=k} delta_{ik'.jk'}
  ! i,j : 多电子态标记
  ! k,k': 电子索引
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 初始化多电子态空间
  ! 输入：
  ! - olap: 输入输出，多粒子重叠矩阵
  ! - olap_sp: 输入，单粒子重叠矩阵
  ! 流程：
  ! 1. 计算多粒子态能量
  ! 2. 计算非绝热耦合矩阵
  ! 3. 输出结果到文件
  subroutine initSpace(olap, olap_sp)
    implicit none

    type(overlap), intent(inout) :: olap
    type(overlap), intent(in)    :: olap_sp

    integer :: i, j, k, l, N, band, ierr
    integer :: horb, lorb
    real(q) :: heig, leig
    integer, dimension(inp%NACELE) :: bi, bj
    logical :: beq = .false.

    ! 打开输出文件用于存储计算结果
    open(unit=35, file='ACEIGTXT', status='unknown', action='write', iostat=ierr)
    open(unit=36, file='ACNATXT', status='unknown', action='write', iostat=ierr)

    ! 初始化数组
    olap%Eig = 0.0_q
    olap%Dij = 0.0_q

    ! 计算激发态能量：将单粒子能量叠加得到多粒子态能量
    do i=1,inp%NBASIS
      do j=1,inp%NACELE
        band = inp%BASIS(j,i) - inp%BMIN + 1
        olap%Eig(i,:) = olap%Eig(i,:) + olap_sp%Eig(band,:)
      end do
    end do

    ! 找出最高和最低能量轨道，用于后续能量参考点的选择
    heig = -huge(0.0_q)
    leig = huge(0.0_q)
    do i=1,inp%NBASIS
      if (olap%Eig(i,1) > heig) then
        horb = i
        heig = olap%Eig(i,1)
      end if
      if (olap%Eig(i,1) < leig) then
        lorb = i
        leig = olap%Eig(i,1)
      end if
    end do
    if (inp%LHOLE) then
      call inp%setHLORB(lorb, horb)
    else
      call inp%setHLORB(horb, lorb)
    end if

    ! 减去基态能量，得到相对能量
    do N=1,inp%NSW
      olap%Eig(:,N) = olap%Eig(:,N) - olap%Eig(inp%LORB,N)
    end do

    ! 计算非绝热耦合矩阵元（反厄米）
    ! 使用Slater-Condon规则计算多电子态之间的耦合
    ! 公式说明：
    !         N                  N
    ! H_ij = SUM  h_{i_k' j_k'} PROD delta_{i_k'' j_k''}    
    !        k'=1               k''=1,
    !                           k''!=k'
    ! 其中：
    ! - N = NACELE: 参与耦合的电子数
    ! - i,j: 多电子态标记
    ! - k',k'': 电子索引
    ! - h_{i_k' j_k'}: 单电子耦合矩阵元
    ! - delta_{i_k'' j_k''}: Kronecker delta函数，确保其他电子在相同轨道
    ! - PROD: 表示对所有k''≠k'的电子取乘积
    ! 
    ! 矩阵存储说明：
    ! - Dij存储为Dij(j,i,N)，其中j,i为态索引，N为时间步
    ! - 利用反厄米性质：Dij(i,j) = -Dij(j,i)
    ! - 只计算上三角矩阵，下三角通过反厄米性质获得
    do i=1,inp%NBASIS
      do j=i+1,inp%NBASIS
        ! 将基态索引转换为单粒子态索引
        bi = inp%BASIS(:,i) - inp%BMIN + 1
        bj = inp%BASIS(:,j) - inp%BMIN + 1
        do k=1,inp%NACELE   !! k'
          beq = .true.
          ! 检查除k外的其他电子是否在相同轨道
          ! 这对应于公式中的PROD delta_{i_k'' j_k''}部分
          do l=1,inp%NACELE !! k'' != k'
            if (l == k) cycle
            if (bi(l) /= bj(l)) then
              beq = .false.
              exit
            end if
          end do
          ! 如果其他电子都在相同轨道，则计算k电子的耦合
          ! 这对应于公式中的SUM h_{i_k' j_k'}部分
          if (beq) then
            olap%Dij(j,i,:) = olap%Dij(j,i,:) + olap_sp%Dij(bj(k),bi(k),:)
          end if
        end do 
        ! 利用反厄米性质设置下三角矩阵
        ! 这确保了耦合矩阵满足反厄米性质：Dij(i,j) = -Dij(j,i)
        olap%Dij(i,j,:) = -olap%Dij(j,i,:)
      end do
    end do

    ! 调试模式下不输出文件
    if (inp%DEBUGLEVEL > 0) return

    ! 输出结果到文件
    do N=1,inp%NSW
      write(unit=35,fmt=*) (olap%Eig(i,N), i=1, olap%NBANDS)
      write(unit=36,fmt=*) ((olap%Dij(j,i,N), j=1, olap%NBANDS), &
                                              i=1, olap%NBANDS)
    end do
  end subroutine

  ! 设置重叠矩阵的边界
  ! 输入：
  ! - olap: 输入输出，重叠矩阵
  ! - start: 输入，起始时间步
  ! 流程：
  ! 1. 设置上下边界
  ! 2. 重新分配内存
  ! 3. 从临时文件读取数据
  subroutine setOlapBoundry(olap, start)
    implicit none
    type(overlap), intent(inout) :: olap
    integer, intent(in) :: start
    integer :: ierr

    if (.NOT. olap%LLONG) return

    ! 设置上下边界
    olap%LBD = start - 1
    olap%UBD = olap%LBD + MAXSIZE + 1

    ! 重新分配内存
    if (allocated(olap%Eig)) then
      deallocate(olap%Eig)
      deallocate(olap%Dij)
    end if
    allocate(olap%Eig(olap%NBANDS, olap%LBD : olap%UBD))
    allocate(olap%Dij(olap%NBANDS, olap%NBANDS, olap%LBD : olap%UBD))

    ! 从临时文件读取数据
    open(unit=42, file='EIGTMP', form='unformatted', status='unknown', &
         access='stream', action='read', iostat=ierr)
    open(unit=43, file='NATMP', form='unformatted', status='unknown',  &
         access='stream', action='read', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "[E] Read EIG and NA from TMP file fails."
      stop
    end if

    ! 从起始位置1读取数据
    read(unit=42, iostat=ierr, pos=1+kind(0.0_q)*olap%NBANDS*(olap%LBD-1)) olap%Eig
    read(unit=43, iostat=ierr, pos=1+kind(0.0_q)*olap%NBANDS*olap%NBANDS*(olap%LBD-1)) olap%Dij

    close(42)
    close(43)
  end subroutine

end module
