! 时间演化模块
! 本模块实现了多种时间演化算法，用于求解含时薛定谔方程
! 主要功能包括：
! 1. Trotter分解算法（二阶精度）
! 2. Euler方法（一阶精度）
! 3. 改进的Euler方法
! 4. 矩阵对角化方法
! 5. 支持BLAS加速计算

module TimeProp
  use prec      ! 提供精度控制
  use constants ! 提供物理常量
  use utils     ! 提供工具函数
  use fileio    ! 提供文件操作接口
  use hamil     ! 提供哈密顿量相关操作

  implicit none

contains 
  
  !! 使用Trotter公式求解含时薛定谔方程
  !! 该算法由Akimov, A. V., 和 Prezhdo, O. V.在J. Chem. Theory Comput. 2014, 10, 2, 789–804提出，由Dr.Li yunhai(liyunhai1016@hotmail.com)修改
  !! 特点：
  !! 1. 提供稳健高效的积分方案
  !! 2. 电子时间步长(POTIM/NELM)可增加到核时间步长（并非所有情况）
  !! 3. 使用小NELM前需检查收敛性
  !! 4. 要求哈密顿量的非对角元素为实数（不含虚数单位）
  subroutine Trotter(ks, edt)
    implicit none

    type(TDKS), intent(inout) :: ks
    real(kind=q), intent(in)  :: edt
    integer :: jj, kk
    complex(kind=q) :: phi, cos_phi, sin_phi, cjj, ckk

    !! 根据Liouville-Trotter算法演化波函数
    !! exp[i*(L_{ij}+L_i)*dt/2]
    !! = exp(i*L_{ij}*dt/2) * exp(i*L_i*dt/2) * exp(i*L_i*dt/2) * exp(i*L_{ij}*dt/2)
    !! = exp(i*L_{ij}*dt/2) * exp(i*L_i*dt) * exp(i*L_{ij}*dt/2)

    !! 第一部分：L_{ij}项
    !! 注意：非绝热耦合矩阵结构已更改
    !! 现在Ham_c(i,j)存储为ks%ham_c(j,i)
    !! 作者：Weibin

    do jj = 1, ks%ndim
      do kk = jj+1, ks%ndim
        phi = 0.5_q * edt * con%miI * ks%ham_c(kk, jj) / con%hbar
        cos_phi = cos(phi)
        sin_phi = sin(phi)
        cjj = ks%psi_c(jj)
        ckk = ks%psi_c(kk)
        ks%psi_c(jj) =  cos_phi * cjj + sin_phi * ckk
        ks%psi_c(kk) = -sin_phi * cjj + cos_phi * ckk
      end do
    end do
    
    !! L_i项
    !! 注意：非绝热耦合矩阵结构已更改
    !! 现在Ham_c(i,j)存储为ks%ham_c(j,i)
    !! 作者：Weibin

    do jj = 1, ks%ndim
      phi = edt * con%miI * ks%ham_c(jj, jj) / con%hbar
      ks%psi_c(jj) = ks%psi_c(jj) * exp(phi)
    end do

    !! 第二部分：L_{ij}项
    do jj = ks%ndim, 1, -1
      do kk = ks%ndim, jj+1, -1
        phi = 0.5_q * edt * con%miI * ks%ham_c(kk, jj) / con%hbar
        cos_phi = cos(phi)
        sin_phi = sin(phi)
        cjj = ks%psi_c(jj)
        ckk = ks%psi_c(kk)
        ks%psi_c(jj) =  cos_phi * cjj + sin_phi * ckk
        ks%psi_c(kk) = -sin_phi * cjj + cos_phi * ckk
      end do
    end do
  end subroutine
  
  !! Trotter算法的矩阵形式实现
  subroutine TrotterMat(ks, edt)
    implicit none
    type(TDKS), intent(inout) :: ks
    real(kind=q), intent(in)  :: edt

    integer :: jj, kk
    complex(kind=q) :: cos_phi, sin_phi
    complex(kind=q), dimension(ks%ndim,ks%ndim) :: identity, ham_l, ham_r
    complex(kind=q), dimension(ks%ndim) :: cjj, ckk

    !! 输入Ham_c(j,i)，输出Ham_c(i,j)
    identity = zDIAG(ks%ndim, [(con%uno, jj=1,ks%ndim)])
    !! 计算相位
    ks%ham_c(:,:) = 0.5_q * edt * con%miI * ks%ham_c(:,:) / con%hbar

    ham_l = identity
    !! 第一部分：L_{ij}项
    do jj = 1, ks%ndim
      do kk = jj+1, ks%ndim
        cos_phi = cos(ks%ham_c(kk, jj))
        sin_phi = sin(ks%ham_c(kk, jj))
        cjj(:) = cos_phi * ham_l(:,jj) - sin_phi * ham_l(:,kk)
        ckk(:) = sin_phi * ham_l(:,jj) + cos_phi * ham_l(:,kk)
        ham_l(:,jj) = cjj(:)
        ham_l(:,kk) = ckk(:)
        ! ham_r = identity
        ! ham_r(jj,jj) =  cos_phi
        ! ham_r(kk,kk) =  cos_phi
        ! ham_r(jj,kk) =  sin_phi
        ! ham_r(kk,jj) = -sin_phi
        ! ham_l = matmul(ham_l, ham_r)
      end do
    end do

    ! ham_r = identity
    do jj = 1, ks%ndim
      ham_l(:,jj) = ham_l(:,jj) * exp(2 * ks%ham_c(jj,jj))
      ! ham_r(jj,jj) = exp(2 * ks%ham_c(jj,jj))
    end do
    ! ham_l = matmul(ham_l, ham_r)

    do jj = ks%ndim, 1, -1
      do kk = ks%ndim, jj+1, -1
        cos_phi = cos(ks%ham_c(kk, jj))
        sin_phi = sin(ks%ham_c(kk, jj))
        cjj(:) = cos_phi * ham_l(:,jj) - sin_phi * ham_l(:,kk)
        ckk(:) = sin_phi * ham_l(:,jj) + cos_phi * ham_l(:,kk)
        ham_l(:,jj) = cjj(:)
        ham_l(:,kk) = ckk(:)
        ! ham_r = identity
        ! ham_r(jj,jj) =  cos_phi
        ! ham_r(kk,kk) =  cos_phi
        ! ham_r(jj,kk) =  sin_phi
        ! ham_r(kk,jj) = -sin_phi
        ! ham_l = matmul(ham_l, ham_r)
      end do
    end do

    ks%ham_c = ham_l
  end subroutine

  !! 一阶Euler方法
  subroutine Euler(ks, edt)  
    implicit none
    type(TDKS), intent(inout) :: ks
    real(kind=q), intent(in) :: edt

    ks%psi_c = ks%psi_c &
               & - con%I * edt / con%hbar * matmul(transpose(ks%ham_c), ks%psi_c)
  end subroutine

  !! 改进的Euler方法
  !! 第一步应使用Euler方法
  subroutine EulerMod(ks, edt)
    implicit none
    type(TDKS), intent(inout) :: ks
    real(kind=q), intent(in) :: edt

    !! 检查内存分配
    if ((inp%DEBUGLEVEL <= 0 .AND. (.NOT. allocated(ks%psi_p)) .OR. (.NOT. allocated(ks%psi_n)))) then
      write(*,*) '[E] Did not allocate memory for psi_p and psi_n'
      stop
    end if

    ! if (.NOT. init) then
      ks%psi_n = ks%psi_p &
                 & - 2 * con%I * edt / con%hbar * matmul(transpose(ks%ham_c), ks%psi_c)
    ! else
    !   ks%psi_n = ks%psi_c &
    !              & - con%I * edt / con%hbar * matmul(transpose(ks%ham_c), ks%psi_c)
    ! end if
    ks%psi_p = ks%psi_c
    ks%psi_c = ks%psi_n
  end subroutine

  !! Euler方法的矩阵形式
  subroutine EulerMat(ks, edt)
    implicit none
    type(TDKS), intent(inout) :: ks
    real(kind=q), intent(in) :: edt
  
    ks%ham_c = con%miI * edt / con%hbar * ks%ham_c
  end subroutine

  !! 改进Euler方法的矩阵形式
  subroutine EulerModMat(ks, edt)
    implicit none
    type(TDKS), intent(inout) :: ks
    real(kind=q), intent(in) :: edt
  
    ks%ham_c = 2 * con%miI * edt / con%hbar * ks%ham_c
  end subroutine

  !! 基于对角化的时间演化方法
  !! 使用MKL库进行矩阵对角化
  subroutine DiagonizeMat(ks, edt)
    implicit none
    type(TDKS), intent(inout) :: ks
    real(kind=q), intent(in) :: edt

    integer :: INFO, LWORK
    real(kind=q),    dimension(ks%ndim)         :: eig
    complex(kind=q), dimension(ks%ndim,ks%ndim) :: exp_eig
    complex(kind=q), dimension(2*ks%ndim-1)     :: WORK
    real(kind=q),    dimension(3*ks%ndim-2)     :: RWORK

    complex(kind=q), dimension(ks%ndim,ks%ndim) :: C1, C2
    ! complex(kind=q), dimension(ks%ndim)         :: Y

#ifndef ENABLEMKL
    write(*,*) "[E] Diagnization Algorithm is not implemented without MKL."
    stop
#else
    !! 使用LAPACK的zHEEV进行矩阵对角化
    !! 语法为zHEEV(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,RWORK,INFO)
    !! 参考：https://netlib.org/lapack/explore-html
    LWORK = 2*ks%ndim - 1
    call zHEEV('V', 'U',                  &
              ks%ndim, ks%ham_c, ks%ndim, &
              eig,                        &
              WORK, LWORK, RWORK,         &
              INFO)
    if (INFO == 0) then
      exp_eig = zDIAG(ks%ndim, exp(con%miI*edt/con%hbar*eig))
    else
      write(*,*) '[E] Hamiltonian Diagonization Error. Code: ', INFO
      stop
    end if

    !! zGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    !! C2 = matmul(matmul(conjg(ks%ham_c), exp_eig), transpose(ks%ham_c)) !! H_ij
    !! JOBZ='C' <=> conjg * transpose
    C1 = con%cero
    C2 = con%cero
    call zGEMM('N', 'N',                                                    &
              ks%ndim, ks%ndim, ks%ndim, con%uno, conjg(ks%ham_c), ks%ndim, &
              exp_eig, ks%ndim,                                             &
              con%cero, C1, ks%ndim)
    call zGEMM('N', 'T',                                                    &
              ks%ndim, ks%ndim, ks%ndim, con%uno, C1, ks%ndim,              &
              ks%ham_c, ks%ndim,                                            &
              con%cero, C2, ks%ndim)

    ks%ham_c = C2
    !! zGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    ! Y = con%cero
    ! call zGEMV('N',                                   &
    !           ks%ndim, ks%ndim, con%uno, C2, ks%ndim, &
    !           ks%psi_c, 1,                            &
    !           con%cero, Y, 1)
    ! ks%psi_c = Y
#endif
  end subroutine

#ifndef ENABLEMKL
  !! 不使用BLAS的哈密顿量-波函数乘法
  function HamPsi(ham, psi, jobz)
    implicit none
    complex(kind=q), intent(in), dimension(:,:) :: ham
    complex(kind=q), intent(in), dimension(:)   :: psi
    character(len=1), intent(in) :: jobz
    
    complex(kind=q), dimension(size(psi)) :: HamPsi

    select case (jobz)
    case ('N')
      HamPsi = matmul(ham, psi)
    case ('T')
      HamPsi = matmul(transpose(ham), psi)
    end select
  end function
#else
  !! 使用BLAS的哈密顿量-波函数乘法
  function HamPsi(ham, psi, jobz) result(HamPsiBlas)
    implicit none
    complex(kind=q), intent(in), dimension(:,:) :: ham
    complex(kind=q), intent(in), dimension(:)   :: psi
    character(len=1), intent(in) :: jobz
    
    complex(kind=q), dimension(size(psi)) :: HamPsiBlas, Y
    integer :: N

    !! zGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    N = size(psi)
    Y = con%cero
    call zGEMV(jobz,                 &
              N, N, con%uno, ham, N, &
              psi, 1,                &
              con%cero, Y, 1)
    HamPsiBlas = Y
  end function
#endif

end module
