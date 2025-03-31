! 精度和常量定义模块
! 本模块定义了程序中使用的精度类型和物理常量及相互转换

! VASP中使用的精度定义
module prec
  integer, parameter :: q =SELECTED_real_KIND(10)  !! 高精度实数类型(10位有效数字)
  integer, parameter :: qs=SELECTED_real_KIND(5)   !! 单精度实数类型(5位有效数字)
end module


module constants
  use prec

  !! 复数常量定义
  type, private :: constant
    complex(q) :: I = (0.0_q, 1.0_q)      !! 虚数单位 i
    complex(q) :: miI = (0.0_q, -1.0_q)   !! -i
    complex(q) :: cero = (0.0_q, 0.0_q)   !! 复数零
    complex(q) :: uno = (1.0_q, 0.0_q)    !! 复数一
    real(q)    :: hbar = 0.6582119281559802_q  !! 约化普朗克常数 (eV·fs)
  end type
  
  type(constant), protected :: con

  !! 重要物理常数（转换为原子单位）
  !! - AUTOA  = 1 a.u. 对应的埃米
  !! - RYTOEV = 1 Ry 对应的电子伏特
  !! - EVTOJ  = 1 eV 对应的焦耳
  !! - AMTOKG = 1 原子质量单位（质子质量）对应的千克
  !! - BOLKEV = 玻尔兹曼常数 (eV/K)
  !! - BOLK   = 玻尔兹曼常数 (J/K)
  real(q), parameter :: AUTOA=0.529177249_q, RYTOEV=13.605826_q
  real(q), parameter :: CLIGHT = 137.037  !! 光速（原子单位）
  real(q), parameter :: EVTOJ=1.60217733E-19_q, AMTOKG=1.6605402E-27_q, &
                        BOLKEV=8.6173857E-5_q, BOLK=BOLKEV*EVTOJ
  real(q), parameter :: EVTOKCAL=23.06    !! 1 eV 对应的千卡

  !! 元电荷相关常数
  !! FELECT = (电子电荷)/(4π·真空介电常数)，在原子单位中为 e^2
  !! EDEPS = 电子电荷/真空介电常数，在原子单位中为 4πe^2
  !! HSQDTM = (普朗克常数/(2π))^2/(2·电子质量)
  real(q), parameter  :: PI =3.141592653589793238_q, TPI=2*PI
  real(q), parameter  :: FELECT = 2*AUTOA*RYTOEV, EDEPS=4*PI*2*RYTOEV*AUTOA, &
                         HSQDTM = RYTOEV*AUTOA*AUTOA
  complex(q), parameter  :: CITPI = (0._q,1._q)*TPI  !! 2πi

  !! 磁矩相关常数
  !! 矢量场A与动量相乘再乘以e/(2m_ec)得到能量
  !! 磁矩以玻尔磁子为单位
  !! e/(2m_ec) A(r) p(r) = 能量
  !! e/(2m_ec) m_s × (r-r_s)/(r-r_s)^3 hbar ∇ = 
  !! e^2 hbar^2/(2m_e^2 c^2) 1/长度^3 = 能量
  !! 磁矩到能量的转换因子
  !! 已在SI单位中由Gilles de Wijs独立验证
  real(q), parameter :: MAGMOMTOENERGY=1/CLIGHT**2*AUTOA**3*RYTOEV

  !! 连接输入和输出磁矩的无量纲数
  !! AUTOA e^2/(2m_ec^2)
  real(q), parameter :: MOMTOMOM=AUTOA/CLIGHT/CLIGHT/2
  real(q), parameter :: AUTOA2=AUTOA *AUTOA
  real(q), parameter :: AUTOA3=AUTOA2*AUTOA
  real(q), parameter :: AUTOA4=AUTOA2*AUTOA2
  real(q), parameter :: AUTOA5=AUTOA3*AUTOA2
  
  !! 原子单位到德拜的偶极矩转换因子
  real(q), parameter :: AUTDEBYE=2.541746

  integer, parameter :: MAXSIZE = 1E5  !! 最大数组大小 (2^23)
  
end module constants
