! 并行计算模块
! 本模块用于处理并行计算任务

module parallel
  implicit none
contains

  !! divideTasks
  !! 输入：
  !!   iprog: 当前进程编号
  !!   ngroup: 进程组数量
  !!   ntask: 总任务数
  !! 输出：
  !!   lower: 当前进程负责的起始任务编号
  !!   upper: 当前进程负责的结束任务编号
  !!   nout: 当前进程负责的任务数量
  !! 说明：
  !!   将ntask个任务均匀分配给ngroup个进程
  !!   当任务数不能被进程数整除时，前remainder个进程多分配一个任务
  !!   例如：ntask=8, ngroup=5时
  !!   进程0: 任务1-2 (2个任务)
  !!   进程1: 任务3-4 (2个任务)
  !!   进程2: 任务5-6 (2个任务)
  !!   进程3: 任务7 (1个任务)
  !!   进程4: 任务8 (1个任务)
  subroutine divideTasks(iprog, ngroup, ntask, lower, upper, nout)
    integer, intent(in)    :: iprog, ntask
    integer, intent(inout) :: ngroup
    integer, intent(out)   :: lower, upper, nout

    integer :: quotient, remainder

    !! 确保进程组数不超过任务数
    ngroup = MIN(ngroup, ntask)
    !! 计算每个进程的基础任务数和余数
    quotient  = ntask / ngroup
    remainder = MOD(ntask, ngroup)

    !! 分配任务：
    !! 1. 前remainder个进程各分配quotient+1个任务
    !! 2. 剩余进程各分配quotient个任务
    if (iprog + 1 > remainder) then
      nout = quotient
      lower = remainder + iprog*quotient + 1
      upper = lower + nout - 1
    else
      nout = quotient + 1
      lower = iprog*(quotient+1) + 1
      upper = lower + nout - 1
    end if
  end subroutine

end module parallel