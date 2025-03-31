#!/usr/bin/env python
############################################################
import os, re
import numpy as np
from glob import glob


############################################################
def EnergyFromPro(infile="PROCAR"):
    """
    从PROCAR文件中提取能带能量。
    """
    print(infile)
    assert os.path.isfile(infile), "%s cannot be found!" % infile
    FileContents = [line for line in open(infile) if line.strip()]

    # 当能带数量太大时，";"和实际能带数之间可能没有空格。这是由Homlee Guo发现的一个bug。
    # 这里，#kpts、#bands和#ions都是整数
    nkpts, nbands, nions = [
        int(xx) for xx in re.sub("[^0-9]", " ", FileContents[1]).split()
    ]

    energies = np.asarray(
        [line.split()[-4] for line in FileContents if "occ." in line], dtype=float
    )
    nspin = energies.shape[0] // (nkpts * nbands)
    energies.resize(nspin, nkpts, nbands)

    return energies


def parallel_energy(runDirs, nproc=None):
    """
    并行计算指定目录中的局域化能量。
    """
    import multiprocessing

    nproc = multiprocessing.cpu_count() if nproc is None else nproc
    pool = multiprocessing.Pool(processes=nproc)

    results = []
    for rd in runDirs:
        res = pool.apply_async(EnergyFromPro, (rd + "/PROCAR",))
        results.append(res)

    en = []
    for ii in range(len(results)):
        tmp_en = results[ii].get()
        en.append(tmp_en)

    return np.array(en)


############################################################
# 计算空间局域化
############################################################
nsw = 2000  # 如果all_en.npy存在则无效
nproc = 2  # 并行进程数
prefix = "NAMD/run/"  # 运行目录前缀
runDirs = [prefix + "{:04d}".format(ii + 1) for ii in range(nsw)]

# 针对Gamma点版本，无自旋
if os.path.isfile("all_en.npy"):
    energies = np.load("all_en.npy")
else:
    energies = parallel_energy(runDirs, nproc=nproc)
    print(energies.shape)
    energies = energies[:, 0, 0, :]
    np.save("all_en.npy", energies)

############################################################
# 加载FSSH结果文件
############################################################

########################################
bmin = 310  # 起始能带序号（从1开始）
bmax = 324  # 结束能带序号（从1开始）
namdTime = 1000  # NAMD模拟时间步数
potim = 1.0  # 时间步长
Eref = "EVBM"  # 能量参考点 'ECBM'（导带底）或'EVBM'（价带顶）
inpFiles = glob("MPSHPROP.*")  # 表面跳跃（MP）属性文件
########################################
### po.npy: 3D数组，维度为3*namdtime*(nband+2)，包含初始时间|平均能量|能带
if not os.path.isfile("po.npy"):

    iniTimes = [int(F.split(".")[-1]) - 1 for F in inpFiles]

    dat = np.array([np.loadtxt(F) for F in inpFiles])
    dat = np.average(dat, axis=0)  # 系综平均

    Ci_t = dat[:, 2:]
    En_t = np.zeros_like(Ci_t)
    Time = np.zeros_like(Ci_t)

    energies = energies[:, (bmin - 1) : bmax]
    EREF = {"ECBM": np.average(energies[:, 0]), "EVBM": np.average(energies[:, -1])}

    for start in iniTimes:
        end = start + namdTime
        En_t += energies[start:end, :] - EREF[Eref]
    else:
        En_t /= len(iniTimes)  # 时间平均能量

    for i in range(bmax - bmin + 1):
        Time[:, i] = np.arange(namdTime)

    # 时间、平均能量，用于绘图
    # dat[:,1] -= EREF[Eref] # MP
    x, y = dat[:, :2].T

    np.save("po.npy", (Time, En_t, Ci_t))
    np.savetxt("average_energy.dat", dat[:, :2], fmt="%8.4f")

else:
    inp = np.load("po.npy")
    Time = inp[0, :]
    En_t = inp[1, :]
    Ci_t = inp[2, :]
    x, y = np.loadtxt("average_energy.dat", unpack=True)

############################################################
# 图形绘制设置
############################################################
import matplotlib as mpl

mpl.use("agg")
mpl.rcParams["axes.unicode_minus"] = False  # 解决负号显示问题

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure()
fig.set_size_inches(4.8, 3.0)  # 设置图像尺寸

ax = plt.subplot()
divider = make_axes_locatable(ax)
ax_cbar = divider.append_axes("right", size="5%", pad=0.02)  # 为色条创建额外的轴

############################################################
# 绘制平均能量随时间的变化曲线
############################################################
(line,) = ax.plot(x, y, ls="--", color="blue", lw=1.5, alpha=0.6)

# 使用散点图展示能量分布和占据概率
kmap = ax.scatter(
    Time,
    En_t,
    c=Ci_t / 2,  # /2是针对MP的处理
    cmap="hot_r",  # 色图选择
    vmin=0,
    vmax=1,
    s=15,
    alpha=0.8,
    lw=0.0,
)
cbar = plt.colorbar(
    kmap, cax=ax_cbar, orientation="vertical", ticks=np.linspace(0, 1, 6, endpoint=True)
)

ax.legend(
    [
        line,
    ],
    [
        "Average Hole Energy",
    ],  # 图例标签
    fancybox=True,
    loc="upper right",
    framealpha=0.7,
    fontsize=9,
)

############################################################
# 设置图表范围和标签
############################################################
ax.set_xlim(0, namdTime)
# ax.set_ylim(-1.0, 1.0)

ax.set_xlabel("Time [fs]", fontsize="small", labelpad=5)
ax.set_ylabel("Energy [eV]", fontsize="small", labelpad=5)

########################################
plt.tight_layout(pad=0.2)
plt.savefig("kpoen.png", dpi=360)  # 保存高分辨率图像
