import numpy as np
import matplotlib.pyplot as plt

plt.style.use("plot_style.txt")


def plot_single_artery_validation(
    Tt,
    B3D,
    FV,
    P,
    Q,
    A,
    xlims,
    shift,
    figname,
    figformat,
    fignum,
    title,
):
    T = P[:, 0]
    t0 = np.argwhere(T > T[-1] - Tt)
    t0 = t0[0][0]

    P = P[t0:, :]
    Q = Q[t0:, :]
    A = A[t0:, :]

    T = np.linspace(0.0, Tt, P.shape[0])

    R = np.sqrt(A / np.pi)
    R -= np.min(R)
    dP = P[:, 1] - P[:, -1]

    t0FV = np.argwhere(FV[:, 0] > FV[-1, 0] - Tt - 0.05)
    t0FV = t0FV[0][0]
    tFV = np.linspace(0.0, Tt, len(FV[t0FV:-50, 0]))

    t3d = np.linspace(19, 19.955, B3D.shape[0])

    waves = [P[:, 3] * 1e-3, Q[:, 3] * 1e6, R[:, 3] * 1e3, dP * 1e-3]
    titles = [
        "(a) Pressure",
        "(b) Flow",
        "(c) Radius Change",
        "(d) Pressure Difference",
    ]
    ylabels = [
        "$P$ (kPa)",
        "$Q$ (ml$\cdot$s$^{-1}$)",
        "$\Delta r$ (mm)",
        "$\Delta P$ (kPa)",
    ]

    fig = plt.figure(fignum, figsize=(12, 3))
    fig.clf()

    for i in range(4):
        ax = fig.add_subplot(1, 4, i + 1)

        if figname == "cca":
            ax.plot(tFV / T[-1], FV[t0FV:-50, i + 1], color="slategrey", label="1D")
            ax.plot(
                B3D[:, 0] / T[-1],
                B3D[:, i + 1],
                linestyle=":",
                color="gray",
                label="3D",
            )
            ax.plot(
                T[:-shift] / T[-1],
                waves[i][shift:],
                "k",
                linestyle="--",
                label="openBF",
            )
            ax.plot(T[-shift:] / T[-1], waves[i][:shift], "k", linestyle="--")
        else:
            ax.plot(
                FV[:, 0] / FV[len(tFV), 0], FV[:, i + 1], color="slategrey", label="1D"
            )
            ax.plot(
                B3D[:, 0] / B3D[-1, 0],
                B3D[:, i + 1],
                linestyle=":",
                color="gray",
                label="3D",
            )
            ax.plot(T / T[-1], waves[i], "k", linestyle="--", label="openBF")

        ax.set_xlim(0, 1.0)
        ax.set_ylim(xlims[i])
        ax.set_xlabel("$t\,/\,T_c$")

        if i == 1:
            legend_handle, legend_labels = ax.get_legend_handles_labels()
        ax.set_ylabel(ylabels[i])
        # ax.set_title(titles[i])

    fig.legend(
        legend_handle,
        legend_labels,
        loc="lower center",
        bbox_to_anchor=(0.5, 0),
        ncol=3,
    )
    # plt.suptitle(title)

    plt.tight_layout(rect=[0, 0.09, 1, 1])

    plt.savefig(f"{figname}.{figformat}", dpi=300, format=figformat)


def plot_iliac_bifurcation_validation(fignum, figname, figformat):
    B3Da = np.loadtxt("ibif/ibif_ref/B3D_AortaMid.txt")
    B3Di = np.loadtxt("ibif/ibif_ref/B3D_IliacMid.txt")
    B3Dj = np.loadtxt("ibif/ibif_ref/B3D_Junction.txt")

    P = np.loadtxt("ibif/ibif_results/parent_P.last")
    Q = np.loadtxt("ibif/ibif_results/parent_Q.last")
    A = np.loadtxt("ibif/ibif_results/parent_A.last")
    R = np.sqrt(A[:, 3] / np.pi)
    R0 = 0.83183427183602578 * 1e-2 + 0.000125
    dR = R - np.min(R)  # R0

    Rb = np.sqrt(A[:, -1] / np.pi)
    R0b = 0.83183427183602578 * 1e-2 + 0.000125
    dRb = Rb - R0b

    P1 = np.loadtxt("ibif/ibif_results/d1_P.last")
    Q1 = np.loadtxt("ibif/ibif_results/d1_Q.last")
    A1 = np.loadtxt("ibif/ibif_results/d1_A.last")
    R1 = np.sqrt(A1[:, 3] / np.pi)
    R01 = 0.58844816698429576 * 1e-2 + 0.00003
    dR1 = R1 - np.min(R1)  # R01

    FV = np.loadtxt("ibif/ibif_ref/FV_aorta-bif.dat")
    FVb = np.loadtxt("ibif/ibif_ref/FV_aorta_bif-bif.dat")
    FVi = np.loadtxt("ibif/ibif_ref/FV_iliac-bif.dat")

    fig = plt.figure(fignum, figsize=(9, 6))
    fig.clf()

    t3d = np.linspace(0.0, 1.0, B3Da.shape[0])

    qs1d = [FV, FVb, FVi]
    qs3d = [B3Da, B3Dj, B3Di]
    qsobf = [
        [P[:, 3] / 1000, Q[:, 3] * 1e6, dR * 1e3],
        [P[:, -1] / 1000, Q[:, -1] * 1e6, dRb * 1e3],
        [P1[:, 3] / 1000, Q1[:, 3] * 1e6, dR1 * 1e3],
    ]
    ylims = [
        [(8, 18), (-24, 80), (-0.05, 0.9)],
        [(7, 19), (-20, 70), (-0.1, 0.9)],
        [(7, 18), (-10, 30), (-0.07, 0.46)],
    ]
    for i, location in enumerate(["aorta", "junction", "iliac"]):
        for j, q in enumerate(
            [
                "$P$ (kPa)",
                "$Q$ (ml$\cdot$s$^{-1}$)",
                "$\Delta r$ (mm)",
            ]
        ):
            ax = fig.add_subplot(3, 3, 3 * i + j + 1)
            ax.plot(
                (qs1d[i][8800:9900, 0] - 8.8) / 1.1,
                qs1d[i][8800:9900, j + 1],
                color="slategrey",
                label="1D",
            )
            ax.plot(t3d, qs3d[i][:, j + 1], color="gray", linestyle=":", label="3D")
            ax.plot(
                (P[:, 0] - P[0, 0]) / (P[-1, 0] - P[0, 0]),
                qsobf[i][j],
                "k",
                linestyle="--",
                label="openBF",
            )

            ax.set_xlim(0, 1)
            ax.set_ylim(ylims[i][j])
            ax.set_ylabel(q)

            if i == 1 and j == 1:
                legend_handle, legend_labels = ax.get_legend_handles_labels()

            ax.set_xlabel("$t\,/\,T_c$")
            # ax.set_ylabel(ylabels[i])
            if j == 0:
                ax.annotate(
                    location,
                    xy=(0, 0.5),
                    xytext=(-ax.yaxis.labelpad - 5, 0),
                    xycoords=ax.yaxis.label,
                    textcoords="offset points",
                    size="large",
                    ha="right",
                    va="center",
                )

    fig.legend(
        legend_handle,
        legend_labels,
        loc="lower center",
        bbox_to_anchor=(0.5, 0),
        ncol=3,
    )

    plt.tight_layout(rect=[0, 0.09, 1, 1])

    plt.savefig(f"{figname}.{figformat}", dpi=300, format=figformat)


# uta
B3D = np.loadtxt("uta/uta_ref/B3D_aortic.txt")
FV = np.loadtxt("uta/uta_ref/FV_aortic.dat")
P = np.loadtxt("uta/uta_results/upper_thoracic_aorta_P.last")
Q = np.loadtxt("uta/uta_results/upper_thoracic_aorta_Q.last")
A = np.loadtxt("uta/uta_results/upper_thoracic_aorta_A.last")
Tt = 9.549999999999999600e-01
xlims = [(8.5, 18), (-80, 520), (-0.1, 1.7), (-3, 5)]
shift = 0
figname = "uta"
figformat = "eps"
fignum = 1
title = "Upper Thoracic Aorta"

plot_single_artery_validation(
    Tt, B3D, FV, P, Q, A, xlims, shift, figname, figformat, fignum, title
)

# cca
Tt = 1.1
B3D = np.loadtxt("cca/cca_ref/B3D_commoncarotid.txt")
FV = np.loadtxt("cca/cca_ref/FV_carotid.dat")
P = np.loadtxt("cca/cca_results/common_carotid_artery_P.out")
Q = np.loadtxt("cca/cca_results/common_carotid_artery_Q.out")
A = np.loadtxt("cca/cca_results/common_carotid_artery_A.out")
xlims = [(10, 17.5), (3, 13), (-0.01, 0.2), (-0.2, 0.8)]
shift = 5
figname = "cca"
fignum = 2
title = "Common Carotid Artery"

plot_single_artery_validation(
    Tt, B3D, FV, P, Q, A, xlims, shift, figname, figformat, fignum, title
)

# ibif
fignum = 3
figname = "ibif"
plot_iliac_bifurcation_validation(fignum, figname, figformat)
