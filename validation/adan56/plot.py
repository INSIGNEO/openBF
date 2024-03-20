import matplotlib.pyplot as plt
import numpy as np

plt.style.use("../plot_style.txt")

def plot_arteries(arteries, fig_num, j):
    fig = plt.figure(fig_num+2)
    fig.clf()

    for i, artery in enumerate(arteries):
        ax1 = fig.add_subplot(3, 2, i * 2 + 1)
        ax2 = fig.add_subplot(3, 2, i * 2 + 2)

        P = np.loadtxt(f"adan56_results/{artery}_P.last")
        Q = np.loadtxt(f"adan56_results/{artery}_Q.last")

        ax1.plot(ta, P_ref[:, j] * 1e3 / 133.332, color="dimgrey")
        ax2.plot(ta, Q_ref[:, j], label="ref", color="dimgrey")
        ax1.plot(np.linspace(0, 1, P.shape[0]), P[:, 3] / 133.332, label="openBF", linestyle="--", color="k")
        ax2.plot(np.linspace(0, 1, P.shape[0]), Q[:, 3] * 1e6, label="openBF", linestyle="--", color="k")
        
        ax1.set_title(artery)
        plt.legend()
        j += 1

    plt.tight_layout()


P_ref = np.loadtxt("adan56_ref/FVM_P_adan.txt")
Q_ref = np.loadtxt("adan56_ref/FVM_Q_adan.txt")

t_adan = P_ref[:,0]
t0_index = np.argwhere(t_adan == t_adan[-1]-1.)
P_ref = P_ref[t0_index[0][0]:, :]
Q_ref = Q_ref[t0_index[0][0]:, :]
ta = np.linspace(0, 1, len(P_ref[:,0]))


fig = plt.figure(4, figsize=(10,12))
fig.clf()

mapping = dict(zip(["aortic_arch_I", "thoracic_aorta_III", "abdominal_aorta_V",
    "common_carotid_R", "renal_R", "common_iliac_R",
    "internal_carotid_R", "radial_R", "internal_iliac_R",
    "posterior_interosseous_R", "femoral_R_I", "anterior_tibial_R",],
    ["aortic_arch", "thoracic_aorta", "abdominal_aorta", "common_carotid",
    "renal", "common_iliac", "internal_carotid", "radial", "internal_iliac",
    "posterior_interosseous", "femoral", "anterior_tibial"]))

for i, artery in enumerate(["aortic_arch_I", "thoracic_aorta_III", "abdominal_aorta_V",
    "common_carotid_R", "renal_R", "common_iliac_R",
    "internal_carotid_R", "radial_R", "internal_iliac_R",
    "posterior_interosseous_R", "femoral_R_I", "anterior_tibial_R",]):

    axP = fig.add_subplot(6, 4, (2*i)+1)
    axQ = fig.add_subplot(6, 4, (2*i)+2)

    axP.plot(ta, P_ref[:, i+1] * 1e3 / 133.332, color="dimgrey")
    axQ.plot(ta, Q_ref[:, i+1], label="Boileau2015-FV1D", color="dimgrey")

    P = np.loadtxt(f"../../adan56avg_results/{artery}_P.last")
    Q = np.loadtxt(f"../../adan56avg_results/{artery}_Q.last")
    axP.plot(np.linspace(0, 1, P.shape[0]), P[:, 3] / 133.332, label="openBF", linestyle="--", color="k")
    axQ.plot(np.linspace(0, 1, P.shape[0]), Q[:, 3] * 1e6, label="openBF", linestyle="--", color="k")

    axP.set_ylabel("$P$ (kPa)")
    axQ.set_ylabel("$Q$ (ml$\cdot$s$^{-1}$)")

    for ax in [axP, axQ]:
        ax.set_xlim(0,1)
        ax.set_xlabel("$t\,/\,T_c$")

    axP.set_title(mapping[artery].replace("_", " ").title())

    if i == 1:
        legend_handle, legend_labels = ax.get_legend_handles_labels()

fig.legend(
        legend_handle,
        legend_labels,
        loc="lower center",
        bbox_to_anchor=(0.5, 0),
        ncol=3,
    )
plt.tight_layout(rect=[0, 0.03, 1, 1])
plt.savefig("adan56.eps", dpi=300, format="eps")

