import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


P = np.loadtxt("bifurcation_results/P_P.out")
Q = np.loadtxt("bifurcation_results/P_Q.out")
A = np.loadtxt("bifurcation_results/P_A.out")
R = np.sqrt(A[:,3]/np.pi)
R0 = 0.83183427183602578*1e-2 +0.000125
dR = R - R0

Rb = np.sqrt(A[:,-1]/np.pi)
R0b = 0.83183427183602578*1e-2 +0.000125
dRb = Rb - R0b

P1 = np.loadtxt("bifurcation_results/d1_P.out")
Q1 = np.loadtxt("bifurcation_results/d1_Q.out")
A1 = np.loadtxt("bifurcation_results/d1_A.out")
R1 = np.sqrt(A1[:,3]/np.pi)
R01 = 0.58844816698429576*1e-2 +0.00003
dR1 = R1 - R01


sns.set_style("white")
sns.set_context("notebook", font_scale=1, rc={"lines.linewidth": 2})
# plt.locator_params(nbins=5)
fig = plt.figure(1)
fig.clf()

ax1 = fig.add_subplot(331)
ax1.plot(P[:,0], P[:,3]/1000, 'k')
ax1.set_xlim(8.8, 9.9)
ax1.set_ylim(8, 18)
ax1.set_title("(a) Pressure")

pad = 5
ax1.annotate("aorta", xy=(0, 0.5), xytext=(-ax1.yaxis.labelpad - pad, 0),
                xycoords=ax1.yaxis.label, textcoords='offset points',
                size='large', ha='right', va='center')

plt.setp( ax1.get_xticklabels(), visible=False)
ax1.set_ylabel("$P$  $($kPa$)$")
ax1.legend(loc="best")

ax2 = fig.add_subplot(332)
ax2.plot(Q[:,0], Q[:,3]*1e6, 'k')
ax2.set_xlim(8.8, 9.9)
ax2.set_ylim(-24, 80)
ax2.set_title("(b) Flow")
plt.setp( ax2.get_xticklabels(), visible=False)
ax2.set_ylabel("$Q$  $($ml$\cdot$s$^{-1})$")

ax3 = fig.add_subplot(333)
ax3.plot(A[:,0], dR*1e3, 'k')
ax3.set_xlim(8.8, 9.9)
ax3.set_ylim(-0.05, 0.9)
ax3.set_title("(c) Radius Change")
plt.setp( ax3.get_xticklabels(), visible=False)
ax3.set_ylabel("$\Delta r$  $($mm$)$")

ax4 = fig.add_subplot(334)
ax4.plot(P[:,0], P[:,-1]/1000, 'k')
ax4.set_xlim(8.8, 9.9)
ax4.set_ylim(7, 19)

ax4.annotate("bifurcation", xy=(0, 0.5), xytext=(-ax4.yaxis.labelpad - pad, 0),
                xycoords=ax4.yaxis.label, textcoords='offset points',
                size='large', ha='right', va='center')

plt.setp( ax4.get_xticklabels(), visible=False)
ax4.set_ylabel("$P$  $($kPa$)$")

ax5 = fig.add_subplot(335)
ax5.plot(Q[:,0], Q[:,-1]*1e6, 'k')
ax5.set_xlim(8.8, 9.9)
ax5.set_ylim(-20, 70)
plt.setp( ax5.get_xticklabels(), visible=False)
ax5.set_ylabel("$Q$  $($ml$\cdot$s$^{-1})$")

ax6 = fig.add_subplot(336)
ax6.plot(A[:,0], dRb*1e3, 'k')
ax6.set_xlim(8.8, 9.9)
ax6.set_ylim(-0.1, 0.9)
plt.setp( ax6.get_xticklabels(), visible=False)
ax6.set_ylabel("$\Delta r$  $($mm$)$")

ax7 = fig.add_subplot(337)
ax7.plot(P1[:,0], P1[:,3]/1000, 'k')
ax7.set_xlim(8.8, 9.9)
ax7.set_ylim(7, 18)

ax7.annotate("iliac", xy=(0, 0.5), xytext=(-ax7.yaxis.labelpad - pad, 0),
                xycoords=ax7.yaxis.label, textcoords='offset points',
                size='large', ha='right', va='center')

ax7.set_ylabel("$P$  $($kPa$)$")
ax7.set_xlabel("$t$  $($s$)$")

ax8 = fig.add_subplot(338)
ax8.plot(Q1[:,0], Q1[:,3]*1e6, 'k')
ax8.set_xlim(8.8, 9.9)
ax8.set_ylim(-10, 30)
ax8.set_ylabel("$Q$  $($ml$\cdot$s$^{-1})$")
ax8.set_xlabel("$t$  $($s$)$")

ax9 = fig.add_subplot(339)
ax9.plot(A1[:,0], dR1*1e3, 'k')
ax9.set_xlim(8.8, 9.9)
ax9.set_ylim(-0.07, 0.46)
ax9.set_ylabel("$\Delta r$  $($mm$)$")
ax9.set_xlabel("$t$  $($s$)$")

sns.despine()

'''
def quin(t):
	return 1e-6 * np.exp(-1e4*(t-0.05)**2)

t = np.linspace(0., 2., 1000)
quint = quin(t)

quint += 1e-10

data = np.vstack((t, quint))

np.savetxt("tutorial_inlet.dat", np.transpose(data))

plt.figure(2)
plt.clf()
plt.plot(t, quint)


PI=3.141592653589793
period T=1.1
time t
inlet flow rate in ml/s

T = 1.1
t = np.linspace(0., T, 100)
PI = np.pi
inflow=10e5*(7.9853e-06+2.6617e-05*np.sin(2*PI*t/T+0.29498)+2.3616e-05*np.sin(4*PI*t/T-1.1403)-1.9016e-05*np.sin(6*PI*t/T+0.40435)-8.5899e-06*np.sin(8*PI*t/T-1.1892)-2.436e-06*np.sin(10*PI*t/T-1.4918)+1.4905e-06*np.sin(12*PI*t/T+1.0536)+1.3581e-06*np.sin(14*PI*t/T-0.47666)-6.3031e-07*np.sin(16*PI*t/T+0.93768)-4.5335e-07*np.sin(18*PI*t/T-0.79472)-4.5184e-07*np.sin(20*PI*t/T-1.4095)-5.6583e-07*np.sin(22*PI*t/T-1.3629)+4.9522e-07*np.sin(24*PI*t/T+0.52495)+1.3049e-07*np.sin(26*PI*t/T-0.97261)-4.1072e-08*np.sin(28*PI*t/T-0.15685)-2.4182e-07*np.sin(30*PI*t/T-1.4052)-6.6217e-08*np.sin(32*PI*t/T-1.3785)-1.5511e-07*np.sin(34*PI*t/T-1.2927)+2.2149e-07*np.sin(36*PI*t/T+0.68178)+6.7621e-08*np.sin(38*PI*t/T-0.98825)+1.0973e-07*np.sin(40*PI*t/T+1.4327)-2.5559e-08*np.sin(42*PI*t/T-1.2372)-3.5079e-08*np.sin(44*PI*t/T+0.2328))
data = np.vstack((t, inflow*1e-6))

np.savetxt("tutorial_inlet.dat", np.transpose(data))


plt.figure(2)
plt.clf()
plt.plot(t, inflow)
'''

plt.draw()
plt.show()
