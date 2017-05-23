import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style("ticks")
sns.set_palette("husl")

fig = plt.figure(1, figsize=(8,8))
fig.clf()
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)

for pext in [100, 10, 0, -10, -100]:
	p = np.loadtxt("single-artery-pext-%s_results/A1_P.last" % (pext))
	q = np.loadtxt("single-artery-pext-%s_results/A1_Q.last" % (pext))
	a = np.loadtxt("single-artery-pext-%s_results/A1_A.last" % (pext))
	u = np.loadtxt("single-artery-pext-%s_results/A1_u.last" % (pext))

	ax1.plot(p[:,0], p[:,3]/133.332, label=pext)
	ax2.plot(q[:,0], q[:,3]*1e6)
	ax3.plot(a[:,0], np.sqrt(a[:,3]/np.pi)*1e3)
	ax4.plot(u[:,0], u[:,3])

ax1.set_xlabel("time (s)")
ax2.set_xlabel("time (s)")
ax3.set_xlabel("time (s)")
ax4.set_xlabel("time (s)")

ax1.set_ylabel("pressure (mmHg)")
ax2.set_ylabel("flow (ml/s)")
ax3.set_ylabel("lumen radius (mm)")
ax4.set_ylabel("velocity (m/s)")

ax1.legend(loc="best")

plt.tight_layout()

sns.despine()

plt.draw()
plt.show()