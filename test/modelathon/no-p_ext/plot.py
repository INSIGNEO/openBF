import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

arteries = ["femoral-a","post_tibial-a","medial_plantar-a","lateral_plantar-a",
		"plantar_digital-a","ant_tibial-a","dorsalis-a","metatarsal-a","digital-a"]

veins = ["metatarsal-v","medial_plantar_II-v","medial_plantar_I-v","post_tibial_II-v",
		"post_tibial_I-v","popliteal-v","deep_peroneal-v","medial_marginal-v",
		"ant_tibial_III-v","ant_tibial_II-v","ant_tibial_I-v","great_saphenous-v"]

sns.set_style("ticks")
sns.set_palette("husl")

fig_idx = 1
for fldr in ["modelathon_results"]:

	fig = plt.figure(fig_idx, figsize=(8,8))
	fig.clf()
	ax1 = fig.add_subplot(221)
	ax2 = fig.add_subplot(222)
	ax3 = fig.add_subplot(223)
	ax4 = fig.add_subplot(224)

	for a in arteries:

		pA = np.loadtxt("%s/%s_P.last" % (fldr, a))
		qA = np.loadtxt("%s/%s_Q.last" % (fldr, a))
		aA = np.loadtxt("%s/%s_A.last" % (fldr, a))
		uA = np.loadtxt("%s/%s_u.last" % (fldr, a))

		ax1.plot(pA[:,0], pA[:,3]/133.332, label=a)
		ax2.plot(qA[:,0], qA[:,3]*1e6)
		ax3.plot(aA[:,0], np.sqrt(aA[:,3]/np.pi)*1e3)
		ax4.plot(uA[:,0], uA[:,3])

	for a in veins:
		pA = np.loadtxt("%s/%s_P.last" % (fldr, a))
		qA = np.loadtxt("%s/%s_Q.last" % (fldr, a))
		aA = np.loadtxt("%s/%s_A.last" % (fldr, a))
		uA = np.loadtxt("%s/%s_u.last" % (fldr, a))

		ax1.plot(pA[:,0], pA[:,3]/133.332, '--', label=a)
		ax2.plot(qA[:,0], qA[:,3]*1e6, '--')
		ax3.plot(aA[:,0], np.sqrt(aA[:,3]/np.pi)*1e3, '--')
		ax4.plot(uA[:,0], uA[:,3], '--')

	ax1.set_xlabel("time (s)")
	ax1.set_ylabel("pressure (mmHg)")
	ax1.legend(loc="best")

	ax2.set_xlabel("time (s)")
	ax2.set_ylabel("flow (ml/s)")

	plt.tight_layout()

	fig_idx += 1

plt.draw()
plt.show()
