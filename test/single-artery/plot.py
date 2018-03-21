import numpy as np
import matplotlib.pyplot as plt

pA = np.loadtxt("single-artery_results/A1_P.out")

fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)

ax.plot(pA[:,0], pA[:,-1]/133.332, label="artery")

ax.set_xlabel("time (s)")
ax.set_ylabel("pressure (mmHg)")

plt.legend(loc="best")

plt.draw()
plt.show()
