import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


Q1 = np.loadtxt("aspirator_results/v1_Q.last")
Q2 = np.loadtxt("aspirator_results/v2_Q.last")
Q3 = np.loadtxt("aspirator_results/v3_Q.last")


sns.set_style("white")
sns.set_context("notebook", font_scale=1, rc={"lines.linewidth": 2})
# plt.locator_params(nbins=5)
fig = plt.figure(1)
fig.clf()

ax = fig.add_subplot(111)

ax.plot(Q1[:,0], Q1[:,1], 'r', label="Q1")
ax.plot(Q2[:,0], Q2[:,1], 'k', linestyle="--", label="Q2")
ax.plot(Q3[:,0], Q3[:,-1], 'c', label="Q3")
plt.legend()
plt.tight_layout()
plt.draw()
plt.show()

print(Q1[0, 1])
print(Q2[0, 1])
print(Q3[-1,-1])
