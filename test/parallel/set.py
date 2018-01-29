import os

for i in range(1, 25):
    os.mkdir(str(i))
    os.system("cp bifurcation* {0}/".format(i))
    os.system("cp run_parallel_openBF.sh {0}/".format(i))
