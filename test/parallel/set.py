import os

for i in range(1, 25):
    os.mkdir(str(i))
    os.system("cp bifurcation* {0}/".format(i))
