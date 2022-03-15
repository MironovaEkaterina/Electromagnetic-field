import matplotlib.pyplot as plt
import numpy as np

c = 3e+10
n = 120
s1 = []
s2 = []
s3 = []
Ex = np.zeros( (n+1, n+1) )

f = open("result.txt", "r")
s=f.readlines()
for i in range(len(s)):
    pair=s[i].split(";")
    s2.append(int(pair[0]))
    s3.append(int(pair[1]))
    s1.append(float(pair[2]))
f.close()

for i in range(len(s2)):
    Ex[s2[i]][s3[i]]=s1[i]

#sctr = plt.matshow(Ex, cmap = 'Greys')
sctr = plt.scatter(s2,s3,c=s1, cmap = 'Greys')
plt.colorbar(sctr, format='%f')
plt.show()

 