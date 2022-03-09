import matplotlib.pyplot as plt

points_Ey=[]
points_x=[]
points_Ey_real=[]

f = open('result.txt', 'r')
s=f.readlines()
for i in range(len(s)):
    pair=s[i].split(";")
    points_Ey.append(float(pair[0]))
    points_x.append(float(pair[1]))

f2 = open('result_real.txt', 'r')
s=f2.readlines()
for i in range(len(s)):
    pair=s[i].split(";")
    points_Ey_real.append(float(pair[0]))

plt.plot(points_x,points_Ey_real,'bo')
plt.plot(points_x,points_Ey,'ro')
plt.show()
f.close()