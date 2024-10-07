import numpy as np
import matplotlib.pyplot as plt

m1x, s1x, m1y, s1y = -2, 2, 0, 1

data1 = np.random.normal(m1x, s1x, 25)
data2 = np.random.normal(m1y, s1y, 25)

m2x, s2x, m2y, s2y = 2, 1, 2, 1

data3 = np.random.normal(m2x, s2x, 25)
data4 = np.random.normal(m2y, s2y, 25)

# m3x, s3x, m3y, s3y = 0, 1, -2, 1

# data5 = np.random.normal(m3x, s3x, 25)
# data6 = np.random.normal(m3y, s3y, 25)

plt.plot(data1, data2, 'o', color='black')
plt.plot(data3, data4, 'o', color='blue')
# plt.plot(data5, data6, 'o', color='red')
plt.style.use('fivethirtyeight')
plt.title('Two Clusters of Data')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
#Show lengend
plt.show()
plt.legend(['Cluster 1', 'Cluster 2'])

