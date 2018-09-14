import matplotlib.pyplot as plt
import numpy as np

plt.figure(1)

for t in range(1,5):
    print(t)
    plt.subplot(2,2,t)
    I=np.arange(5)
    plt.plot(I, t*I, c='b')
    plt.plot(range(t), c='r')
    print(t*I)
plt.show()