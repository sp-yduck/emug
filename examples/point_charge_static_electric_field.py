import numpy as np
import matplotlib.pyplot as plt

from emug.charge import PointCharge


if __name__ == "__main__":
    p = PointCharge(1e-12)

    n = 200
    x = np.array([[i]*n for i in range(-n//2,n//2)])*0.5*1e-3
    y = np.array([[i for i in range(-n//2,n//2)] for k in range(n)])*0.5*1e-3

    name_tags = ['Ex [N/C]', 'Ey [N/C]', 'Ez [N/C]']
    tags = [0,1,2]

    plt.figure(figsize=(16, 4))
    for k, (tag, name_tag) in enumerate(zip(tags, name_tags)):
        plt.subplot(1,3,k+1)
        mat1 = p.e(x,y,0.050)[tag]
        plt.colorbar(plt.imshow(mat1))
        plt.title(f"{name_tag}", fontsize=18)
        plt.xticks(np.arange(0,200, 50), np.arange(0,100,25))
        plt.yticks(np.arange(0,200, 50), np.arange(0,100,25))
        plt.xlabel('y axis [mm]')
        plt.ylabel('x axis [mm]')
        plt.scatter(100,100, color='red')
    plt.show()