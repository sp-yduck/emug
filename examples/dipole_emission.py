import numpy as np
import matplotlib.pyplot as plt

from emug.dipole import ElectricDipole
    

if __name__ == "__main__":
    dipole = ElectricDipole(moment=1e-4, freq=2e+9, dipole_direction='x')
    n = 200
    x = np.array([[i]*n for i in range(-n//2,n//2)])*0.5*1e-3
    y = np.array([[i for i in range(-n//2,n//2)] for k in range(n)])*0.5*1e-3

    name_tags = ['Hx [A/m]', 'Hy [A/m]', 'Hz [A/m]', 'Ex [V/m]', 'Ey [V/m]', 'Ez [V/m]']
    tags = [3,4,5,0,1,2]

    plt.figure(figsize=(16, 9))
    for k, (tag, name_tag) in enumerate(zip(tags, name_tags)):
        plt.subplot(2,3,k+1)
        mat1 = dipole.em(x,y,-0.050,0)[tag]
        plt.colorbar(plt.imshow(np.abs(mat1)))
        plt.title(f"{name_tag}", fontsize=18)
        plt.xticks(np.arange(0,200, 50), np.arange(0,100,25))
        plt.yticks(np.arange(0,200, 50), np.arange(0,100,25))
        plt.xlabel('y axis [mm]')
        plt.ylabel('x axis [mm]')
        plt.scatter(100,100, color='red')
    plt.show()