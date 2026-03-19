import numpy as np
from math import sqrt

def Read_cube_data(E_unit, cube_file):
    f=open(cube_file, 'r')
    f.readline()
    f.readline()

    atoms=f.readline().split()          
    for i in range(4):
        atoms[i]=int(float(atoms[i]))

    grid=[]
    for i in range(3):
        line=f.readline().split()
        line[0]=int(line[0])
        for j in range(3):
            line[j+1]=float(line[j+1])*0.529177249
        grid.append(line)

    grid_size=[sqrt(grid[0][1]**2+grid[0][2]**2+grid[0][3]**2)
               , sqrt(grid[1][1]**2+grid[1][2]**2+grid[1][3]**2)
               , sqrt(grid[2][1]**2+grid[2][2]**2+grid[2][3]**2)
                ]                       #Grid size (distance between points in each direction) in Angstrom
    f.close()

    """read potential data"""
    if E_unit==1:                       #kJ/mol
        conversion=1.0  
    elif E_unit==2:                     #kcal/mol   
        conversion=4.184
    elif E_unit==3:                     #Ry
        conversion=1312.7497
    elif E_unit==4:                     #eV
        conversion=96.4853
    elif E_unit==5:                     #Hartree
        conversion=627.5096

    cube_data=np.loadtxt(cube_file, skiprows=6+atoms[0])
    shift_cube_data=(cube_data-np.min(cube_data))*conversion

    data_matrix=np.zeros((grid[0][0], grid[1][0], grid[2][0]))
    line=0
    for a in range(grid[0][0]):
        for b in range(grid[1][0]):
            for c in range(grid[2][0]):
                data_matrix[a][b][c]=shift_cube_data[line]
                line+=1
    return atoms, grid, grid_size, data_matrix