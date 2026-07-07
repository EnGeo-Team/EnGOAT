import numpy as np
import logging
import data                                                                                             #Python file containing all the data-containing matrices used in the simulation
logger = logging.getLogger(__name__)

def get_topological_descriptors(N_levels):

    """Surface area"""
    a=np.array([data.grid[0][1], data.grid[0][2], data.grid[0][3]])                                     #Unit cell vectors [Å]
    b=np.array([data.grid[1][1], data.grid[1][2], data.grid[1][3]])
    c=np.array([data.grid[2][1], data.grid[2][2], data.grid[2][3]])

    ab=np.linalg.norm(np.cross(a, b))                                                                   #Surface areas of grid point sides [Å**2]
    ac=np.linalg.norm(np.cross(a, c))
    bc=np.linalg.norm(np.cross(b, c))

    border_points=np.argwhere(data.Cluster_matrix.Levels==0)                                            #Points above E_cutoff

    Total_area=0
    for point in border_points:
        x_neighbours, y_neighbours, z_neighbours=PBC3D_neighbours(point)
        for neighbour_point in x_neighbours:
            if data.Cluster_matrix.Levels[neighbour_point]>0:
                Total_area+=bc
        for neighbour_point in y_neighbours:
            if data.Cluster_matrix.Levels[neighbour_point]>0:
                Total_area+=ac
        for neighbour_point in z_neighbours:
            if data.Cluster_matrix.Levels[neighbour_point]>0:
                Total_area+=ab
    
    Accessible_area=0
    for tunnel_system in data.Tunnel_system.Tunnel_system_list:
        Tunnel_area=0
        for point in border_points:
            x_neighbours, y_neighbours, z_neighbours=PBC3D_neighbours(point)
            for neighbour_point in x_neighbours:
                if data.Cluster_matrix.Levels[neighbour_point]>0 and data.Cluster_matrix.IDs[neighbour_point] in tunnel_system.cluster_family:
                    Tunnel_area+=bc
            for neighbour_point in y_neighbours:
                if data.Cluster_matrix.Levels[neighbour_point]>0 and data.Cluster_matrix.IDs[neighbour_point] in tunnel_system.cluster_family:
                    Tunnel_area+=ac
            for neighbour_point in z_neighbours:
                if data.Cluster_matrix.Levels[neighbour_point]>0 and data.Cluster_matrix.IDs[neighbour_point] in tunnel_system.cluster_family:
                    Tunnel_area+=ab
        tunnel_system.surface_area=Tunnel_area
        Accessible_area+=Tunnel_area

    """Volume and histograms"""
    abc=abs(np.dot(a, np.cross(b, c)))                                                              #Volume of a grid point  [Å**3]

    Accessible_points=np.argwhere(data.Cluster_matrix.Levels!=0)
    Total_volume=len(Accessible_points)*abc
    Total_volume_fraction=len(Accessible_points)/(data.grid[0][0]*data.grid[1][0]*data.grid[2][0])

    Accessible_volume=0
    for tunnel_system in data.Tunnel_system.Tunnel_system_list:
        N_grid_points=0                                                                             #Calculate the volume fraction the tunnel system occupies within the unit cell
        for cluster_ID in tunnel_system.cluster_family:
            N_grid_points+=np.count_nonzero(data.Cluster_matrix.IDs==cluster_ID)
        tunnel_system.volume=N_grid_points*abc
        tunnel_system.V_fraction=N_grid_points/(data.grid[0][0]*data.grid[1][0]*data.grid[2][0])
        Accessible_volume+=tunnel_system.volume

        histogram=np.zeros(N_levels)
        tunnel_system_mask=(data.Tunnel_matrix.IDs==tunnel_system.ID)
        for level in range(1, N_levels+1):
            histogram[level-1]=np.sum((data.level_matrix==level) & tunnel_system_mask)
        tunnel_system.histogram=histogram
        
    print(Total_volume)
    return Total_area, Accessible_area, Total_volume, Total_volume_fraction, Accessible_volume

def Get_basin_area(basin_ID, Basin_matrix):

    # Unit cell vectors
    a = np.array([data.grid[0][1], data.grid[0][2], data.grid[0][3]])
    b = np.array([data.grid[1][1], data.grid[1][2], data.grid[1][3]])
    c = np.array([data.grid[2][1], data.grid[2][2], data.grid[2][3]])

    # Face areas
    ab = np.linalg.norm(np.cross(a, b))   # z-faces
    ac = np.linalg.norm(np.cross(a, c))   # y-faces
    bc = np.linalg.norm(np.cross(b, c))   # x-faces

    basin_points = np.argwhere(Basin_matrix == basin_ID)

    area = 0.0

    for point in basin_points:

        x_neighbours, y_neighbours, z_neighbours = PBC3D_neighbours(point)

        for neighbour in x_neighbours:
            if Basin_matrix[neighbour] != basin_ID:
                area += bc

        for neighbour in y_neighbours:
            if Basin_matrix[neighbour] != basin_ID:
                area += ac

        for neighbour in z_neighbours:
            if Basin_matrix[neighbour] != basin_ID:
                area += ab

    return area


def PBC3D_neighbours(point):                                                                            #Subroutine that finds PBC neighbours of a given point
    i, j, k = point
    Nx, Ny, Nz = data.grid[0][0], data.grid[1][0], data.grid[2][0]
    x_neighbours = [((i + 1) % Nx, j, k), ((i - 1) % Nx, j, k)]
    y_neighbours = [(i, (j + 1) % Ny, k), (i, (j - 1) % Ny, k)]
    z_neighbours = [(i, j, (k + 1) % Nz), (i, j, (k - 1) % Nz)]
    return x_neighbours, y_neighbours, z_neighbours