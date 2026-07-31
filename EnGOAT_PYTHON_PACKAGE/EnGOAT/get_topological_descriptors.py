import numpy as np
import logging
from scipy import constants
from collections import defaultdict
import data                                                                                             #Python file containing all the data-containing matrices used in the simulation
logger = logging.getLogger(__name__)

def get_topological_descriptors(N_steps, E_step, T_list):

    V_UC = abs(np.linalg.det(data.grid["vectors"]))*np.prod(data.grid["shape"])                         #Unit cell volume

    for cluster in data.Cluster.Cluster_list:
        if cluster.active:
            cluster.A = Get_area(cluster.ID, data.Cluster_matrix.IDs)
            cluster.A_rel = cluster.A/V_UC
            cluster.V = Get_volume(cluster.ID, data.Cluster_matrix.IDs)
            cluster.V_rel = cluster.V/V_UC

            mask = (data.Cluster_matrix.IDs == cluster.ID)
            cluster.V_B = Get_Boltzmann_volume(mask, T_list)
            cluster.histogram = Get_histogram(mask, E_step, N_steps)

    for transition_state in data.Transition_state.Transition_state_list:
        transition_state.A = Get_area(transition_state.ID, data.TS_matrix.IDs)
        transition_state.A_rel = transition_state.A/V_UC
        transition_state.V = Get_volume(transition_state.ID, data.TS_matrix.IDs)
        transition_state.V_rel = transition_state.V/V_UC

        mask = (data.TS_matrix.IDs == transition_state.ID)
        transition_state.V_B = Get_Boltzmann_volume(mask, T_list)
        transition_state.histogram = Get_histogram(mask, E_step, N_steps)

    V_tot = 0
    V_acc = 0
    A_tot = 0
    A_acc = 0

    for tunnel_system in data.Tunnel_system.Tunnel_system_list:

        cluster_family = tunnel_system.cluster_family
        TS_list = tunnel_system.transition_state_list

        V=0
        V_B = {T: 0.0 for T in T_list}
        
        histogram = defaultdict(float)

        for cluster in data.Cluster.Cluster_list:
            if cluster.ID in cluster_family:
                V += cluster.V
                for T in T_list:
                    V_B[T] += cluster.V_B[T]

                for edge, count in cluster.histogram.items():
                    histogram[edge] += count
        
        for transition_state in data.Transition_state.Transition_state_list:
            if transition_state.ID in TS_list:
                V += transition_state.V
                for T in T_list:
                    V_B[T] += transition_state.V_B[T]

                for edge, count in transition_state.histogram.items():
                    histogram[edge] += count
        
        A = Get_group_area(list(cluster_family), list(TS_list), data.Cluster_matrix.IDs, data.TS_matrix.IDs)

        tunnel_system.A = A
        tunnel_system.V = V
        tunnel_system.A_rel = A/V_UC
        tunnel_system.V_rel = V/V_UC
        tunnel_system.V_B = V_B
        tunnel_system.histogram = histogram

        V_tot += V
        V_acc += V
        A_tot += A
        A_acc += A


    for isolated_group in data.Isolated_group.Isolated_group_list:
        cluster_family = isolated_group.cluster_family
        TS_list = isolated_group.transition_state_list

        V=0
        V_B = {T: 0.0 for T in T_list}
        
        histogram = defaultdict(float)

        for cluster in data.Cluster.Cluster_list:
            if cluster.ID in cluster_family:
                V += cluster.V
                for T in T_list:
                    V_B[T] += cluster.V_B[T]

                for edge, count in cluster.histogram.items():
                    histogram[edge] += count
        
        for transition_state in data.Transition_state.Transition_state_list:
            if transition_state.ID in TS_list:
                V += transition_state.V
                for T in T_list:
                    V_B[T] += transition_state.V_B[T]

                for edge, count in transition_state.histogram.items():
                    histogram[edge] += count

        A = Get_group_area(list(cluster_family), list(TS_list), data.Cluster_matrix.IDs, data.TS_matrix.IDs)
        isolated_group.A = A
        isolated_group.V = V
        isolated_group.A_rel = A/V_UC
        isolated_group.V_rel = V/V_UC
        isolated_group.V_B = V_B
        isolated_group.histogram = histogram

        V_tot += V
        A_tot += A

    return A_tot, A_acc, V_tot, V_acc, V_UC



def Get_Boltzmann_volume(mask, T_list):

    V_B = {}
    for T in T_list:
        Beta=1/(constants.R*T) 
        V_B[T] = np.sum(np.exp(-1000*Beta*(data.Energy_matrix[mask])))
    return V_B
                                    

def Get_histogram(mask, E_step, N_steps):
    voxel_V = abs(np.linalg.det(data.grid["vectors"]))

    levels = data.level_matrix[mask]

    # Always create N_steps bins, including empty ones
    counts = np.bincount(levels, minlength=N_steps)

    volumes = counts * voxel_V

    # Generate N_steps energy edges
    edges = np.arange(N_steps) * E_step

    return {
        float(edge): float(volume)
        for edge, volume in zip(edges[1:], volumes[1:])
    }  # skip the first bin


def Get_volume(ID, ID_matrix):
    voxel_V = abs(np.linalg.det(data.grid["vectors"]))
    N_points = np.sum(ID_matrix == ID)
    V = N_points*voxel_V
    return V

def Get_area(ID, ID_matrix):

    # Unit cell vectors
    a=data.grid["vectors"][0, :]
    b=data.grid["vectors"][1, :]
    c=data.grid["vectors"][2, :]

    # Face areas
    ab = np.linalg.norm(np.cross(a, b))   # z-faces
    ac = np.linalg.norm(np.cross(a, c))   # y-faces
    bc = np.linalg.norm(np.cross(b, c))   # x-faces

    points = np.argwhere(ID_matrix == ID)

    area = 0.0

    for point in points:

        x_neighbours, y_neighbours, z_neighbours = PBC3D_neighbours(point)

        for neighbour in x_neighbours:
            if ID_matrix[neighbour] != ID:
                area += bc

        for neighbour in y_neighbours:
            if ID_matrix[neighbour] != ID:
                area += ac

        for neighbour in z_neighbours:
            if ID_matrix[neighbour] != ID:
                area += ab

    return area

def Get_group_area(B_IDs, TS_IDs, Basin_matrix, TS_matrix):

    # Unit cell vectors
    a=data.grid["vectors"][0, :]
    b=data.grid["vectors"][1, :]
    c=data.grid["vectors"][2, :]

    # Face areas
    ab = np.linalg.norm(np.cross(a, b))   # z-faces
    ac = np.linalg.norm(np.cross(a, c))   # y-faces
    bc = np.linalg.norm(np.cross(b, c))   # x-faces

    basin_points = np.argwhere(np.isin(Basin_matrix, B_IDs))
    ts_points = np.argwhere(np.isin(TS_matrix, TS_IDs))
    
    points = np.vstack((basin_points, ts_points))

    area = 0.0

    for point in points:

        x_neighbours, y_neighbours, z_neighbours = PBC3D_neighbours(point)

        for neighbour in x_neighbours:
            if (Basin_matrix[tuple(neighbour)] not in B_IDs and TS_matrix[tuple(neighbour)] not in TS_IDs):
                area += bc

        for neighbour in y_neighbours:
            if (Basin_matrix[tuple(neighbour)] not in B_IDs and TS_matrix[tuple(neighbour)] not in TS_IDs):
                area += ac

        for neighbour in z_neighbours:
            if (Basin_matrix[tuple(neighbour)] not in B_IDs and TS_matrix[tuple(neighbour)] not in TS_IDs):
                area += ab

    return area

def PBC3D_neighbours(point):                                                                            #Subroutine that finds PBC neighbours of a given point
    i, j, k = point
    Nx, Ny, Nz = data.grid["shape"]
    x_neighbours = [((i + 1) % Nx, j, k), ((i - 1) % Nx, j, k)]
    y_neighbours = [(i, (j + 1) % Ny, k), (i, (j - 1) % Ny, k)]
    z_neighbours = [(i, j, (k + 1) % Nz), (i, j, (k - 1) % Nz)]
    return x_neighbours, y_neighbours, z_neighbours