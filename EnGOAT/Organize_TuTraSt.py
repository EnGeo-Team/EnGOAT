import numpy as np
from collections import defaultdict
import os
import logging
import data                                                                                             #Python file containing all the data-containing matrices used in the simulation
import get_process_cross_vector                                                                         #Fortran file containing subroutine that finds the cross vector of a cluster1->TS->cluster2 transition process
logger = logging.getLogger(__name__)

def Organize_transition_states():

    cluster_pair_map=defaultdict(list)                                                                  #Create a dictionary "cluster_pair_map", which is organised as; cluster_pair: list of TS_points on the interface between the two clusters
    for TS_point in data.TS_point.TS_point_list:
        key=frozenset(TS_point.clusters)
        cluster_pair_map[key].append(TS_point)
    
    ID=0
    for clusters, cluster_pair_interface in cluster_pair_map.items():                                   #For each cluster pair, loop through the interface points and organize them into transition states
        for point in cluster_pair_interface:
            if point.organized==False:
                ID+=1                                                                                   #Initiate a new transition state and asign points to it
                TS_points=[point.coordinates]
                TS_energies=[]
                point.organized=True
                point_index=0
                while point_index<len(TS_points):
                    data.TS_matrix.IDs[TS_points[point_index]]=ID
                    TS_energies.append(data.Energy_matrix[TS_points[point_index]])
                    neighbours=PBC3D_diagonal_neighbours(TS_points[point_index])
                    for neighbour_point in cluster_pair_interface:
                        if neighbour_point.coordinates in neighbours and neighbour_point.coordinates not in TS_points:
                                TS_points.append(neighbour_point.coordinates)
                                neighbour_point.organized=True
                    point_index=point_index+1
                TS_points=np.array(TS_points)
                TS_energies=np.array(TS_energies)
                E_min=TS_energies.min()

                """Save the transition state between given clusters containing above-defined TS_points. cluster pair is saved in a tuple to be able to differentiate betwee forward and backward transitions"""
                data.Transition_state.Transition_state_list.append(data.Transition_state(ID, tuple(clusters), TS_points, TS_energies, E_min))
    return None

def get_isolated_processes(cluster_family):

    transition_state_list=[]                                                                        #Create a list of all transition states between clusters of the current cluster family
    for transition_state in data.Transition_state.Transition_state_list:
        if bool(set(transition_state.clusters)&cluster_family):                                     #Checks whether either of the clusters is in the cluster family (either due to "-1" elements for the same cluster transition states)
            transition_state_list.append(transition_state)
    
    process_list=[]                                                                                 #Save all possible transitions between clusters that can occur within the tunnel system
    for transition_state in transition_state_list:
        if -1 in transition_state.clusters:                                                         #Case for a single cluster tunnel (double check if it works)
            cluster1_ID=next(x for x in transition_state.clusters if x != -1)
            cluster1=data.Cluster.Cluster_list[cluster1_ID-1]
            start_point=np.asfortranarray(cluster1.center, dtype=np.int64)
            for point in data.TS_point.TS_point_list:
                if point.coordinates==tuple(transition_state.TS_points[0]):
                    process_cross_vector=np.array(point.cross_vector)
            dE=transition_state.E_min-cluster1.E_min

            processforward=data.Process(cluster1_ID, cluster1_ID, transition_state.ID, start_point, start_point, dE, process_cross_vector)
            processbackward=data.Process(cluster1_ID, cluster1_ID, transition_state.ID, start_point, start_point, dE, -process_cross_vector)
            
        else:                                                                                       #Case for a regular cluster1--TS->cluster2 transition
            cluster1_ID=transition_state.clusters[0]
            cluster2_ID=transition_state.clusters[1]
            cluster1=data.Cluster.Cluster_list[cluster1_ID-1]
            cluster2=data.Cluster.Cluster_list[cluster2_ID-1]
            start_point=np.asfortranarray(cluster1.center, dtype=np.int64)
            end_point=np.asfortranarray(cluster2.center, dtype=np.int64)
            C1_C2_TS=np.asfortranarray((cluster1.ID, cluster2.ID, transition_state.ID), dtype=np.int64)

            if cluster1.boundary or cluster2.boundary:
                temp_cluster_matrix = np.array(data.Cluster_matrix.IDs, copy=True)
                for transition_state_2 in transition_state_list:
                    if C1_C2_TS[0] in transition_state_2.clusters and C1_C2_TS[1] in transition_state_2.clusters and C1_C2_TS[2] != transition_state_2.ID:
                        mask = data.TS_matrix.IDs == transition_state_2.ID
                        temp_cluster_matrix[mask] = 0
                process_cross_vector=get_process_cross_vector.get_process_cross_vector(start_point, end_point, C1_C2_TS, temp_cluster_matrix)   #np.array([0,0,0])#
            else:
                process_cross_vector=np.zeros(3).astype(int)
            
            dE_1=transition_state.E_min-cluster1.E_min
            dE_2=transition_state.E_min-cluster2.E_min

            processforward=data.Process(C1_C2_TS[0], C1_C2_TS[1], C1_C2_TS[2], start_point, end_point, dE_1, process_cross_vector)
            #processbackward=data.Process(C1_C2_TS[1], C1_C2_TS[0], C1_C2_TS[2], end_point, start_point, dE_2, -process_cross_vector)
            #logger.debug(f"Transition from {processforward.start_cluster} to {processforward.end_cluster}, dE={processforward.dE}")
            #logger.debug(f"starting point -> end point: {processforward.start_point}->{processforward.end_point}. Crossing vector: {processforward.process_cross_vector}")
            if tuple(process_cross_vector)==(100, 100, 100):                                        #Safety check in case the fortran subroutine didn't find the crossing vector
                    logger.warning("Endpoint for this process was not found!")
        process_list.append(processforward)
        #process_list.append(processbackward)

    return process_list

    
def Organize_tunnel_systems(level_max, organize_isolated):

    if organize_isolated == True:
        os.makedirs("TuTraSt_data", exist_ok=True) 
        iso_file = open(os.path.join("TuTraSt_data", "isolated_processes.dat"), mode = "w")
        iso_file.write(f"#Start basin \t end basin \tTransition state \t Energy barrier \t Start point [a, b, c] \tEnd point [a, b, c] \t process vector (a, b, c) \t Basin family\n")

    ID=1
    for cluster_family in data.Cluster.Cluster_families:
        family_tunnel_list=[]                                                                           #Find all tunnels that were identified in clusters belonging to the current cluster family
        for tunnel in data.Tunnel.Tunnel_set:
            if tunnel.cluster in cluster_family:
                family_tunnel_list.append(tunnel)
        if len(family_tunnel_list)==0:                                                                  #Skip cluster families for which no tunnel was identified
            
            if organize_isolated == True:
                process_list = get_isolated_processes(cluster_family)
                for process in process_list:
                    a_start, b_start, c_start = process.start_point
                    a_end, b_end, c_end = process.end_point
                    a_cross, b_cross, c_cross = process.process_cross_vector
                    start_string=f"[{a_start}, {b_start}, {c_start}]"
                    end_string=f"[{a_end}, {b_end}, {c_end}]"
                    cross_string=f"({a_cross}, {b_cross}, {c_cross})"
                    iso_file.write(f"{process.start_cluster:>14} \t {process.end_cluster:>10} \t{process.transition_state:>15} \t {process.dE:>12.6f} \t {start_string:>20} {end_string:>20} {cross_string:>20} {cluster_family}\n")
            continue

        height=level_max                                                                                #The energy level of the lowest energy point within the tunnel system
        for cluster_ID in cluster_family:                                                               #Make an entry in the Tunnel matrix
            cluster_mask=(data.Cluster_matrix.IDs==cluster_ID)
            data.Tunnel_matrix.IDs[cluster_mask]=ID
            if data.Cluster.Cluster_list[cluster_ID-1].E_min<height:
                height=data.Cluster.Cluster_list[cluster_ID-1].E_min

        directions=set()                                                                                #Find the directions of the tunnel system
        for tunnel in family_tunnel_list:
            directions.add(tunnel.direction)
        directions=np.array(list(directions))
        direction_breakthroughs=np.zeros(len(directions))                                               #Find the energy level at which the breakthrough in a given direction first occurs
        direction_breakthroughs[:]=level_max
        for idx, direction in enumerate(directions):
            for tunnel in family_tunnel_list:
                if tuple(tunnel.direction)==tuple(direction):
                    if direction_breakthroughs[idx]>tunnel.level:
                        direction_breakthroughs[idx]=tunnel.level

        transition_state_list=[]                                                                        #Create a list of all transition states between clusters of the current cluster family
        for transition_state in data.Transition_state.Transition_state_list:
            if bool(set(transition_state.clusters)&cluster_family):                                     #Checks whether either of the clusters is in the cluster family (either due to "-1" elements for the same cluster transition states)
                transition_state_list.append(transition_state)
        
        process_list=[]                                                                                 #Save all possible transitions between clusters that can occur within the tunnel system
        for transition_state in transition_state_list:
            if -1 in transition_state.clusters:                                                         #Case for a single cluster tunnel (double check if it works)
                cluster1_ID=next(x for x in transition_state.clusters if x != -1)
                cluster1=data.Cluster.Cluster_list[cluster1_ID-1]
                start_point=np.asfortranarray(cluster1.center, dtype=np.int64)
                for point in data.TS_point.TS_point_list:
                    if point.coordinates==tuple(transition_state.TS_points[0]):
                        process_cross_vector=np.array(point.cross_vector)
                dE=transition_state.E_min-cluster1.E_min

                processforward=data.Process(cluster1_ID, cluster1_ID, transition_state.ID, start_point, start_point, dE, process_cross_vector)
                processbackward=data.Process(cluster1_ID, cluster1_ID, transition_state.ID, start_point, start_point, dE, -process_cross_vector)
                
            else:                                                                                       #Case for a regular cluster1--TS->cluster2 transition
                cluster1_ID=transition_state.clusters[0]
                cluster2_ID=transition_state.clusters[1]
                cluster1=data.Cluster.Cluster_list[cluster1_ID-1]
                cluster2=data.Cluster.Cluster_list[cluster2_ID-1]
                start_point=np.asfortranarray(cluster1.center, dtype=np.int64)
                end_point=np.asfortranarray(cluster2.center, dtype=np.int64)
                C1_C2_TS=np.asfortranarray((cluster1.ID, cluster2.ID, transition_state.ID), dtype=np.int64)

                if cluster1.boundary or cluster2.boundary:
                    temp_cluster_matrix = np.array(data.Cluster_matrix.IDs, copy=True)
                    for transition_state_2 in transition_state_list:
                        if C1_C2_TS[0] in transition_state_2.clusters and C1_C2_TS[1] in transition_state_2.clusters and C1_C2_TS[2] != transition_state_2.ID:
                            mask = data.TS_matrix.IDs == transition_state_2.ID
                            temp_cluster_matrix[mask] = 0
                    process_cross_vector=get_process_cross_vector.get_process_cross_vector(start_point, end_point, C1_C2_TS, temp_cluster_matrix)   #np.array([0,0,0])#
                else:
                    process_cross_vector=np.zeros(3).astype(int)
                
                dE_1=transition_state.E_min-cluster1.E_min
                dE_2=transition_state.E_min-cluster2.E_min

                processforward=data.Process(C1_C2_TS[0], C1_C2_TS[1], C1_C2_TS[2], start_point, end_point, dE_1, process_cross_vector)
                processbackward=data.Process(C1_C2_TS[1], C1_C2_TS[0], C1_C2_TS[2], end_point, start_point, dE_2, -process_cross_vector)
                #logger.debug(f"Transition from {processforward.start_cluster} to {processforward.end_cluster}, dE={processforward.dE}")
                #logger.debug(f"starting point -> end point: {processforward.start_point}->{processforward.end_point}. Crossing vector: {processforward.process_cross_vector}")
                if tuple(process_cross_vector)==(100, 100, 100):                                        #Safety check in case the fortran subroutine didn't find the crossing vector
                        logger.warning("Endpoint for this process was not found!")
            process_list.append(processforward)
            process_list.append(processbackward)

        """Save the tunnel system containing the above-calculated parameters"""
        data.Tunnel_system.Tunnel_system_list.append(data.Tunnel_system(ID, cluster_family, height, directions, direction_breakthroughs, transition_state_list, process_list))
        ID+=1
    return None

def PBC3D_diagonal_neighbours(point):                                                                   #Subroutine that finds PBC neighbours of a given point
    i, j, k = point
    Nx, Ny, Nz = data.grid[0][0], data.grid[1][0], data.grid[2][0]
    neighbours = []
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            for dz in [-1, 0, 1]:
                if dx == 0 and dy == 0 and dz == 0:
                    continue  # Skip center point
                ni = (i + dx) % Nx
                nj = (j + dy) % Ny
                nk = (k + dz) % Nz
                neighbours.append((ni, nj, nk))
    return neighbours