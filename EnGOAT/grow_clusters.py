import numpy as np
import logging
import data                                                                                                                             #Python file containing all the data-containing matrices used in the simulation
import find_neighbours                                                                                                                  #Fortran subroutine that returns PBC neighbours of a given boundary point of a cluster
logger = logging.getLogger(__name__)

def grow_clusters(level, E_step):                                           

    Filling_matrix=np.zeros(data.Cluster.N_clusters[0])                                                                                 #Matrix keeping track of which clusters have been completely filled (1) or not (0) at a current level
    while sum(Filling_matrix)<data.Cluster.N_clusters[0]:                                           
        for current_cluster in data.Cluster.Cluster_list:                                           
            is_filled=1                                                                                                                 #A flag checking whether a cluster has been filled. If a new point is added in the neighbour search or two clusters are combined, the flag becomes 0
            if (Filling_matrix[current_cluster.ID-1]==0) and (current_cluster.active):                                                  #checks if the cluster has not yet been filled or removed

                mask=current_cluster.cluster_points[:, 6]==1                                                                            #Find boundary points of the current cluster
                boundary_points=np.where(mask)[0]       

                for current_point_index in boundary_points:         
                    current_point=np.asfortranarray(current_cluster.cluster_points[current_point_index], dtype=np.int64)                #Find neighbours of the boundary point and get their crossing information
                    neighbour_points=find_neighbours.find_neighbours(level, current_point, 
                                                                     data.cross_matrix.i, data.cross_matrix.j, data.cross_matrix.k)
                    boundary=0    
                    for neighbour in neighbour_points:
                        i_n=neighbour[0]
                        j_n=neighbour[1]
                        k_n=neighbour[2]
                        
                        if data.TS_matrix.Levels[current_point[0], current_point[1], current_point[2]]==0:
                            if data.Cluster_matrix.Levels[i_n, j_n, k_n]==0:
                                """If the neighbour has not yet been asigned to a cluster and is at a current level, asign it to the current cluster at the current level."""
                                if data.level_matrix[i_n, j_n, k_n]==level:                             
                                    is_filled=0                                                                         
                                    data.Cluster_matrix.Levels[i_n, j_n, k_n]=level
                                    data.Cluster_matrix.IDs[i_n, j_n, k_n]=current_cluster.ID
                                    data.cross_matrix.i[i_n, j_n, k_n]=neighbour[3]
                                    data.cross_matrix.j[i_n, j_n, k_n]=neighbour[4]
                                    data.cross_matrix.k[i_n, j_n, k_n]=neighbour[5]
                                    current_cluster.cluster_points=np.vstack([current_cluster.cluster_points, neighbour])
                                    current_cluster.point_energies=np.append(current_cluster.point_energies, data.Energy_matrix[i_n, j_n, k_n])

                                else:
                                    boundary=1
                            else:
                                """If the neighbour is already asigned to a cluster, make the following decision:"""
                                if data.Cluster_matrix.IDs[i_n, j_n, k_n]==current_cluster.ID:
                                    """If the neighbour was assigned to the current cluster, then check the crossing information: 
                                        If the boundry was not crossed, do nothing. 
                                        If it was, there is a tunnel"""
                                    i_crossing=neighbour[3]-data.cross_matrix.i[i_n, j_n, k_n] 
                                    j_crossing=neighbour[4]-data.cross_matrix.j[i_n, j_n, k_n]
                                    k_crossing=neighbour[5]-data.cross_matrix.k[i_n, j_n, k_n]
                                    if (i_crossing, j_crossing, k_crossing)==(0, 0, 0):
                                        pass
                                    else:
                                        tunnel_direction = (abs(i_crossing), abs(j_crossing), abs(k_crossing)) 
                                        data.Tunnel.Tunnel_set.add(data.Tunnel(tunnel_direction, current_cluster.ID, level))                                                                        #Add the tunnel to the tunnel list
                                        data.Tunnel.total_breakthrough_dimension=data.Tunnel.total_breakthrough_dimension+np.array(tunnel_direction)
                                        i_current=current_cluster.cluster_points[current_point_index][0]
                                        j_current=current_cluster.cluster_points[current_point_index][1]
                                        k_current=current_cluster.cluster_points[current_point_index][2]
                                        if data.TS_matrix.Levels[i_n, j_n, k_n]==0 and data.TS_matrix.Levels[i_current, j_current, k_current]==0:                                                   #Check whether the energy barrier between the crossing point and the cluster energy minimum exceeds the threshold value
                                            if data.Energy_matrix[i_n, j_n, k_n]-current_cluster.E_min>E_step or data.Energy_matrix[i_current, j_current, k_current]-current_cluster.E_min>E_step:  #If yes, add either the current or the neighbour point to the TS list
                                                if data.Energy_matrix[i_n, j_n, k_n]>=data.Energy_matrix[i_current, j_current, k_current]:          
                                                    data.TS_matrix.Levels[i_n, j_n, k_n]=data.level_matrix[i_n, j_n, k_n]
                                                    data.TS_point.TS_point_list.append(data.TS_point((i_n, j_n, k_n), current_cluster.ID, -1, (i_crossing, j_crossing, k_crossing)))                #The "second cluster" information to form the TS is in this case -1 in order for the TS not to be deleted in case some clusters combine with the current cluster later. Crossing information is also stored in the TS point for single cluster tunnels.
                                                    for idx in range(len(current_cluster.cluster_points)):
                                                        if (i_n, j_n, k_n)==(current_cluster.cluster_points[idx][0], current_cluster.cluster_points[idx][1], current_cluster.cluster_points[idx][2]):
                                                            current_cluster.cluster_points[idx][8]=1                                              
                                                else:
                                                    data.TS_matrix.Levels[i_current, j_current, k_current]=data.level_matrix[i_current, j_current, k_current]
                                                    data.TS_point.TS_point_list.append(data.TS_point((i_current, j_current, k_current), current_cluster.ID, -1, (i_crossing, j_crossing, k_crossing)))
                                                    current_cluster.cluster_points[current_point_index][8]=1 
                                            else:                                                                                                                                                   #If not, there is infinite diffusion in the given direction (output an error message)
                                                logger.warning(f"Zero energy barrier tunnel found! Infinite diffusion in {tunnel_direction} direction!")

                                else:
                                    """If the neighbour is asigned to another cluster, check the energy barrier dE for a transition from one cluster to the other:
                                        If dE is smaller than a threshold value, COMBINE the two clusters (=the two clusters become one)
                                        If not, the two clusters MERGE (=they retain their identity, but are added to a common family of clusters. They are separated by a transition state point)"""

                                    i_crossing=neighbour[3]-data.cross_matrix.i[i_n, j_n, k_n]
                                    j_crossing=neighbour[4]-data.cross_matrix.j[i_n, j_n, k_n]
                                    k_crossing=neighbour[5]-data.cross_matrix.k[i_n, j_n, k_n]
                                    neighbour_cluster=data.Cluster.Cluster_list[data.Cluster_matrix.IDs[neighbour[0], neighbour[1], neighbour[2]]-1]

                                    dE1=data.Energy_matrix[i_n, j_n, k_n]-max(current_cluster.E_min, neighbour_cluster.E_min)                                                                       #Minimum energy barrier of transition from one cluster to the other via the examined neighbour point
                                    dE2=data.Energy_matrix[current_point[0], current_point[1], current_point[2]]-max(current_cluster.E_min, neighbour_cluster.E_min) 
                                    dE=max(dE1, dE2)
                                    if dE<E_step:
                                        #logger.debug(f"Combining clusters {current_cluster.ID} and {neighbour_cluster.ID}. Energy difference = {dE}")
                                        is_filled=0
                                        for cluster_family in data.Cluster.Cluster_families:                                                                                                        #Identify the cluster families the current and the neighbour cluster belong to
                                            if current_cluster.ID in cluster_family:
                                                current_cluster_family=cluster_family
                                            if neighbour_cluster.ID in cluster_family:
                                                neighbour_cluster_family=cluster_family

                                        if current_cluster_family == neighbour_cluster_family:
                                            """If the two combining clusters belong to the same cluster family, check crossing information for possible tunnels"""
                                            if (i_crossing, j_crossing, k_crossing)==(0, 0, 0):
                                                pass   
                                            else: 
                                                tunnel_direction = (abs(i_crossing), abs(j_crossing), abs(k_crossing)) 
                                                data.Tunnel.Tunnel_set.add(data.Tunnel(tunnel_direction, current_cluster.ID, level))                                                                #Add the tunnel to the tunnel list 
                                                data.Tunnel.total_breakthrough_dimension=data.Tunnel.total_breakthrough_dimension+np.array(tunnel_direction)

                                            for neighbour_cluster_point in neighbour_cluster.cluster_points:                                                                                        #Add the points of the combining cluster to the current clusters point list
                                                
                                                data.cross_matrix.i[neighbour_cluster_point[0], neighbour_cluster_point[1], neighbour_cluster_point[2]]+=i_crossing
                                                data.cross_matrix.j[neighbour_cluster_point[0], neighbour_cluster_point[1], neighbour_cluster_point[2]]+=j_crossing
                                                data.cross_matrix.k[neighbour_cluster_point[0], neighbour_cluster_point[1], neighbour_cluster_point[2]]+=k_crossing
                                                neighbour_cluster_point[3]+=i_crossing
                                                neighbour_cluster_point[4]+=j_crossing
                                                neighbour_cluster_point[5]+=k_crossing
                                                
                                                current_cluster.cluster_points=np.vstack([current_cluster.cluster_points, neighbour_cluster_point.reshape(1, -1)])
                                                data.Cluster_matrix.IDs[neighbour_cluster_point[0], neighbour_cluster_point[1], neighbour_cluster_point[2]]=current_cluster.ID                      #Fix the minID matrix inputs
                                            for neighbour_cluster_point_energy in neighbour_cluster.point_energies:
                                                current_cluster.point_energies=np.append(current_cluster.point_energies, neighbour_cluster_point_energy)
                                            neighbour_cluster.active=False
                                            if current_cluster.E_min>neighbour_cluster.E_min:
                                                current_cluster.E_min=neighbour_cluster.E_min
                                                current_cluster.center=neighbour_cluster.center                                                                                                     #Update the minimum energy value of the combined cluster
                                            current_cluster_family.discard(neighbour_cluster.ID)

                                            for tunnel in data.Tunnel.Tunnel_set:
                                                if tunnel.cluster == neighbour_cluster.ID:
                                                    tunnel.cluster = current_cluster.ID                                                                                                             #Update the cluster family list of included clusters (discard the cluster that has been combined with the current cluster)

                                        else:
                                            """If the two combining clusters belong to different cluster families, it is impossible for a tunnel to exist. 
                                            Check the crossing information and if there is need, update the crossing information of the combining cluster's cluster family to reflect as if it was grown from the starting point of the current cluster"""
                                            if (i_crossing, j_crossing, k_crossing)==(0, 0, 0):
                                                pass
                                            else:                                                                                                                                                   #the crossing info of points from the combining cluster's cluster family has to be updated to reflect as if the cluster has been grown from the current cluster origin
                                                for combining_cluster in neighbour_cluster_family:
                                                    for combining_cluster_point in data.Cluster.Cluster_list[combining_cluster-1].cluster_points:
                                                        data.cross_matrix.i[combining_cluster_point[0], combining_cluster_point[1], combining_cluster_point[2]]+=i_crossing
                                                        data.cross_matrix.j[combining_cluster_point[0], combining_cluster_point[1], combining_cluster_point[2]]+=j_crossing
                                                        data.cross_matrix.k[combining_cluster_point[0], combining_cluster_point[1], combining_cluster_point[2]]+=k_crossing
                                                        combining_cluster_point[3]+=i_crossing
                                                        combining_cluster_point[4]+=j_crossing
                                                        combining_cluster_point[5]+=k_crossing

                                            for neighbour_cluster_point in neighbour_cluster.cluster_points:                                                                                        #Add the points of the combining cluster to the current clusters point list
                                                current_cluster.cluster_points=np.vstack([current_cluster.cluster_points, neighbour_cluster_point.reshape(1, -1)])
                                                data.Cluster_matrix.IDs[neighbour_cluster_point[0], neighbour_cluster_point[1], neighbour_cluster_point[2]]=current_cluster.ID                      #Fix the minID matrix inputs
                                            for neighbour_cluster_point_energy in neighbour_cluster.point_energies:
                                                current_cluster.point_energies=np.append(current_cluster.point_energies, neighbour_cluster_point_energy)
                                            neighbour_cluster.active=False                                                                                                                          #Deactivate the combining cluster
                                            if current_cluster.E_min>neighbour_cluster.E_min:
                                                current_cluster.E_min=neighbour_cluster.E_min
                                                current_cluster.center=neighbour_cluster.center                                                                                                     #Update the minimum energy value and the center point of the combined cluster

                                            merged_cluster_family=current_cluster_family|neighbour_cluster_family                                                                                   #Merge the cluster families
                                            data.Cluster.Cluster_families.remove(current_cluster_family)
                                            data.Cluster.Cluster_families.remove(neighbour_cluster_family)
                                            merged_cluster_family.discard(neighbour_cluster.ID)
                                            data.Cluster.Cluster_families.append(merged_cluster_family)

                                        TS_remove_list=[]
                                        for TS in data.TS_point.TS_point_list:                                                                                                                      #Update the transition states
                                            if neighbour_cluster.ID in TS.clusters:                                                                                                                 #Update the cluster indices (replace the ID of combining cluster in TSs with the current cluster ID)
                                                TS.clusters.remove(neighbour_cluster.ID)
                                                TS.clusters.add(current_cluster.ID)
                                            if len(TS.clusters)==1:                                                                                                                                 #If there is a transition state that includes only the one ccombining cluster, delete it
                                                data.TS_matrix.Levels[TS.coordinates[0], TS.coordinates[1], TS.coordinates[2]]=0
                                                for cluster_point in current_cluster.cluster_points:
                                                    if TS.coordinates==(cluster_point[0], cluster_point[1], cluster_point[2]):
                                                        cluster_point[8]=0
                                                TS_remove_list.append(TS)
                                        for TS in TS_remove_list:
                                            data.TS_point.TS_point_list.remove(TS)
                                        
                                        for tunnel in data.Tunnel.Tunnel_set:
                                            if tunnel.cluster == neighbour_cluster.ID:
                                                tunnel.cluster = current_cluster.ID  

                                    else:
                                        """if dE>E_step, the two clusters MERGE"""
                                        #logger.debug(f"Merging clusters {current_cluster.ID} and {neighbour_cluster.ID}. Energy difference = {dE}")  
                                        boundary=1
                                        for cluster_family in data.Cluster.Cluster_families:                                                                                                        #Identify the cluster families the current and the neighbour cluster belong to
                                            if current_cluster.ID in cluster_family:
                                                current_cluster_family=cluster_family
                                            if neighbour_cluster.ID in cluster_family:
                                                neighbour_cluster_family=cluster_family

                                        if current_cluster_family == neighbour_cluster_family:
                                            """If the two merging clusters belong to the same cluster family, check crossing information for possible tunnels"""                 
                                            if (i_crossing, j_crossing, k_crossing)==(0, 0, 0):
                                                pass   
                                            else: 
                                                tunnel_direction = (abs(i_crossing), abs(j_crossing), abs(k_crossing))  
                                                data.Tunnel.Tunnel_set.add(data.Tunnel(tunnel_direction, current_cluster.ID, level))                                                                #Add the tunnel to the tunnel list 
                                                data.Tunnel.total_breakthrough_dimension=data.Tunnel.total_breakthrough_dimension+np.array(tunnel_direction)

                                            """Identify the transition state between the two merging clusters (the two candidates are the current and the neighbouring point)"""
                                            i_current=current_cluster.cluster_points[current_point_index][0]
                                            j_current=current_cluster.cluster_points[current_point_index][1]
                                            k_current=current_cluster.cluster_points[current_point_index][2]
                                            if data.TS_matrix.Levels[i_n, j_n, k_n]==0 and data.TS_matrix.Levels[i_current, j_current, k_current]==0:                                               #If none of the two candidates is already the transition state
                                                if data.Energy_matrix[i_n, j_n, k_n]>=data.Energy_matrix[i_current, j_current, k_current]:                                                          #if neighbour is higher in E than current point, neighbour is TS
                                                    data.TS_matrix.Levels[i_n, j_n, k_n]=data.level_matrix[i_n, j_n, k_n]
                                                    data.TS_point.TS_point_list.append(data.TS_point((i_n, j_n, k_n), current_cluster.ID, neighbour_cluster.ID))
                                                    for idx in range(len(neighbour_cluster.cluster_points)):
                                                        if (i_n, j_n, k_n)==(neighbour_cluster.cluster_points[idx][0], neighbour_cluster.cluster_points[idx][1], neighbour_cluster.cluster_points[idx][2]):
                                                            neighbour_cluster.cluster_points[idx][8]=1                                                                                              #Set the TS info of neighbouring point in the list of cluster points of the neighbouring cluster to 1
                                                else:                                                                                                                                               #if current point is higher in E than neighbour, current point is TS
                                                    data.TS_matrix.Levels[i_current, j_current, k_current]=data.level_matrix[i_current, j_current, k_current]
                                                    data.TS_point.TS_point_list.append(data.TS_point((i_current, j_current, k_current), current_cluster.ID, neighbour_cluster.ID))
                                                    current_cluster.cluster_points[current_point_index][8]=1 

                                        else:
                                            """If the two merging clusters belong to different cluster families, it is impossible for a tunnel to exist. 
                                            Check the crossing information and if there is need, update the crossing information of the merging cluster's cluster family to reflect as if it was grown from the starting point of the current cluster"""                                                         
                                            merged_cluster_family=current_cluster_family|neighbour_cluster_family

                                            data.Cluster.Cluster_families.remove(current_cluster_family)
                                            data.Cluster.Cluster_families.remove(neighbour_cluster_family)
                                            data.Cluster.Cluster_families.append(merged_cluster_family)

                                            """Identify the transition state between the two merging clusters (the two candidates are the current and the neighbouring point)"""
                                            i_current=current_cluster.cluster_points[current_point_index][0]
                                            j_current=current_cluster.cluster_points[current_point_index][1]
                                            k_current=current_cluster.cluster_points[current_point_index][2]
                                            if data.TS_matrix.Levels[i_n, j_n, k_n]==0 and data.TS_matrix.Levels[i_current, j_current, k_current]==0:                                              #If none of the two candidates is already the transition state
                                                if data.Energy_matrix[i_n, j_n, k_n]>=data.Energy_matrix[i_current, j_current, k_current]:                                                         #if neighbour is higher in E than current point, neighbour is TS
                                                    data.TS_matrix.Levels[i_n, j_n, k_n]=data.level_matrix[i_n, j_n, k_n]
                                                    data.TS_point.TS_point_list.append(data.TS_point((i_n, j_n, k_n), current_cluster.ID, neighbour_cluster.ID))
                                                    for idx in range(len(neighbour_cluster.cluster_points)):
                                                        if (i_n, j_n, k_n)==(neighbour_cluster.cluster_points[idx][0], neighbour_cluster.cluster_points[idx][1], neighbour_cluster.cluster_points[idx][2]):
                                                            neighbour_cluster.cluster_points[idx][8]=1                                                                                             #Set the TS info of neighbouring point in the list of cluster points of the neighbouring cluster to 1
                                                else:                                                                                                                                              #if current point is higher in E than neighbour, current point is TS
                                                    data.TS_matrix.Levels[i_current, j_current, k_current]=data.level_matrix[i_current, j_current, k_current]
                                                    data.TS_point.TS_point_list.append(data.TS_point((i_current, j_current, k_current), current_cluster.ID, neighbour_cluster.ID))
                                                    current_cluster.cluster_points[current_point_index][8]=1    
                                       
                                            if (i_crossing, j_crossing, k_crossing)==(0, 0, 0):
                                                pass   
                                            else:                                                                                                                                                  #the crossing info of points from the merging cluster's cluster family has to be updated to reflect as if the cluster family has been grown from the current cluster origin
                                                for merging_cluster in neighbour_cluster_family:
                                                    for merging_cluster_point in data.Cluster.Cluster_list[merging_cluster-1].cluster_points:
                                                        data.cross_matrix.i[merging_cluster_point[0], merging_cluster_point[1], merging_cluster_point[2]]+=i_crossing
                                                        data.cross_matrix.j[merging_cluster_point[0], merging_cluster_point[1], merging_cluster_point[2]]+=j_crossing
                                                        data.cross_matrix.k[merging_cluster_point[0], merging_cluster_point[1], merging_cluster_point[2]]+=k_crossing                                                            
                                                        merging_cluster_point[3]+=i_crossing
                                                        merging_cluster_point[4]+=j_crossing
                                                        merging_cluster_point[5]+=k_crossing
                        else:
                            boundary=1    
                            
                    current_cluster.cluster_points[current_point_index][6]=boundary                                                     #Update the boundary information of the current point
            Filling_matrix[current_cluster.ID-1]=is_filled                                                                              #Update the filling matrix of the current cluster
    
    data.Tunnel.total_breakthrough_dimension = (data.Tunnel.total_breakthrough_dimension != 0).astype(int)                              #Save the total dimension of the tunnel systems discovered up until this level 
    return None