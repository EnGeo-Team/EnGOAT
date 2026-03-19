import numpy as np
import data
import logging
import initiate_cluster                                                                                                         #Fortran file containing the subroutine for finding new clusters at a current level
logger = logging.getLogger(__name__)

def find_new_clusters(level):

    unexplored_x, unexplored_y, unexplored_z=np.where(data.Cluster_matrix.Levels==0)                                            #Find all unexplored points
    
    for unexplored_point in range(len(unexplored_x)):       
        i=unexplored_x[unexplored_point]        
        j=unexplored_y[unexplored_point]        
        k=unexplored_z[unexplored_point]        
        if data.level_matrix[i,j,k]==level:                                                                                     #If the point is at the current level, it is a new cluster
            N_clusters_old=np.array([data.Cluster.N_clusters[0]], dtype=np.int64)               
            E_min=initiate_cluster.find_new_clusters(level, data.Cluster.N_clusters, i, j, k,       
                                               data.level_matrix, data.Cluster_matrix.Levels, data.Cluster_matrix.IDs,         
                                               data.cross_matrix.i, data.cross_matrix.j, data.cross_matrix.k,       
                                               data.Energy_matrix)                                                              #Call Fortran subroutine to build new clusters
            if N_clusters_old[0]!=data.Cluster.N_clusters[0]:                                                                   #Check whether a new cluster has been found
                Cluster_points=np.loadtxt("temp1.dat", dtype=np.int64)                      
                point_energies=np.loadtxt("temp2.dat", dtype=np.float64)                        
                if Cluster_points.ndim == 1:                        
                    Cluster_points=Cluster_points.reshape(1, -1)                                                                #if only one point in a cluster, reshape into a 2D matrix for compatibility
                    point_energies=point_energies.reshape(1, -1)                        
                Cluster_points[:,0]=Cluster_points[:,0]-1                                                                       #Fortran <-> python indexing
                Cluster_points[:,1]=Cluster_points[:,1]-1               
                Cluster_points[:,2]=Cluster_points[:,2]-1
                Emin_point=np.argmin(point_energies)
                center=(Cluster_points[Emin_point,0], Cluster_points[Emin_point,1], Cluster_points[Emin_point,2])
                data.Cluster.Cluster_list.append(data.Cluster(data.Cluster.N_clusters[0], Cluster_points, E_min, center, point_energies))
                data.Cluster.Cluster_families.append({data.Cluster.N_clusters[0]})
                #logger.debug(f"New cluster was found. Cluster ID: {int(data.Cluster.N_clusters[0])}")
    
    return None