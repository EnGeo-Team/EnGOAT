import numpy as np
import logging
import Read_cube_data                                                   #Python file containing subroutine to extract data from the "".cube" file
logger = logging.getLogger(__name__)

def create_matrices(E_step, E_cutoff, E_unit, cube_file):

    global Energy_matrix, level_matrix                                  #3D matrix that holds 1) energies (float numbers) and 2) energy levels (integers)
    global atoms, grid, grid_size                                       #1) number of atoms and origin from the .cube file, 2) grid information from the .cube file, and 2) a tuple containing grid size in angstrom in the unit cell vector direcitons
    global Cluster_matrix, TS_matrix, Tunnel_matrix                     #1) ID_matrix with information about cluster IDs and levels the points belong to, 2) ID_matrix with information about transition state IDs and levels the points belong to, and 3) ID_matrix with information about tunnel system IDs and levels the points belong to. Also serves as accessible volume mask
    global cross_matrix                                                 #Object containing matrices holding crossing information

    logger.info("Reading cube data...")                                 #Read .cube file data
    atoms, grid, grid_size, Energy_matrix = Read_cube_data.Read_cube_data(E_unit, cube_file) 
    logger.info("Cube data read successfully.")

    """Make a level matrix, as well as the minID and cross matrices"""
    Energy_matrix[Energy_matrix>E_cutoff]=E_cutoff
    Energy_matrix[np.isinf(Energy_matrix)]=E_cutoff
    level_matrix=np.ceil(Energy_matrix/E_step+10**(-8))                 #+10**(-8) Makes energy levels start from 1 to avoid 0 issues

    grid_shape=Energy_matrix.shape                                      #Number of points in x, y, and z directions
    Cluster_matrix = ID_matrix(grid_shape)
    TS_matrix = ID_matrix(grid_shape)
    Tunnel_matrix = ID_matrix(grid_shape)
    cross_matrix = CROSS_matrix(grid_shape)

    """Make matrices fortran compatible"""
    Energy_matrix=np.asfortranarray(Energy_matrix, dtype=np.float64)
    level_matrix = np.asfortranarray(level_matrix, dtype=np.int64)   
    Cluster_matrix.Levels = np.asfortranarray(Cluster_matrix.Levels, dtype=np.int64)
    Cluster_matrix.IDs = np.asfortranarray(Cluster_matrix.IDs, dtype=np.int64)
    TS_matrix.Levels = np.asfortranarray(TS_matrix.Levels, dtype=np.int64)
    TS_matrix.IDs = np.asfortranarray(TS_matrix.IDs, dtype=np.int64)    
    Tunnel_matrix.Levels = np.asfortranarray(Tunnel_matrix.Levels, dtype=np.int64)
    Tunnel_matrix.IDs = np.asfortranarray(Tunnel_matrix.IDs, dtype=np.int64)    
    cross_matrix.i = np.asfortranarray(cross_matrix.i, dtype=np.int64)
    cross_matrix.j = np.asfortranarray(cross_matrix.j, dtype=np.int64)
    cross_matrix.k = np.asfortranarray(cross_matrix.k, dtype=np.int64)

class ID_matrix:                                                        #Matrix to keep track of explored points (stores information of which cluster/TS a point belongs to and at which level it was discovered)
    def __init__(self, grid_shape):     
        self.Levels=np.zeros(grid_shape)                                #Energy level at which each grid point was first explored
        self.IDs=np.zeros(grid_shape)                                   #Cluster ID of each grid point. Later (in the organize_transition_states subroutine), the values of minIDmatrix.Clusters are set to -1 for TS points! 

class CROSS_matrix:                                                     #Matrix to keep track of crossing information     
    def __init__(self, grid_shape):                     
        self.i=np.zeros(grid_shape)                                     #Crossing in x direction (values: -1, 0, 1)
        self.j=np.zeros(grid_shape)                                     #Crossing in y direction (values: -1, 0, 1)
        self.k=np.zeros(grid_shape)                                     #Crossing in z direction (values: -1, 0, 1)

class Cluster:                                                          #Each cluster formed is an object with the following attributes:
    def __init__(self, ID, cluster_points, E_min, center, point_energies):                      
        self.ID=ID                                                      #Identity number
        self.cluster_points=cluster_points                              #list of points belonging to the current cluster. Each point is described as an array [x, y, z, Vx, Vy, Vz, B, level, TS], where xyz are point coordinates, Vxyz is the crossing vector, B is boundary info (0 or 1), and TS is a flag indicating whether a point is a TS (1) or not (0)
        self.E_min=E_min                                                #Energy of the point with the lowest energy value in the cluster
        self.active=True                                                #Flag to indicate whether a cluster is still active or has been combined with another and removed
        self.point_energies=point_energies                              #List of energies of the corresponding points in the cluster (is filled after all points are identified)
        self.center=center                                              #The coordinates of the minimum energy point saved in a tuple (i, j, k)
        self.boundary=False                                             #A flag indicating whether the cluster is situated at the boundary of the unit cell
        self.Boltzmann_weighted_V_fraction=None                         #Boltzmann weighted volume fraction of a cluster at a given temperature. Used for constructing a graph network output.
    Cluster_list=[]                                                     #List of all identified clusters
    Cluster_families=[]                                                 #List of sets, containing cluster IDs of merged cluster groups -> cluster families. Total power of all sets is always equal to the number of ACTIVE clusters
    N_clusters=np.asfortranarray([0], dtype=np.int64)                   #Number of clusters found on this level. A 1x1 array to ensure compatibility with fortran

class TS_point:                                                         #Each Transition state point is an object with the following attributes:
    def __init__(self, coordinates, cluster1, cluster2, cross_vector=None):        
        self.coordinates=coordinates                                    #Transition state point coordinates, saved in a tuple (i, j, k)
        self.clusters={cluster1, cluster2}                              #Cluster IDs merged through a given transition state point, saved in a set
        self.organized=False                                            #A flag used in the process of organizing TS points into transition states
        self.cross_vector=cross_vector                                  #For single cluster tunnels (transition states within the same cluster), crossing information is recorded here.
    TS_point_list=[]                                                    #List of all identified transition state points
    def __eq__(self, other):                                            #Ensure that equality checking works correctly
        if not isinstance(other, TS_point):
            return False
        return (self.coordinates == other.coordinates)

class Transition_state:                                                 #Individual transition state points are organized together in transition states, whch have the following attributes:
    def __init__(self, ID, clusters, TS_points, point_energies, E_min):
        self.ID=ID                                                      #Identity number
        self.clusters=clusters                                          #Cluster IDs merged through a given transition state, saved in a tuple (cluster1, cluster2)
        self.TS_points=TS_points                                        #List of all points belonging to a given transition state plane between the two clusters (saved as tuples (i, j, k))
        self.point_energies=point_energies                              #List of all energies of the corresponding points in the TS
        self.E_min=E_min                                                #Energy value of the lowest energy point on the transition state surface
        self.Process_cross_vector = None                                #PBC crossing information
    Transition_state_list=[]                                            #List of all identified transition states

class Tunnel:                                                           #Unique breakthroughs recorded are saved as objects with the following attributes:
    def __init__(self, direction, cluster, level):     
        self.direction=direction                                        #Direction of the tunnel
        self.cluster=cluster                                            #Cluster, for which this tunnel was identified. Later, the cluster family can be restored form that information
        self.level=level                                                #Level at which the tunnel was identified
    Tunnel_set=set()                                                    #Set containing all identified tunnels
    total_breakthrough_dimension=np.zeros(3)                            #np.array containing the total dimension of all breakthroughs recorded so far (1 if any breakthroughs were recorded along a given axis and 0 if not)

    def __eq__(self, other):                                            #Ensure that duplicate checking for Tunnel_set works correctly
        if not isinstance(other, Tunnel):
            return False
        return (self.direction == other.direction and
                self.cluster == other.cluster and
                self.level == other.level)
    def __hash__(self):
        return hash((self.direction, self.cluster, self.level))

class Tunnel_system:                                                    #All identified breakthroughs are organized into their corresponding Tunnel systems -> merged cluster families
    def __init__(self, ID, cluster_family, height, directions, direction_breakthroughs, transition_state_list, process_list):
        self.ID=ID                                                      #Identity number
        self.cluster_family=cluster_family                              #Cluster family forming a tunnel
        self.height=height                                              #Energy value of the lowest energy point within the tunnel system
        self.directions=directions                                      #Directions of the tunnel system (np.array)
        self.direction_breakthroughs=direction_breakthroughs            #Levels at which a breakthrough in a given direction first occured (np.array)
        self.dimension=np.linalg.matrix_rank(directions)                #Dimensionality of the tunnel systems = number of linearly independent directions
        self.transition_state_list=transition_state_list                #list of all transition states (objects, see above) between clusters of the tunnel system
        self.process_list=process_list                                  #List of all possible processes (objects, see below) occuring in the tunnel system
        self.surface_area=None                                          #Surface area at cutoff energy of a tunnel
        self.volume=None                                                #Volume under the cutoff energy of the tunnel
        self.V_fraction=None                                            #Volume fraction of the tunnel system within the unit cell
        self.Boltzmann_weighted_V_fraction=None                         #Boltzmann weighted volume fraction of the tunnel at a given temperature. Used as a weight for the contribution of a given tunnel to the overall diffusion in kMC
        self.histogram=None                                             #Histogram of points in the tunnel system at a given energy level
    Tunnel_system_list=[]                                               #List of all tunnel systems identified in the structure

class Process:                                                          #Each transition C1--TS->C2 is stored as an object
    def __init__(self, start_cluster, end_cluster, transition_state, start_point, end_point, dE, process_cross_vector):
        self.start_cluster=start_cluster                                #Starting cluster
        self.end_cluster=end_cluster                                    #Destination cluster
        self.transition_state=transition_state                          #Given transition state between the two clusters
        self.start_point=start_point                                    #Center point of the starting cluster
        self.end_point=end_point                                        #Center point of the destination cluster
        self.dE=dE                                                      #Energy difference between the TS minimum energy point and starting cluster minimum energy point
        self.process_cross_vector=process_cross_vector                  #Crossing information of the C1--TS->C2 transition
        self.k=None                                                     #Rate constant of the C1--TS->C2 transition
        self.distance=None                                              #PBC distance between C1 and C2, expressed as a vector [a, b, c], where a, b, and c have the directions of the unit cell vectors specified in the cube file