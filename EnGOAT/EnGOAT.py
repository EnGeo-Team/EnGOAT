import numpy as np
from scipy import constants
import time
import os
import json
import logging
from datetime import datetime
import data                                                                                                                         #Python file containing all the data-containing matrices used in the simulation
import find_new_clusters                                                                                                            #Python file containing the subroutine for finding new clusters at current level
import grow_clusters                                                                                                                #Python file containing the subroutine for growing clusters at a current level
import Organize_TuTraSt                                                                                                             #Python file containing the subroutine for grouping individual transition state points into transition state planes between pairs of clusters
import get_topological_descriptors                                                                                                  #Python file containing the subroutine for calculating surface areas, volumes, boltzmann weighted volumes, ... of the whole unit cell and each tunnel system individually
import PBC_minimax                                                                                                                  #A python file containing a function to find a minimum energy path for a given offset in the graph network
import kMC                                                                                                                          #Python file containing the subroutine for performing kinetic Monte Carlo simulation on a tunnel system to obtain the diffusion coefficients

"""Initiate the log file and timer"""
console_handler = logging.StreamHandler()
file_handler = logging.FileHandler("output.log", mode="w")
console_handler.setLevel(logging.INFO)
file_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(message)s")
console_handler.setFormatter(formatter)
file_handler.setFormatter(formatter)
logger = logging.getLogger()
logger.setLevel(logging.INFO)
logger.addHandler(console_handler)
logger.addHandler(file_handler)

logger.info("*" * 80)                                                                                                               #Create header
logger.info("")                                                                                                     
logger.info(f"{'E  n  G  O  A  T':^80}") 
logger.info("") 
logger.info("*" * 80)   
logger.info(f"{'The TuTraSt Algorithm Based Energy GeOmetry Analysis Toolkit':^80}")  
logger.info("") 
logger.info(f"{datetime.now().strftime('%Y-%m-%d'):^80}")   
logger.info("*" * 80)   
logger.info("") 

#"""read user defined input"""    #Legacy input file                             
#input_data = np.loadtxt('input.param', comments='%')   
#E_unit=int(input_data[0])                                                                                                          #Energy unit:  1=kJ/mol, 2=kcal/mol, 3=Ry, 4=eV, 5=Hartee.
#N_Temp=int(input_data[1])                                                                                                          #Number of temperatures
#T_list=[]                                                                                                                          #List of temperatures
#for i in range(N_Temp):                        
#    T_list.append(input_data[2+i])                         
#run_kMC=int(input_data[2+N_Temp])                                                                                                  #Run kMC? 0=no, 1=yes
#plot_msd=int(input_data[3+N_Temp])                                                                                                 #Plot MSD? 0=no, 1=yes.
#N_kMC_steps=int(input_data[4+N_Temp])                                                                                              #Number of kMC steps
#print_every=int(input_data[5+N_Temp])                                                                                              #Print trajectory every n steps.                            #redundant
#N_kMC_runs=int(input_data[6+N_Temp])                                                                                               #Number of kMC simulations to run for averaging.
#N_particles=int(input_data[7+N_Temp])                                                                                              #Number of particles for kMC.                               #redundant
#per_tunnel=int(input_data[8+N_Temp])                                                                                               #%number of particles given: 0=in total, 1=per tunnel.      #redundant
#particle_mass=input_data[9+N_Temp]                                                                                                 #Mass of the diffusing particle.
#E_step=input_data[10+N_Temp]                                                                                                       #Energy step in kJ/mol
#E_cutoff=input_data[11+N_Temp]                                                                                                     #Energy cutoff on kJ/mol
#coarsening_factor=(int(input_data[12+N_Temp]), int(input_data[13+N_Temp]), int(input_data[14+N_Temp]))                             #Coarsening factor in x, y, z direction (not implemented)
#   
#cube_file='grid_UC.cube'                                                                                                           #Cube file from which the energy matrix is read

lines = []  
with open('input.param', 'r') as f:    
    for line in f:  
        line = line.split('%')[0].strip()   
        if line:    
            lines.append(line)  
E_unit          = int(lines[0])                                                                                                     #Energy unit:  1=kJ/mol, 2=kcal/mol, 3=Ry, 4=eV, 5=Hartee.
E_step          = float(lines[1])                                                                                                   #Energy step in kJ/mol
E_cutoff        = float(lines[2])                                                                                                   #Energy cutoff on kJ/mol
autocorrect     = int(lines[3])                                                                                                     #Dynamic Energy step correction? 0=no, 1=yes
run_kMC         = int(lines[4])                                                                                                     #Run kMC? 0=no, 1=yes
N_Temp          = int(lines[5])                                                                                                     #Number of temperatures
T_list = [float(lines[6 + i]) for i in range(N_Temp)]                                                                               #List of temperatures
N_kMC_steps     = int(lines[6 + N_Temp])                                                                                            #Number of kMC steps
N_kMC_runs      = int(lines[7 + N_Temp])                                                                                            #Number of kMC simulations to run for averaging.
particle_mass   = float(lines[8 + N_Temp])                                                                                          #Mass of the diffusing particle.
cube_file       = lines[9 + N_Temp]                                                                                                 #Cube file name from which the energy matrix is read

logger.info(">>> INITIALIZATION <<<")   
logger.info("-"*80) 
logger.info(f"{'Input file:':<25} {'input.param':>20}") 
logger.info(f"{'Cube file:':<25} {cube_file:>20}")  
logger.info("-"*80) 
logger.info("") 

timer = time.process_time()                                                                                                         #Timer for total time

logger.info(">>> READING USER DEFINED INPUT <<<")
timer_writing = time.process_time()                                                                                            
"""Read grid data and create matrices and classes (module 'data.py'):
Energy_matrix: 3D matrix that holds energy values (floats) of points on the grid
level_matrix: Energy_matrix, but the energy values of points are rounded up to the nearest energy level (integers)
grid, grid_size: hold information of the size of the grid
Cluster_matrix(object of the ID_matrix class): matrix that is filled during the course of the program, holds information of clusters (levels and IDs)
TS_matrix(object of the ID_matrix class): matrix that is filled during the course of the program, holds information of transition states (levels and IDs)
cross_matrix(object of the CROSS_matrix class): matrix that is filled during the course of the program, holds information about the boundary crossing for identifying tunnels
Cluster(class): defines cluster as an object and holds the list of all clusters
TS_point(class): defines transition state points as objects and holds the list of all transition state points
Transition_state(class): After all points with energy lower than E_cutoff have been checked, neighbouring TS_points between a pair of clusters are mergd into a transition state surface, stored as an object
Tunnel(class): defines tunnels as objects and holds the list of all tunnels
Tunnel_system(class): After all points with energy lower than E_cutoff have been checked, all identified breakthroughs are organized into Tunnel_systems (=merged cluster families), stored as an object
Process(class): Objcet containing information about a given transition (cluster 1 --TS-> cluster 2) from one cluster to the other through a given transition state
"""
data.create_matrices(E_step, E_cutoff, E_unit, cube_file)

logger.info("-" * 80)
logger.info(f"{'Number of atoms:':<25} {data.atoms[0]:>10}")
logger.info(f"{'Origin (x,y,z):':<25} ({data.atoms[1]:>8.3f}, {data.atoms[2]:>8.3f}, {data.atoms[3]:>8.3f})")
logger.info("-" * 80)
logger.info(f"{'Grid points (Na,Nb,Nc):':<25} ({data.grid[0][0]:>8}, {data.grid[1][0]:>8}, {data.grid[2][0]:>8})")
logger.info(f"{'a-axis vector:':<25} ({data.grid[0][1]:>8.3f}, {data.grid[0][2]:>8.3f}, {data.grid[0][3]:>8.3f})")
logger.info(f"{'b-axis vector:':<25} ({data.grid[1][1]:>8.3f}, {data.grid[1][2]:>8.3f}, {data.grid[1][3]:>8.3f})")
logger.info(f"{'c-axis vector:':<25} ({data.grid[2][1]:>8.3f}, {data.grid[2][2]:>8.3f}, {data.grid[2][3]:>8.3f})")
logger.info("-" * 80)
logger.info(f"Grid successfully loaded in {time.process_time() - timer_writing:.2f} s")
logger.info("")        

logger.info(">>> POTENTIAL ENERGY GRID ANALYSIS <<<")
timer_analysis = time.process_time()              

def TuTraSt_analysis(autocorrect, E_step):                                                                                          #Defined as a function to be able to be restarted for the E_step adjustment in the while loop below!
    N_levels=int(E_cutoff/E_step)                                                                                                   #Number of levels in the level matrix
    level_max=int(np.max(data.level_matrix))                                                                                        #Maximum level in the level matrix
    level_min=int(np.min(data.level_matrix))                                                                                        #Minimum level in the level matrix
    E_volume=np.zeros(level_max)                                                                                                    #Energy volume (number of grid points below each energy level)
    logger.info(f"{'Number of levels:':<25} {int(N_levels):>10}")
    logger.info(f"{'E step:':<25} {E_step:>10}")
    logger.info(f"{'E cutoff:':<25} {E_cutoff:>10}")
    logger.info("-" * 80)
    logger.info(f"{'Level':>6} {'Energy (kJ/mol)':>17} {'Volume (pts)':>15} {'Breakthrough':>15} {'Δt(s)':>10} {'Total(s)':>10}")
    logger.info("-" * 80)

    """Explore the potential energy grid one energy level at a time.
    First, grow all existing clusters at the current level: identify all cluster combinings and merges, transition states, and breakthroughs (tunnels)
    Then, find new clusters at the current level.
    """
    breakthrough_level = np.zeros(3, dtype=int)                                                                                     #Array that stores the levels at which the breakthrough in a, b, and c direction has first occured
    breakthrough_level[:]=level_max
    breakthrough="/"                                                                                                                #Unit cell directions in which the breakthrough occurs, just for display
    for level in range(level_min, level_max):                                                                                       #Loop over all levels except the highest one    level_max OR level_max-1, check
        E_volume[level]=np.sum(data.level_matrix<=level)                                                                            #Volume (number of grid points) at this energy level
        timer_level = time.process_time()                                                                                           #Timer for checking each level

        if level>level_min:
            grow_clusters.grow_clusters(level, E_step)                                                                              #Grow existing clusters at the current level
            
        find_new_clusters.find_new_clusters(level)                                                                                  #Find new clusters at the current level
        
        first_breakthrough = (breakthrough == "/")                                                                                  #Flag indicating if a first breakthrough in the system occurred at a given level (part 2 below)
        if sum(data.Tunnel.total_breakthrough_dimension)!=0:
            breakthrough=str(data.Tunnel.total_breakthrough_dimension)
            if data.Tunnel.total_breakthrough_dimension[0]==1 and level<breakthrough_level[0]:
                breakthrough_level[0]=level
            if data.Tunnel.total_breakthrough_dimension[1]==1 and level<breakthrough_level[1]:
                breakthrough_level[1]=level
            if data.Tunnel.total_breakthrough_dimension[2]==1 and level<breakthrough_level[2]:
                breakthrough_level[2]=level
        logger.info(
            f"{level:6d} "
            f"{E_step * level:17.2f} "
            f"{int(E_volume[level]):15d} "
            f"{breakthrough:>15} "
            f"{time.process_time()-timer_level:10.2f} "
            f"{time.process_time()-timer_analysis:10.2f}"
        )
        first_breakthrough = breakthrough != "/" and first_breakthrough

        restart=False
        if first_breakthrough and autocorrect:                                                                                      #If the breakthrough occurs between levels 10 and 20 proceed, otherwise restart the search with adjusted E_step (except if autocorrect flag is set to 0 (dont adjust), or 2 (stop at first breakthrough) by user)
            if autocorrect == 1:    
                if 10 <= level <= 20:
                    restart=False
                else:
                    pass
                    restart=True
                    if level<10:
                        E_step=E_step/(20//level)
                    if level>20:
                        E_step=E_step*(level//10)
                    break
            elif autocorrect == 2:
                restart=False
                break
    return restart, E_step, N_levels, level_max, level_min, breakthrough_level

restart=True
while restart:                                                                                                                      #E_step adjustment routine
    E_step_old=E_step
    restart, E_step, N_levels, level_max, level_min, breakthrough_level = TuTraSt_analysis(autocorrect, E_step)
    if restart:
        data.Cluster.Cluster_list=[]
        data.Cluster.N_clusters=np.asfortranarray([0], dtype=np.int64)
        data.Cluster.Cluster_families=[]
        data.TS_point.TS_point_list=[]
        data.Tunnel.Tunnel_set=set()
        data.Tunnel.total_breakthrough_dimension=np.zeros(3)
        logger.info("-"*80)
        logger.info(f"Adjusting energy step from {E_step_old} -> {E_step}")
        data.create_matrices(E_step, E_cutoff, E_unit, cube_file)
        logger.info("Restarting potential energy grid analysis:")
        logger.info("-"*80)

os.remove("temp1.dat")                                                                                                              #Remove the temporary files created by initiate_cluster.f90
os.remove("temp2.dat")
logger.info("-" * 80)
logger.info(f"Potential grid analyzed successfully in {time.process_time()-timer_analysis:.2f} s")
logger.info("")

logger.info(">>> TUTRAST ANALYSIS <<<")
timer_organization=time.process_time()

"""Remove all 1-point clusters (noise)"""                       
for cluster in data.Cluster.Cluster_list:                       
    if cluster.active:                      
        if len(cluster.cluster_points)==1:                      
            i_remove=cluster.cluster_points[0][0]                       
            j_remove=cluster.cluster_points[0][1]                       
            k_remove=cluster.cluster_points[0][2]                       
            data.Cluster_matrix.IDs[i_remove, j_remove, k_remove]=0                                                                 #Remove the minID matrix input of this cluster, set ID to 0

            for TS in data.TS_point.TS_point_list[:]:                                                                               #Remove all transition states of this cluster
                data.TS_matrix.Levels[TS.coordinates[0], TS.coordinates[1], TS.coordinates[2]] = 0                                  #Remove the TS_matrix entry of the transition state
                if cluster.ID in TS.clusters:                       
                    TS.clusters.remove(cluster.ID)                      
                    remaining_ID, = TS.clusters                     
                    for point in data.Cluster.Cluster_list[remaining_ID-1].cluster_points:                                          #Unset the point in the neighbouring cluster as TS
                        if (point[0], point[1], point[2]) == TS.coordinates:                        
                            point[8] = 0                        
                    data.TS_point.TS_point_list.remove(TS)                                                                          #Remove the transition state from TS_point_list                                             
            for cluster_family in data.Cluster.Cluster_families:                                                                    #Remove this cluster from all cluster families
                if cluster.ID in cluster_family:                        
                    cluster_family.remove(cluster.ID)                       
            cluster.active=False                                                                                                    #Inactivate the cluster

"""Identify clusters at the boundary of the unit cell"""        
border_mask = np.zeros(data.Cluster_matrix.IDs.shape, dtype=bool)        
border_mask[0, :, :] = True     
border_mask[-1, :, :] = True        
border_mask[:, 0, :] = True     
border_mask[:, -1, :] = True        
border_mask[:, :, 0] = True     
border_mask[:, :, -1] = True        
border_indices = np.argwhere(border_mask)                                                                                           # Get indices of all border points

boundary_clusters=set()
for i, j, k in border_indices:                      
    boundary_clusters.add(data.Cluster_matrix.IDs[i][j][k])                      
boundary_clusters.remove(0)                     
for ID in boundary_clusters:                        
    data.Cluster.Cluster_list[ID-1].boundary=True                                                                                   #Set the boundary flag of all clusters at the boundary to True

"""Organize transition state points into transition state surfaces, and identify tunnel systems"""                      
Organize_TuTraSt.Organize_transition_states()                                                                                       #A python subroutine that organizes transition state points into transition states, fills the TS_matrix.IDs, and sets the minIDmatrix.Clusters values of TS points to -1
Organize_TuTraSt.Organize_tunnel_systems(level_max)                                                                                 #A python subroutine that organizes recorded breakthroughs into cluster-family-wide tunnel systems          
Total_area, Accessible_area, Total_volume, Total_volume_fraction, Accessible_volume = get_topological_descriptors.get_topological_descriptors(int(N_levels))

logger.info(f"{'Total volume:':<25} {Total_volume:10.4f} Å³ {Total_volume_fraction:10.4f} %")
logger.info(f"{'Accessible volume:':<25} {Accessible_volume:10.4f} Å³ {Accessible_volume/Total_volume*Total_volume_fraction:10.4f} %")
logger.info(f"{'Total surface area:':<25} {Total_area:10.4f} Å² {Total_area/(Total_volume/Total_volume_fraction):10.4f} Å²/Å³")
logger.info(f"{'Accessible surface area:':<25} {Accessible_area:10.4f} Å² {Accessible_area/(Total_volume/Total_volume_fraction):10.4f} Å²/Å³")
logger.info("")
logger.info(f"{len(data.Tunnel_system.Tunnel_system_list)} tunnel systems were identified:")
logger.info("-" * 80)

os.makedirs("tunnel data", exist_ok=True)                                                                                           #Create a folder where tunnel systems will be saved as graph networks, where basins will be defined with their corresponding boltzmann weighted probabilities and a process list of possible transitions defined each with the corresponding rate
Ea_abc = np.zeros(3)                                                                                                                #minimum energy barrier along a given crystal lattice direction
Ea_abc[:] = E_cutoff
margin = E_step                                                                                                                     #E tolerance for finding the best energy path along a direction (the paths that differ for at most this energy will be evaluated based on length)
for tunnel_system in data.Tunnel_system.Tunnel_system_list:
    logger.info(f"Tunnel system #{tunnel_system.ID}")
    logger.info("-" * 80)
    logger.info(f"{'Volume:':<25} {tunnel_system.volume:10.4f} Å³ {tunnel_system.V_fraction:10.4f} %")
    logger.info(f"{'Surface area:':<25} {tunnel_system.surface_area:10.4f} Å² {tunnel_system.surface_area/(Total_volume/Total_volume_fraction):10.4f} Å²/Å³")
    logger.info(f"{'Lowest energy point:':<25} {tunnel_system.height:10.4f} kJ/mol")
    logger.info(f"{'Dimensionality:':<25} {tunnel_system.dimension:10d}")
    logger.info("")
    logger.info(f"{'Breakthrough directions and energies:'}")
    logger.info("  Direction    | Highest E barrier (kJ/mol)")
    logger.info("  -----------------------------------------")

    graph_dict_a={}                                                                                                                 #Create a dictionary of all possible transitions (in all 3 directions)
    graph_dict_b={}
    graph_dict_c={}
    for cluster in tunnel_system.cluster_family:
        process_list_a=[]
        process_list_b=[]
        process_list_c=[]
        for process in tunnel_system.process_list:
            if process.start_cluster==cluster:
                process_list_a.append((str(process.end_cluster), float(process.dE), int(process.process_cross_vector[0]), int(process.transition_state)))
                process_list_b.append((str(process.end_cluster), float(process.dE), int(process.process_cross_vector[1]), int(process.transition_state)))
                process_list_c.append((str(process.end_cluster), float(process.dE), int(process.process_cross_vector[2]), int(process.transition_state)))
        graph_dict_a[str(cluster)]=process_list_a
        graph_dict_b[str(cluster)]=process_list_b
        graph_dict_c[str(cluster)]=process_list_c

    directions=tunnel_system.directions
    energies=tunnel_system.direction_breakthroughs*E_step
    directions, energies = directions[np.argsort(energies)], np.sort(energies)

    path_a=None
    path_b=None
    path_c=None
    for idx, direction in enumerate(directions):                                                                                    #Run the minimax algorithm to get the minimum energy barriers along given directions

        Ea_a=E_cutoff
        Ea_b=E_cutoff
        Ea_c=E_cutoff
        for cluster in tunnel_system.cluster_family:
            if direction[0]!=0:
                E_barrier_a, possible_path_a =PBC_minimax.minimax_periodic(graph_dict_a, str(cluster), direction[0])
                if E_barrier_a <= Ea_a+ margin:
                    if path_a is None or len(path_a) > len(possible_path_a):
                        path_a = possible_path_a
                        Ea_a = E_barrier_a
                    else:
                        continue
            else:
                Ea_a = E_cutoff
            if direction[1]!=0:
                E_barrier_b, possible_path_b = PBC_minimax.minimax_periodic(graph_dict_b, str(cluster), direction[1])
                if E_barrier_b <= Ea_b+ margin:
                    if path_b is None or len(path_b) > len(possible_path_b):
                        path_b = possible_path_b
                        Ea_b = E_barrier_b
                    else:
                        continue
            else:
                Ea_b = E_cutoff
            if direction[2]!=0:
                E_barrier_c, possible_path_c = PBC_minimax.minimax_periodic(graph_dict_c, str(cluster), direction[2])
                if E_barrier_c <= Ea_c+ margin:
                    if path_c is None or len(path_c) > len(possible_path_c):
                        path_c = possible_path_c
                        Ea_c = E_barrier_c
                    else:
                        continue                    
            else:
                Ea_c = E_cutoff
        Ea_abc[0] = min(Ea_a, Ea_abc[0])
        Ea_abc[1] = min(Ea_b, Ea_abc[1])
        Ea_abc[2] = min(Ea_c, Ea_abc[2])

        Ea=min(Ea_a, Ea_b, Ea_c)
        logger.info(f"  {str(direction):12}{Ea:10.4f}")

    #"""Save minimum energy pathways in a b c directions (if they exist) for the current tunnel system"""
    if path_a:
        min_E_path_a = open(os.path.join("tunnel data", f"min_E_path_a_tunnel{tunnel_system.ID}.dat"), mode = "w")
        min_E_path_a.write(f"Minimum energy pathway for tunnel system {tunnel_system.ID} in a direction\n")
        min_E_path_a.write(f"{'Start Basin ID':<15}{'End Basin ID':<15}{'TS ID':<15}{'Start Basin E':<15}{'End Basin E':<15}{'TS E':<15}{'Start Coord.':<15}{'End Coord.':<15}{'PBC Crossing':<15}\n")
        for cluster_idx in range(len(path_a)-1):
            cluster1 = int(path_a[cluster_idx][0])
            cluster2 = int(path_a[cluster_idx+1][0])
            TS_ID = int(path_a[cluster_idx+1][1])
            crossing =  int(path_a[cluster_idx+1][2])
            cluster1_frac_coord = data.Cluster.Cluster_list[cluster1-1].center[0]/data.grid[0][0]
            cluster2_frac_coord = data.Cluster.Cluster_list[cluster2-1].center[0]/data.grid[0][0]
            E_transition = data.Transition_state.Transition_state_list[TS_ID-1].E_min
            min_E_path_a.write(f"{cluster1:<15}{cluster2:<15}{TS_ID:<15}{data.Cluster.Cluster_list[cluster1-1].E_min:<15.4f}{data.Cluster.Cluster_list[cluster2-1].E_min:<15.4f}{E_transition:<15.4f}{cluster1_frac_coord:<15.4f}{cluster2_frac_coord:<15.4f}{crossing:<15}\n")
        min_E_path_a.close()

    if path_b:
        min_E_path_b = open(os.path.join("tunnel data", f"min_E_path_b_tunnel{tunnel_system.ID}.dat"), mode = "w")
        min_E_path_b.write(f"Minimum energy pathway for tunnel system {tunnel_system.ID} in b direction\n")
        min_E_path_b.write(f"{'Start Basin ID':<15}{'End Basin ID':<15}{'TS ID':<15}{'Start Basin E':<15}{'End Basin E':<15}{'TS E':<15}{'Start Coord.':<15}{'End Coord.':<15}{'PBC Crossing':<15}\n")
        for cluster_idx in range(len(path_b)-1):
            cluster1 = int(path_b[cluster_idx][0])
            cluster2 = int(path_b[cluster_idx+1][0])
            TS_ID = int(path_b[cluster_idx+1][1])
            crossing =  int(path_b[cluster_idx+1][2])
            cluster1_frac_coord = data.Cluster.Cluster_list[cluster1-1].center[1]/data.grid[1][0]
            cluster2_frac_coord = data.Cluster.Cluster_list[cluster2-1].center[1]/data.grid[1][0]
            E_transition = data.Transition_state.Transition_state_list[TS_ID-1].E_min
            min_E_path_b.write(f"{cluster1:<15}{cluster2:<15}{TS_ID:<15}{data.Cluster.Cluster_list[cluster1-1].E_min:<15.4f}{data.Cluster.Cluster_list[cluster2-1].E_min:<15.4f}{E_transition:<15.4f}{cluster1_frac_coord:<15.4f}{cluster2_frac_coord:<15.4f}{crossing:<15}\n")
        min_E_path_b.close()

    if path_c:
        min_E_path_c = open(os.path.join("tunnel data", f"min_E_path_c_tunnel{tunnel_system.ID}.dat"), mode = "w")
        min_E_path_c.write(f"Minimum energy pathway for tunnel system {tunnel_system.ID} in c direction\n")
        min_E_path_c.write(f"{'Start Basin ID':<15}{'End Basin ID':<15}{'TS ID':<15}{'Start Basin E':<15}{'End Basin E':<15}{'dE':<15}{'Start Coord.':<15}{'End Coord.':<15}{'PBC Crossing':<15}\n")
        for cluster_idx in range(len(path_c)-1):
            cluster1 = int(path_c[cluster_idx][0])
            cluster2 = int(path_c[cluster_idx+1][0])
            TS_ID = int(path_c[cluster_idx+1][1])
            crossing =  int(path_c[cluster_idx+1][2])
            cluster1_frac_coord = data.Cluster.Cluster_list[cluster1-1].center[2]/data.grid[2][0]
            cluster2_frac_coord = data.Cluster.Cluster_list[cluster2-1].center[2]/data.grid[2][0]
            E_transition = data.Transition_state.Transition_state_list[TS_ID-1].E_min
            min_E_path_c.write(f"{cluster1:<15}{cluster2:<15}{TS_ID:<15}{data.Cluster.Cluster_list[cluster1-1].E_min:<15.4f}{data.Cluster.Cluster_list[cluster2-1].E_min:<15.4f}{E_transition:<15.4f}{cluster1_frac_coord:<15.4f}{cluster2_frac_coord:<15.4f}{crossing:<15}\n")
        min_E_path_c.close()

    logger.info("-" * 80)


logger.info(f"TuTraSt analysis performed successfully in {time.process_time()-timer_organization:.2f} s")
logger.info("")

if len(data.Tunnel_system.Tunnel_system_list)!=0 and run_kMC==1:                                                                    #If there is at least one tunnel system and user input for running kMC is 'yes', run the kMC simulations
    """For each temperature, compute reaction rate constants and run kMC simulations"""
    os.makedirs("msd_files", exist_ok=True)                                                                                         #Create a folder where the msd files will be saved
    kappa=0.5                                                                                                                       #Value for an ideal transition state
    mean_grid_size=sum(data.grid_size)*10**(-10)/3                                                                                  #Mean unit cell dimension in meters                    
    data.Cluster_matrix.IDs[data.TS_matrix.IDs != 0] = -1                                                                           #Set the minIDmatrix.Clusters values to -1 for TS points to avoid double counting TS points during the Boltzmann integration
    D_temperatures=np.zeros((len(T_list), len(data.Tunnel_system.Tunnel_system_list), 3, 2))                                        #An array saving the average diffusion coefficients at different temperatures of different tunnel systems and their standard deviations in a b c directions
    D_xyz_temperatures=np.zeros((len(T_list), len(data.Tunnel_system.Tunnel_system_list), 3, 2))                                    #An array saving the average diffusion coefficients at different temperatures of different tunnel systems and their standard deviations in x y z directions
    D_3D_temperatures=np.zeros((len(T_list), len(data.Tunnel_system.Tunnel_system_list), 2))                                        #An array saving the average diffusion coefficient at different temperatures of different tunnel systems and their standard deviations
    for T_idx, T in enumerate(T_list):
        timer_kMC=time.process_time()
        logger.info(f">>> RUNNING kMC SIMULATIONS (T = {T} K) <<<")
        logger.info(f"{'Number of runs:':<25} {N_kMC_runs:>10}")
        logger.info(f"{'Number of steps:':<25} {N_kMC_steps:>10}")              
        Beta=1/(constants.R*T)                      
        prefactor=kappa*np.sqrt(1./(Beta*2.*np.pi*particle_mass))

        Accessible_points_mask=(data.Cluster_matrix.Levels!=0)
        Total_boltzmann_weighted_volume=np.sum(np.exp(-1000*Beta*(data.Energy_matrix[Accessible_points_mask])))                     #Get the boltzmann weighted statistical sum of the accessible volume at temperature T

        for tunnel_idx, tunnel_system in enumerate(data.Tunnel_system.Tunnel_system_list):
            Tunnel_points_mask=(data.Tunnel_matrix.IDs==tunnel_system.ID)                                                           #Get the Boltzmann weighted volume fraction of the tunnel within the system (probability of finding a Li+ in a given tunnel)
            tunnel_system.Boltzmann_weighted_V_fraction=np.sum(np.exp(-1000*Beta*(data.Energy_matrix[Tunnel_points_mask])))/Total_boltzmann_weighted_volume
            logger.info("") 
            logger.info(f"Tunnel system {tunnel_system.ID}, weight = {tunnel_system.Boltzmann_weighted_V_fraction}")
            logger.info("-" * 80)
            logger.info(f"{'Run':>8} {'D_a (m²/s)':>15} {'D_b (m²/s)':>15} {'D_c (m²/s)':>15} {'Δt (s)':>10} {'total (s)':>10}")
            logger.info("-" * 80)   
            for process in tunnel_system.process_list:
                TS_mask=(data.TS_matrix.IDs==process.transition_state)                                                              #Find all points belonging to the current transition state
                TS_sum=np.sum(np.exp(-1000*Beta*(data.Energy_matrix[TS_mask]-data.Energy_matrix[tuple(process.start_point)])))      #Sum the boltzmann weights of all transition state points (offset by the minimum energy in the cluster)
                cluster_mask=(data.Cluster_matrix.IDs==process.start_cluster)                                                       #Do the same for all points belonging to the current starting cluster
                cluster_sum=np.sum(np.exp(-1000*Beta*(data.Energy_matrix[cluster_mask]-data.Energy_matrix[tuple(process.start_point)])))+TS_sum
                k=prefactor*TS_sum/cluster_sum/(mean_grid_size)
                process.k=k
                distance=process.end_point+np.array(data.grid)[:,0]*process.process_cross_vector-process.start_point                #Get the PBC distance between the two clusters involved in the transition
                distance=distance*data.grid_size
                process.distance=distance
            
            """Save the graph networks and boltzmann volume distribution of basins for a given tunnel system at a given temperature in the folder 'tunnel graph networks'"""
            graph_dict={}                                                                                                           #Create a dictionary of all possible transitions (in all 3 directions)
            basin_hist_values = []
            for cluster in tunnel_system.cluster_family:
                cluster_mask = (data.Cluster_matrix.IDs==cluster)
                cluster_Boltzmann_weighted_V_fraction = np.sum(np.exp(-1000*Beta*(data.Energy_matrix[cluster_mask])))/Total_boltzmann_weighted_volume
                basin_hist_values.append(cluster_Boltzmann_weighted_V_fraction)
                process_list = []
                for process in tunnel_system.process_list:
                    if process.start_cluster == cluster:
                        process_list.append((str(process.end_cluster), float(process.k), tuple(int(x) for x in process.process_cross_vector)))
                if len(process_list):                                                                                                #Exclude basins with no transition processes (happens due to 1 point cluster deletion)
                    graph_dict[str(cluster)] = (float(cluster_Boltzmann_weighted_V_fraction), process_list)
            graph_network_file = open(os.path.join("tunnel data", f"Tunnel{tunnel_system.ID}T{T}.json"), mode = "w")
            json.dump(graph_dict, graph_network_file, indent=None)
            graph_network_file.close()

            D, D_xyz, D_3D = kMC.kMC(tunnel_system, N_kMC_runs, N_kMC_steps, T, tunnel_idx)                                          #Run kMC simulations
            D_temperatures[T_idx, tunnel_idx, :, :] = D*tunnel_system.Boltzmann_weighted_V_fraction                                  #Add the weighted D contribution of a given tunnel system to the list
            D_xyz_temperatures[T_idx, tunnel_idx, :, :] = D_xyz*tunnel_system.Boltzmann_weighted_V_fraction
            D_3D_temperatures[T_idx, tunnel_idx, :] = D_3D*tunnel_system.Boltzmann_weighted_V_fraction
            logger.info("-" * 80)
            logger.info(f"Average self diffusion coefficient in a b c: \nD_a = {D[0,0]:>.3e} ± {D[0,1]:>.2e}  D_b = {D[1,0]:>.3e} ± {D[1,1]:>.2e}  D_c = {D[2,0]:>.3e} ± {D[2,1]:>.2e}")
            logger.info(f"Average self diffusion coefficient in x y z: \nD_x = {D_xyz[0,0]:>.3e} ± {D_xyz[0,1]:>.2e}  D_y = {D_xyz[1,0]:>.3e} ± {D_xyz[1,1]:>.2e}  D_z = {D_xyz[2,0]:>.3e} ± {D_xyz[2,1]:>.2e}")
            logger.info(f"Average Isotropic self diffusion coefficient: \nD = {D_3D[0]:>.3e} ± {D_3D[1]:>.2e}")
            logger.info("-"*80)          
            logger.info(f"Completed kMC simulations for tunnel system {tunnel_system.ID} at T = {T} K. \nTotal time elapsed: {time.process_time()-timer_kMC:.2f} s")
            logger.info("")
        logger.info("-" * 80)
        logger.info(f"Joint diffusion coefficients across all tunnel systems at T = {T} K:")
        logger.info("In the a b c directions:")
        logger.info(f"D_a = {np.sum(D_temperatures[T_idx, :, 0, 0]):>.3e} ± {np.sum(D_temperatures[T_idx, :, 0, 1]):>.2e}  D_b = {np.sum(D_temperatures[T_idx, :, 1, 0]):>.3e} ± {np.sum(D_temperatures[T_idx, :, 1, 1]):>.2e}  D_c = {np.sum(D_temperatures[T_idx, :, 2, 0]):>.3e} ± {np.sum(D_temperatures[T_idx, :, 2, 1]):>.2e}")
        logger.info("In the x y z directions:")
        logger.info(f"D_x = {np.sum(D_xyz_temperatures[T_idx, :, 0, 0]):>.3e} ± {np.sum(D_xyz_temperatures[T_idx, :, 0, 1]):>.2e}  D_y = {np.sum(D_xyz_temperatures[T_idx, :, 1, 0]):>.3e} ± {np.sum(D_xyz_temperatures[T_idx, :, 1, 1]):>.2e}  D_z = {np.sum(D_xyz_temperatures[T_idx, :, 2, 0]):>.3e} ± {np.sum(D_xyz_temperatures[T_idx, :, 2, 1]):>.2e}")
        logger.info("Isotropic self diffusion coefficient:")
        logger.info(f"D = {np.sum(D_3D_temperatures[T_idx, :, 0]):>.3e} ± {np.sum(D_3D_temperatures[T_idx, :, 1]):>.2e}")
        logger.info("")

logger.info(">>> SUMMARY <<<")
logger.info("-" * 80)
logger.info("Geometric descriptors")
logger.info("-" * 80)
logger.info("Volume")
logger.info(f"  {'Total:':<15} {Total_volume:10.4f} Å³ {Total_volume_fraction:10.4f} %")
logger.info(f"  {'Accessible:':<15} {Accessible_volume:10.4f} Å³ {Accessible_volume/Total_volume*Total_volume_fraction:10.4f} %")
logger.info("")
logger.info("Surface area")
logger.info(f"  {'Total:':<15} {Total_area:10.4f} Å² {Total_area/(Total_volume/Total_volume_fraction):10.4f} Å²/Å³")
logger.info(f"  {'Accessible:':<15} {Accessible_area:10.4f} Å² {Accessible_area/(Total_volume/Total_volume_fraction):10.4f} Å²/Å³")
logger.info("-" * 80)
logger.info("Breakthrough energies")
logger.info("-" * 80)
if breakthrough_level[0]!=level_max and Ea_abc[0]<E_cutoff:
    logger.info(f"{'Direction a:    level'} {breakthrough_level[0]:>3}    E = {breakthrough_level[0]*E_step:10.3f} kJ/mol     ΔE = {Ea_abc[0]:10.3f} kJ/mol")
else:
    logger.info(f"{'Direction a:    level'} {'/':>3}    E = {'/':10} kJ/mol     ΔE = {'/':10} kJ/mol")
if breakthrough_level[1]!=level_max and Ea_abc[0]<E_cutoff: 
    logger.info(f"{'Direction b:    level'} {breakthrough_level[1]:>3}    E = {breakthrough_level[1]*E_step:10.3f} kJ/mol     ΔE = {Ea_abc[1]:10.3f} kJ/mol")
else:
    logger.info(f"{'Direction b:    level'} {'/':>3}    E = {'/':10} kJ/mol     ΔE = {'/':10} kJ/mol")
if breakthrough_level[2]!=level_max and Ea_abc[0]<E_cutoff: 
    logger.info(f"{'Direction c:    level'} {breakthrough_level[2]:>3}    E = {breakthrough_level[2]*E_step:10.3f} kJ/mol     ΔE = {Ea_abc[2]:10.3f} kJ/mol")
else:
    logger.info(f"{'Direction c:    level'} {'/':>3}    E = {'/':10} kJ/mol     ΔE = {'/':10} kJ/mol")
logger.info("-" * 80)
    
if len(data.Tunnel_system.Tunnel_system_list)!=0 and run_kMC==1:
    logger.info("Diffusion coefficients")
    logger.info("-" * 80)
    logger.info(f"{'T (K)':>8}  {'D_a (m²/s) ± σ':>22}  {'D_b (m²/s) ± σ':>22}  {'D_c (m²/s) ± σ':>21}")
    logger.info("-" * 80)
    for i, T in enumerate(T_list):
        logger.info(f"{T:8.1f}  "
                    f"{np.sum(D_temperatures[i, :, 0, 0]):12.3e} ± {np.sum(D_temperatures[i, :, 0, 1]):8.2e}"
                    f"{np.sum(D_temperatures[i, :, 1, 0]):12.3e} ± {np.sum(D_temperatures[i, :, 1, 1]):8.2e}"
                    f"{np.sum(D_temperatures[i, :, 2, 0]):12.3e} ± {np.sum(D_temperatures[i, :, 2, 1]):8.2e}")
    logger.info("-" * 80)
elif len(data.Tunnel_system.Tunnel_system_list)==0 and run_kMC==1:
    logger.info("Diffusion coefficients")
    logger.info("-" * 80)
    logger.info(f"{'T (K)':>8}  {'D_a (m²/s) ± σ':>22}  {'D_b (m²/s) ± σ':>22}  {'D_c (m²/s) ± σ':>21}")
    logger.info("-" * 80)
    for i, T in enumerate(T_list):
        logger.info(f"{T:8.1f}  "
                    f"{0:12.3e} ± {0:8.2e}"
                    f"{0:12.3e} ± {0:8.2e}"
                    f"{0:12.3e} ± {0:8.2e}")
    logger.info("-" * 80)

logger.info("")
logger.info(f"{'Total runtime:':<30} {time.process_time()-timer:.2f} s")
logger.info(f"{'Log file:':<30} output.log")
logger.info("-" * 80)
logger.info("")
logger.info("Normal termination of EnGOAT.")
logger.info("*" * 80)
console_handler.close()
file_handler.close()

"""Save simulation data"""

os.makedirs("numpy matrices", exist_ok=True)                                                                                        #Create a folder where the numpy matrices generated during the TuTraSt analysis will be stored
np.save(os.path.join("numpy matrices", "Level_matrix.npy"), data.level_matrix)
np.save(os.path.join("numpy matrices", "Basin_matrix.npy"), data.Cluster_matrix.IDs)
np.save(os.path.join("numpy matrices", "TS_matrix.npy"), data.TS_matrix.IDs)
np.save(os.path.join("numpy matrices", "Tunnel_matrix.npy"), data.Tunnel_matrix.IDs)
np.save(os.path.join("numpy matrices", "Energy_matrix.npy"), data.Energy_matrix)
                                                       
for tunnel_system in data.Tunnel_system.Tunnel_system_list:                                                                         #Save histogram files
    histogram_file=open(os.path.join("tunnel data", f"tunnel{tunnel_system.ID}_hist.dat"), mode="w")
    histogram_file.write(f"#Level No.\tLevel E [kJ/mol]\tN points\n")
    for levelm1, N_points in enumerate(tunnel_system.histogram):
        histogram_file.write(f"{levelm1+1:>9d}\t{(levelm1+1)*E_step:>13.2f}\t{int(N_points):>12d}\n")
    histogram_file.close()

a=np.array([data.grid[0][1], data.grid[0][2], data.grid[0][3]])      
b=np.array([data.grid[1][1], data.grid[1][2], data.grid[1][3]])
c=np.array([data.grid[2][1], data.grid[2][2], data.grid[2][3]])                               
voxel_V=abs(np.dot(a, np.cross(b, c))) 
os.makedirs("TuTraSt data", exist_ok=True)                                                                                          #Create a folder where the numpy matrices generated during the TuTraSt analysis will be stored
cluster_file=open(os.path.join("TuTraSt data", "basin_data.dat"), mode="w")
cluster_file.write(f"#{'Basin ID':20} {'Tunnel system':20} {'Center (a, b, c)':20} {'E_min [kJ]':20} {'Volume [Å^3]':20}")
for T in T_list:
    cluster_file.write(f" V_Boltz T = {T:7.2f}   ")
cluster_file.write(f"\n")
for cluster in data.Cluster.Cluster_list:
    if cluster.active:
        ID = int(cluster.ID)
        a_center, b_center, c_center = cluster.center
        Tunnelsystem = "/"
        for tunnel_system in data.Tunnel_system.Tunnel_system_list:
            if ID in tunnel_system.cluster_family:
                Tunnelsystem = tunnel_system.ID
        Volume = np.count_nonzero(data.Cluster_matrix.IDs == ID)*voxel_V
        cluster_file.write(f"{ID:<20} {Tunnelsystem:<20} ({int(a_center):4}, {int(b_center):4}, {int(c_center):4})    {cluster.E_min:<20.4f} {Volume:<20.4f}")
        for T in T_list:
            Beta=1/(constants.R*T)
            cluster_mask = (data.Cluster_matrix.IDs==ID)
            cluster_Boltzmann_weighted_V_fraction = np.sum(np.exp(-1000*Beta*(data.Energy_matrix[cluster_mask])))
            cluster_file.write(f" {cluster_Boltzmann_weighted_V_fraction:<20.6f}  ")
        cluster_file.write(f"\n")
cluster_file.close()

TS_file=open(os.path.join("TuTraSt data", "TS_data.dat"), mode="w")
TS_file.write(f"{'#ID':<15}{'E_min [kJ/mol]':<15}{'Basin 1':<15}{'Basin 2':<15}\n")
for TS in data.Transition_state.Transition_state_list:
    C1_ID, C2_ID = TS.clusters
    TS_file.write(f"{int(TS.ID):<15d}{TS.E_min:<15.4f}{int(C1_ID):<15}{int(C2_ID):<15}\n")
TS_file.close()

for tunnel_system in data.Tunnel_system.Tunnel_system_list:
    tunnel_file=open(os.path.join("TuTraSt data", f"tunnel{tunnel_system.ID}_data.dat"), mode="w")
    tunnel_file.write(f"#Basin family:\n")
    family_string="{"
    for clusterID in tunnel_system.cluster_family:
        family_string+=f"{int(clusterID)}, "
    family_string=family_string[:-2]+"}"
    tunnel_file.write(f"{family_string}\n\n")
    tunnel_file.write(f"#Process list\n")
    tunnel_file.write(f"#Start basin \t end basin \tTransition state \t Energy barrier \t Start point [a, b, c] \tEnd point [a, b, c] \t process vector (a, b, c)\n")
    for process in tunnel_system.process_list:
        a_start, b_start, c_start = process.start_point
        a_end, b_end, c_end = process.end_point
        a_cross, b_cross, c_cross = process.process_cross_vector
        start_string=f"[{a_start}, {b_start}, {c_start}]"
        end_string=f"[{a_end}, {b_end}, {c_end}]"
        cross_string=f"({a_cross}, {b_cross}, {c_cross})"
        tunnel_file.write(f"{process.start_cluster:>14} \t {process.end_cluster:>10} \t{process.transition_state:>15} \t {process.dE:>12.6f} \t {start_string:>20} {end_string:>20} {cross_string:>20} \n")
    tunnel_file.close()

"""
WARNING: basin==cluster in the code!

IDEAS:

make tutrast analysis at intermediate E levels
"""

