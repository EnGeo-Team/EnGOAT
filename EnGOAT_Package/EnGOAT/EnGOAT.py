import numpy as np
from scipy import constants
import time
import os
import json
import logging
from datetime import datetime
from pathlib import Path
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

ROOT = Path(__file__).resolve().parent  #ROOT directory; change once GUI is made -> set to directory of the cube file

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
logger.info(f"{'Origin (x,y,z):':<25} ({data.grid["origin"][0]:>8.3f}, {data.grid["origin"][1]:>8.3f}, {data.grid["origin"][2]:>8.3f})")
logger.info("-" * 80)
logger.info(f"{'Grid points (Na,Nb,Nc):':<25} ({data.grid["shape"][0]:>8}, {data.grid["shape"][1]:>8}, {data.grid["shape"][2]:>8})")
logger.info(f"{'a-axis vector:':<25} ({data.grid["vectors"][0][0]:>8.3f}, {data.grid["vectors"][0][1]:>8.3f}, {data.grid["vectors"][0][2]:>8.3f})")
logger.info(f"{'b-axis vector:':<25} ({data.grid["vectors"][1][0]:>8.3f}, {data.grid["vectors"][1][1]:>8.3f}, {data.grid["vectors"][1][2]:>8.3f})")
logger.info(f"{'c-axis vector:':<25} ({data.grid["vectors"][2][0]:>8.3f}, {data.grid["vectors"][2][1]:>8.3f}, {data.grid["vectors"][2][2]:>8.3f})")
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

logger.info("-" * 80)
logger.info(f"Potential grid analyzed successfully in {time.process_time()-timer_analysis:.2f} s")
logger.info("")

logger.info(">>> TUTRAST ANALYSIS <<<")
timer_organization=time.process_time()                                                                                

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

"""Organize transition state points into transition state surfaces, and identify tunnel systems. Compute the minimum energy pathways of each tunnel system. Save all data into dictionaries and dump them into a .json file"""                      
Organize_TuTraSt.Organize_transition_states()                                                                                       #A python subroutine that organizes transition state points into transition states, fills the TS_matrix.IDs, and sets the minIDmatrix.Clusters values of TS points to -1
Organize_TuTraSt.Organize_tunnel_systems(level_max)                                                                                 #A python subroutine that organizes recorded breakthroughs into cluster-family-wide tunnel systems          

data.Cluster_matrix.IDs[data.TS_matrix.IDs != 0] = -1                                                                               #Set the minIDmatrix.Clusters values to -1 for TS points to avoid double counting TS points

A_tot, A_acc, V_tot, V_acc, V_UC = get_topological_descriptors.get_topological_descriptors(int(N_levels), E_step, T_list)

V_B_tot = {}                                                                                                                        #Total boltzmann weighted volume of the system
for T in T_list:
    Beta=1/(constants.R*T) 
    Accessible_points_mask=(data.Cluster_matrix.Levels!=0)
    V_B_tot[T] = np.sum(np.exp(-1000*Beta*(data.Energy_matrix[Accessible_points_mask])))


#
# Save all data into dictionaries to be stored as .json files
#


metadata = {
    "E_step": E_step,
    "E_cutoff": E_cutoff,
    "N_levels": int(E_cutoff/E_step),
    "origin": data.grid["origin"], 
    "grid_shape": data.grid["shape"].tolist(),
    "grid_vectors": data.grid["vectors"].tolist(),
    "grid_spacing": np.linalg.norm(data.grid["vectors"], axis=1).tolist(),
    "V_voxel": abs(np.linalg.det(data.grid["vectors"])),
    "energy_unit": "kJ/mol",
    "temperature_list": T_list
}

unit_cell = {"A_tot": A_tot, "A_acc": A_acc, "V_tot": V_tot, "V_acc": V_acc, "V_UC": V_UC, "V_B_tot": V_B_tot}

atoms = data.atoms

basin_data = {}
for B in data.Cluster.Cluster_list:
    if B.active:
        basin_data[int(B.ID)] = {"center": tuple(B.center), "E_min": B.E_min, "V": B.V, "A": B.A, "V_rel": B.V_rel, "A_rel": B.A_rel, "histogram": B.histogram, "V_B": B.V_B}

TS_data = {}
for TS in data.Transition_state.Transition_state_list:
    if TS.clusters[1] == -1:                                                                                                        #Treat the same-cluster transition states (make it easier for later)
        basins = (TS.clusters[0], TS.clusters[0])
    else:
        basins = TS.clusters
    TS_data[int(TS.ID)] = {"basins": basins, "start_center": tuple(basin_data[basins[0]]["center"]), "end_center": tuple(basin_data[basins[1]]["center"]), "cross_vector": TS.Process_cross_vector.tolist(), "E_min": TS.E_min, 
                      "V": TS.V, "A": TS.A, "V_rel": TS.V_rel, "A_rel": TS.A_rel, "histogram": TS.histogram, "V_B": TS.V_B}

tunnel_systems = {}
for tunnel in data.Tunnel_system.Tunnel_system_list:
    tunnel_systems[tunnel.ID] = {"basin_list": list(tunnel.cluster_family), "TS_list": list(tunnel.transition_state_list), "E_min": tunnel.height, "dimensionality": tunnel.dimension, 
                      "V": tunnel.V, "A": tunnel.A, "V_rel": tunnel.V_rel, "A_rel": tunnel.A_rel, "histogram": tunnel.histogram, "V_B": tunnel.V_B}
    for B_ID in tunnel.cluster_family:
        basin_data[B_ID]["group_type"] = "tunnel_system"
        basin_data[B_ID]["group_ID"] = tunnel.ID
    
    for TS_ID in tunnel.transition_state_list:
        TS_data[TS_ID]["group_type"] = "tunnel_system"
        TS_data[TS_ID]["group_ID"] = tunnel.ID

isolated_groups = {}
for iso in data.Isolated_group.Isolated_group_list:
    isolated_groups[iso.ID] = {"basin_list": list(iso.cluster_family), "TS_list": list(iso.transition_state_list), "E_min": iso.height, 
                      "V": iso.V, "A": iso.A, "V_rel": iso.V_rel, "A_rel": iso.A_rel, "histogram": iso.histogram, "V_B": iso.V_B}
    for B_ID in iso.cluster_family:
        basin_data[B_ID]["group_type"] = "isolated_group"
        basin_data[B_ID]["group_ID"] = iso.ID
    
    for TS_ID in iso.transition_state_list:
        TS_data[TS_ID]["group_type"] = "isolated_group"
        TS_data[TS_ID]["group_ID"] = iso.ID

logger.info(f"{'Total volume:':<25} {V_tot:10.4f} Å³ {V_tot/V_UC:10.4f} %")
logger.info(f"{'Accessible volume:':<25} {V_acc:10.4f} Å³ {V_acc/V_UC:10.4f} %")
logger.info(f"{'Total surface area:':<25} {A_tot:10.4f} Å² {A_tot/V_UC:10.4f} Å²/Å³")
logger.info(f"{'Accessible surface area:':<25} {A_acc:10.4f} Å² {A_acc/V_tot:10.4f} Å²/Å³")
logger.info("")
logger.info(f"{len(data.Tunnel_system.Tunnel_system_list)} tunnel systems were identified:")
logger.info("-" * 80)

os.makedirs("Graph_networks", exist_ok=True)                                                                                        #Create a folder where tunnel systems will be saved as graph networks, where basins will be defined with their corresponding boltzmann weighted probabilities and a process list of possible transitions defined each with the corresponding rate
margin = E_step                                                                                                                     #E tolerance for finding the best energy path along a direction (the paths that differ for at most this energy will be evaluated based on length)
for tunnel_system in data.Tunnel_system.Tunnel_system_list:
    logger.info(f"Tunnel system #{tunnel_system.ID}")
    logger.info("-" * 80)
    logger.info(f"{'Volume:':<25} {tunnel_system.V:10.4f} Å³ {tunnel_system.V/V_UC:10.4f} %")
    logger.info(f"{'Surface area:':<25} {tunnel_system.A:10.4f} Å² {tunnel_system.A/V_UC:10.4f} Å²/Å³")
    logger.info(f"{'Lowest energy point:':<25} {tunnel_system.height:10.4f} kJ/mol")
    logger.info(f"{'Dimensionality:':<25} {tunnel_system.dimension:10d}")
    logger.info("")
    logger.info(f"{'Breakthrough directions and energies:'}")
    logger.info("  Direction    | Highest E barrier (kJ/mol)")
    logger.info("  -----------------------------------------")

    graph_dict_a={}                                                                                                                 #Create a dictionary of all possible transitions (in all 3 directions)
    graph_dict_b={}
    graph_dict_c={}
    graph_dict = {}                                                                                                                 #Whole graph network (for saving in a file)
    for cluster in tunnel_system.cluster_family:
        process_list_a=[]
        process_list_b=[]
        process_list_c=[]
        process_list = []
        for process in tunnel_system.process_list:
            if process.start_cluster==cluster:
                process_list_a.append((str(process.end_cluster), float(process.dE), int(process.process_cross_vector[0]), int(process.transition_state)))
                process_list_b.append((str(process.end_cluster), float(process.dE), int(process.process_cross_vector[1]), int(process.transition_state)))
                process_list_c.append((str(process.end_cluster), float(process.dE), int(process.process_cross_vector[2]), int(process.transition_state)))
                process_list.append((str(process.end_cluster), float(process.dE), tuple(process.process_cross_vector), str(process.transition_state)))
        graph_dict_a[str(cluster)]=process_list_a
        graph_dict_b[str(cluster)]=process_list_b
        graph_dict_c[str(cluster)]=process_list_c
        graph_dict[str(cluster)] = process_list
    
    filename = os.path.join("Graph_networks", f"TunnelSystem{tunnel_system.ID}_Eenergy_graph.json")
    with open(filename, "w", encoding="utf-8") as f:
        json.dump(
            graph_dict,
            f,
            default=lambda x: x.item(),
            indent=4
        )

    directions=tunnel_system.directions
    energies=tunnel_system.direction_breakthroughs*E_step
    directions, energies = directions[np.argsort(energies)], np.sort(energies)

    #Run the minimax algorithm to get the minimum energy barriers along given directions
    Ea_a=E_cutoff
    Ea_b=E_cutoff
    Ea_c=E_cutoff
    path_a=None
    path_b=None
    path_c=None
    a_break = any(d[0] != 0 for d in directions)
    b_break = any(d[1] != 0 for d in directions)
    c_break = any(d[2] != 0 for d in directions)

    if a_break:   
        for cluster in tunnel_system.cluster_family:
            E_barrier_a, possible_path_a =PBC_minimax.minimax_periodic(graph_dict_a, str(cluster), 1)
            if E_barrier_a <= Ea_a+ margin:
                if path_a is None or len(path_a) > len(possible_path_a):
                    path_a = possible_path_a
                    Ea_a = E_barrier_a
                else:
                    continue
    if b_break:
        for cluster in tunnel_system.cluster_family:
            E_barrier_b, possible_path_b = PBC_minimax.minimax_periodic(graph_dict_b, str(cluster), 1)
            if E_barrier_b <= Ea_b+ margin:
                if path_b is None or len(path_b) > len(possible_path_b):
                    path_b = possible_path_b
                    Ea_b = E_barrier_b
                else:
                    continue
    if c_break:
        for cluster in tunnel_system.cluster_family:
            E_barrier_c, possible_path_c = PBC_minimax.minimax_periodic(graph_dict_c, str(cluster), 1)
            if E_barrier_c <= Ea_c+ margin:
                if path_c is None or len(path_c) > len(possible_path_c):
                    path_c = possible_path_c
                    Ea_c = E_barrier_c

    for idx, direction in enumerate(directions):
        E_abc = np.array([E_cutoff, E_cutoff, E_cutoff])
        if direction[0]:
            E_abc[0] = Ea_a
        if direction[1]:
            E_abc[1] = Ea_b
        if direction[2]:
            E_abc[2] = Ea_c

        Ea=min(E_abc)
        logger.info(f"  {str(direction):12}{Ea:10.4f}")

    #"""Save minimum energy pathways in a b c directions (if they exist) for the current tunnel system"""
    MEPs = {}
    if path_a:
        path = []
        for cluster_idx in range(len(path_a)-1):
            cluster1 = int(path_a[cluster_idx][0])
            cluster2 = int(path_a[cluster_idx+1][0])
            TS_ID = int(path_a[cluster_idx+1][1])
            crossing =  int(path_a[cluster_idx+1][2])
            path.append({"start_basin":cluster1, "end_basin": cluster2, "transition_state": TS_ID, "crossing": crossing})
        MEPs["a"] = {"path": path}
    if path_b:
        path = []
        for cluster_idx in range(len(path_b)-1):
            cluster1 = int(path_b[cluster_idx][0])
            cluster2 = int(path_b[cluster_idx+1][0])
            TS_ID = int(path_b[cluster_idx+1][1])
            crossing =  int(path_b[cluster_idx+1][2])
            path.append({"start_basin":cluster1, "end_basin": cluster2, "transition_state": TS_ID, "crossing": crossing})
        MEPs["b"] = {"path": path}
    if path_c:
        path = []
        for cluster_idx in range(len(path_c)-1):
            cluster1 = int(path_c[cluster_idx][0])
            cluster2 = int(path_c[cluster_idx+1][0])
            TS_ID = int(path_c[cluster_idx+1][1])
            crossing =  int(path_c[cluster_idx+1][2])
            path.append({"start_basin":cluster1, "end_basin": cluster2, "transition_state": TS_ID, "crossing": crossing})
        MEPs["c"] = {"path": path}

    tunnel_systems[tunnel_system.ID]["MEPs"] = MEPs

    logger.info("-" * 80)


logger.info(f"TuTraSt analysis performed successfully in {time.process_time()-timer_organization:.2f} s")
logger.info("")

#Create a master dictionary storing all metadata collected during the TuTraSt analysis
EnGOAT_data = {"metadata": metadata, "unit_cell": unit_cell, "atoms": atoms, "basin_data": basin_data, "TS_data": TS_data, "tunnel_systems": tunnel_systems, "isolated_groups": isolated_groups, "kMC_data": None}

Diffusion = None
"""Run the kMC simulations"""
if len(data.Tunnel_system.Tunnel_system_list)!=0 and run_kMC==1:                                                                    #If there is at least one tunnel system and user input for running kMC is 'yes', run the kMC simulations
    
    Diffusion = {}

    """For each temperature, compute reaction rate constants and run kMC simulations"""
    kappa=0.5                                                                                                                       #Value for an ideal transition state
    mean_voxel_size=sum(data.grid["voxel_size"])*10**(-10)/3                                                                        #Mean unit cell dimension in meters                    
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

        Total_boltzmann_weighted_volume = V_B_tot[T]

        Diffusion[T] = {
            "D_tot": {},
            "tunnel_systems": {}
        }

        for tunnel_idx, tunnel_system in enumerate(data.Tunnel_system.Tunnel_system_list):

            #Get the Boltzmann weighted volume fraction of the tunnel within the system (probability of finding a Li+ in a given tunnel)
            tunnel_system.weight = tunnel_system.V_B[T]/Total_boltzmann_weighted_volume
            
            logger.info("") 
            logger.info(f"Tunnel system {tunnel_system.ID}, weight = {tunnel_system.weight}")
            logger.info("-" * 80)
            logger.info(
    f"{'Run':>8} "
    f"{'D_a (cm^2/s)':>15} "
    f"{'D_b (cm^2/s)':>15} "
    f"{'D_c (cm^2/s)':>15} "
    f"{'dt (s)':>10} "
    f"{'total (s)':>10}"
)
            logger.info("-" * 80)   
            for process in tunnel_system.process_list:
                TS_mask=(data.TS_matrix.IDs==process.transition_state)                                                              #Find all points belonging to the current transition state
                TS_sum=np.sum(np.exp(-1000*Beta*(data.Energy_matrix[TS_mask]-data.Energy_matrix[tuple(process.start_point)])))      #Sum the boltzmann weights of all transition state points (offset by the minimum energy in the cluster)
                cluster_mask=(data.Cluster_matrix.IDs==process.start_cluster)                                                       #Do the same for all points belonging to the current starting cluster
                cluster_sum=np.sum(np.exp(-1000*Beta*(data.Energy_matrix[cluster_mask]-data.Energy_matrix[tuple(process.start_point)])))+TS_sum
                k=prefactor*TS_sum/cluster_sum/(mean_voxel_size)
                process.k=k
                distance=process.end_point+data.grid["shape"]*process.process_cross_vector-process.start_point                      #Get the PBC distance between the two clusters involved in the transition
                distance=distance*data.grid["voxel_size"]
                process.distance=distance
            
            """Save the graph networks of basins for a given tunnel system at a given temperature in the folder 'Tunnel_data'"""
            graph_dict={}                                                                                                           #Create a dictionary of all possible transitions (in all 3 directions)
            for cluster in tunnel_system.cluster_family:
                cluster_Boltzmann_weighted_V_fraction = data.Cluster.Cluster_list[int(cluster)-1].V_B[T]/Total_boltzmann_weighted_volume
                process_list = []
                for process in tunnel_system.process_list:
                    if process.start_cluster == cluster:
                        process_list.append((str(process.end_cluster), float(process.k), tuple(int(x) for x in process.process_cross_vector)))
                graph_dict[str(cluster)] = process_list
            filename = os.path.join(
                "Graph_networks",
                f"TunnelSystem{tunnel_system.ID}_Rate_graph_T{T}.json"
            )            
            with open(filename, "w", encoding="utf-8") as f:
                json.dump(
                    graph_dict,
                    f,
                    default=lambda x: x.item(),
                    indent=4
                )

            kMC_results = kMC.kMC(tunnel_system, N_kMC_runs, N_kMC_steps, T, tunnel_idx)                                            #Run kMC simulations

            Diffusion[T]["tunnel_systems"][tunnel_system.ID] = {
                "weight": float(tunnel_system.weight),
                "directions": kMC_results
            }

            D_temperatures[T_idx, tunnel_idx] = tunnel_system.weight * np.array([
                [kMC_results["a"]["D"], kMC_results["a"]["sd"]],
                [kMC_results["b"]["D"], kMC_results["b"]["sd"]],
                [kMC_results["c"]["D"], kMC_results["c"]["sd"]],
            ])                                                                                                                      #Add the weighted D contribution of a given tunnel system to the list
            D_xyz_temperatures[T_idx, tunnel_idx] = tunnel_system.weight * np.array([
                [kMC_results["x"]["D"], kMC_results["x"]["sd"]],
                [kMC_results["y"]["D"], kMC_results["y"]["sd"]],
                [kMC_results["z"]["D"], kMC_results["z"]["sd"]],
            ])

            D_3D_temperatures[T_idx, tunnel_idx] = tunnel_system.weight * np.array([
                kMC_results["3D"]["D"],
                kMC_results["3D"]["sd"],
            ])
            logger.info("-" * 80)
            logger.info(
                "Average self diffusion coefficient in a b c:\n"
                f"D_a = {kMC_results['a']['D']:.3e} ± {kMC_results['a']['sd']:.2e}  "
                f"D_b = {kMC_results['b']['D']:.3e} ± {kMC_results['b']['sd']:.2e}  "
                f"D_c = {kMC_results['c']['D']:.3e} ± {kMC_results['c']['sd']:.2e}"
            )

            logger.info(
                "Average self diffusion coefficient in x y z:\n"
                f"D_x = {kMC_results['x']['D']:.3e} ± {kMC_results['x']['sd']:.2e}  "
                f"D_y = {kMC_results['y']['D']:.3e} ± {kMC_results['y']['sd']:.2e}  "
                f"D_z = {kMC_results['z']['D']:.3e} ± {kMC_results['z']['sd']:.2e}"
            )

            logger.info(
                "Average isotropic self diffusion coefficient:\n"
                f"D = {kMC_results['3D']['D']:.3e} ± {kMC_results['3D']['sd']:.2e}"
            )
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

        Diffusion[T]["D_tot"] = {
            "a": {
                "D": float(np.sum(D_temperatures[T_idx, :, 0, 0])),
                "sd": float(np.sum(D_temperatures[T_idx, :, 0, 1])),
            },
            "b": {
                "D": float(np.sum(D_temperatures[T_idx, :, 1, 0])),
                "sd": float(np.sum(D_temperatures[T_idx, :, 1, 1])),
            },
            "c": {
                "D": float(np.sum(D_temperatures[T_idx, :, 2, 0])),
                "sd": float(np.sum(D_temperatures[T_idx, :, 2, 1])),
            },
            "x": {
                "D": float(np.sum(D_xyz_temperatures[T_idx, :, 0, 0])),
                "sd": float(np.sum(D_xyz_temperatures[T_idx, :, 0, 1])),
            },
            "y": {
                "D": float(np.sum(D_xyz_temperatures[T_idx, :, 1, 0])),
                "sd": float(np.sum(D_xyz_temperatures[T_idx, :, 1, 1])),
            },
            "z": {
                "D": float(np.sum(D_xyz_temperatures[T_idx, :, 2, 0])),
                "sd": float(np.sum(D_xyz_temperatures[T_idx, :, 2, 1])),
            },
            "3D": {
                "D": float(np.sum(D_3D_temperatures[T_idx, :, 0])),
                "sd": float(np.sum(D_3D_temperatures[T_idx, :, 1])),
            },
        }
EnGOAT_data["kMC_data"] = Diffusion

logger.info(">>> SUMMARY <<<")
logger.info("-" * 80)
logger.info("Geometric descriptors")
logger.info("-" * 80)
logger.info("Volume")
logger.info(f"  {'Total:':<15} {V_tot:10.4f} Å³ {V_tot/V_UC:10.4f} %")
logger.info(f"  {'Accessible:':<15} {V_acc:10.4f} Å³ {V_acc/V_UC:10.4f} %")
logger.info("")
logger.info("Surface area")
logger.info(f"  {'Total:':<15} {A_tot:10.4f} Å² {A_tot/V_UC:10.4f} Å²/Å³")
logger.info(f"  {'Accessible:':<15} {A_acc:10.4f} Å² {A_acc/V_UC:10.4f} Å²/Å³")
logger.info("-" * 80)
logger.info("Breakthrough energies")
logger.info("-" * 80)
if breakthrough_level[0]!=level_max and Ea_a<E_cutoff:
    logger.info(f"{'Direction a:    level'} {breakthrough_level[0]:>3}    E = {breakthrough_level[0]*E_step:10.3f} kJ/mol     ΔE = {Ea_a:10.3f} kJ/mol")
else:
    logger.info(f"{'Direction a:    level'} {'/':>3}    E = {'/':10} kJ/mol     ΔE = {'/':10} kJ/mol")
if breakthrough_level[1]!=level_max and Ea_b<E_cutoff: 
    logger.info(f"{'Direction b:    level'} {breakthrough_level[1]:>3}    E = {breakthrough_level[1]*E_step:10.3f} kJ/mol     ΔE = {Ea_b:10.3f} kJ/mol")
else:
    logger.info(f"{'Direction b:    level'} {'/':>3}    E = {'/':10} kJ/mol     ΔE = {'/':10} kJ/mol")
if breakthrough_level[2]!=level_max and Ea_c<E_cutoff: 
    logger.info(f"{'Direction c:    level'} {breakthrough_level[2]:>3}    E = {breakthrough_level[2]*E_step:10.3f} kJ/mol     ΔE = {Ea_c:10.3f} kJ/mol")
else:
    logger.info(f"{'Direction c:    level'} {'/':>3}    E = {'/':10} kJ/mol     ΔE = {'/':10} kJ/mol")
logger.info("-" * 80)
    
if len(data.Tunnel_system.Tunnel_system_list)!=0 and run_kMC==1:
    logger.info("Diffusion coefficients")
    logger.info("-" * 80)
    logger.info(f"{'T (K)':>8}  {'D_a (cm²/s) ± σ':>22}  {'D_b (cm²/s) ± σ':>22}  {'D_c (cm²/s) ± σ':>21}")
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
    logger.info(f"{'T (K)':>8}  {'D_a (cm²/s) ± σ':>22}  {'D_b (cm²/s) ± σ':>22}  {'D_c (cm²/s) ± σ':>21}")
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

os.makedirs("NumPy_matrices", exist_ok=True)                                                                                        #Create a folder where the numpy matrices generated during the TuTraSt analysis will be stored
np.save(os.path.join("NumPy_matrices", "Level_matrix.npy"), data.level_matrix)
np.save(os.path.join("NumPy_matrices", "Basin_matrix.npy"), data.Cluster_matrix.IDs)
np.save(os.path.join("NumPy_matrices", "TS_matrix.npy"), data.TS_matrix.IDs)
np.save(os.path.join("NumPy_matrices", "Tunnel_matrix.npy"), data.Tunnel_matrix.IDs)
np.save(os.path.join("NumPy_matrices", "Energy_matrix.npy"), data.Energy_matrix)

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        return super().default(obj)

with open("EnGOAT_data.json", "w") as f:
    json.dump(EnGOAT_data, f, cls=NumpyEncoder, indent=4)



"""
WARNING: basin==cluster in the code!
"""

