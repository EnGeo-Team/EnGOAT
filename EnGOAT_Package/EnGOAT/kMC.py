import numpy as np
import logging
import random
import time
import os
import data                                                                                                 #Python file containing all the data-containing matrices used in the simulation
logger = logging.getLogger(__name__)

class basis_site:                                                                                           #Define basis sites as objects
    def __init__(self, ID, center, process_list, k_sum):
        self.ID = ID                                                                                        #Cluster ID of the cluster containing the site
        self.center = center                                                                                #Coordinates of the cluster's center
        self.process_list = process_list                                                                    #List of processes beginning in a given cluster
        self.k_sum = k_sum                                                                                  #Sum of all process rates beginning in a given cluster

def kMC(tunnel_system, N_kMC_runs, N_kMC_steps, T, tunnel_idx):                                             #Performs kMC simulation of a given tunnel system at temperature T
    timer_total=time.process_time()
    msd_all_runs = []                                                                                       #List holding msd plots of every kMC run (in a b c space)
    msd_xyz_all_runs = []                                                                                   #List holding msd plots of every kMC run (in x y z space)
    msd_3D_all_runs = []                                                                                    #List holding msd plots of every kMC run
    D = np.zeros((N_kMC_runs, 3))                                                                           #Diffusion coefficient in the a b c directions
    D_xyz = np.zeros((N_kMC_runs, 3))                                                                       #Diffusion coefficient in the x y z directions
    D_3D = np.zeros((N_kMC_runs))                                                                           #Isotropic diffusion coefficient

    msd_steps = []                                                                                          #Create an array holding log order steps for computing the msd
    for log_order in range(int(np.log10(N_kMC_steps / 10)) + 1):
        msd_steps.extend((10 ** log_order) * np.arange(1, 10))
    msd_steps = np.array(msd_steps, dtype=int)

    basis_sites = []                                                                                        #Create a list of basis sites (objects)
    for cluster_ID in tunnel_system.cluster_family:
        process_list = []
        k_sum=0
        for process in tunnel_system.process_list:
            if process.start_cluster==cluster_ID:
                process_list.append(process)
                k_sum += process.k
        if len(process_list):    
            center = process_list[0].start_point                                                            #Exclude basins with no transition processes (happens due to 1 point cluster deletion)
            basis_sites.append(basis_site(cluster_ID, center, process_list, k_sum))                 
    
    #Define the transition matrix from a b c -> x y z space
    a_vec = data.grid["vectors"][0, :]
    b_vec = data.grid["vectors"][1, :]
    c_vec = data.grid["vectors"][2, :]
    a_vec = a_vec/np.linalg.norm(a_vec)
    b_vec = b_vec/np.linalg.norm(b_vec)
    c_vec = c_vec/np.linalg.norm(c_vec)
    cell_xyz = np.vstack([a_vec*data.grid["shape"][0]*data.grid["voxel_size"][0], b_vec*data.grid["shape"][1]*data.grid["voxel_size"][1], c_vec*data.grid["shape"][2]*data.grid["voxel_size"][2]])                                             
    cell_size_xyz = np.linalg.norm(cell_xyz, axis=0)

    for run in range(N_kMC_runs):
        current_site=basis_sites[0]                                                                         #Chose an arbitrary starting point (first basis site in the list)
        timer_run=time.process_time()
        t=0
        trajectory = np.zeros((N_kMC_steps, 4))                                                             #Initiate an array storing the trajectory at every step of a given run [time, position a, position b, position c]
        for step in range(1, N_kMC_steps):
            pick_process = current_site.k_sum * random.random()
            process_sum = 0
            for process in current_site.process_list:
                process_sum += process.k
                if process_sum >= pick_process:
                    process_event = process
                    break
            t += np.log(1/random.random()) / current_site.k_sum                                             #Move time forward
            trajectory[step, 0] = t
            trajectory[step, 1:] = trajectory[step-1, 1:] + process_event.distance
            for site in basis_sites:
                if site.ID==process_event.end_cluster:
                    current_site=site
                    break
        #Transform the trajectory to the cartesian space
        trajectory_xyz = np.zeros_like(trajectory)
        trajectory_xyz[:, 0] = trajectory[:, 0]
        trajectory_xyz[:, 1] = trajectory[:, 1] * a_vec[0] + trajectory[:, 2] * b_vec[0] + trajectory[:, 3] * c_vec[0]
        trajectory_xyz[:, 2] = trajectory[:, 1] * a_vec[1] + trajectory[:, 2] * b_vec[1] + trajectory[:, 3] * c_vec[1]
        trajectory_xyz[:, 3] = trajectory[:, 1] * a_vec[2] + trajectory[:, 2] * b_vec[2] + trajectory[:, 3] * c_vec[2]
                                                       

        msd = np.zeros((len(msd_steps), 4))                                                                 # columns: [Δt, ⟨Δa²⟩, ⟨Δb²⟩, ⟨Δc²⟩]
        msd_xyz = np.zeros((len(msd_steps), 4))                                                             # columns: [Δt, ⟨Δx²⟩, ⟨Δy²⟩, ⟨Δz²⟩]
        msd_3D = np.zeros((len(msd_steps), 2))                                                              # columns: [Δt, ⟨Δr²⟩]

        for k, j in enumerate(msd_steps):
            diff = trajectory[j:, :] - trajectory[:-j, :]
            msd[k, 0] = np.mean(diff[:, 0])                                                                 # Δt average
            msd[k, 1] = np.mean(diff[:, 1]**2)                                                              # ⟨Δa²⟩
            msd[k, 2] = np.mean(diff[:, 2]**2)                                                              # ⟨Δb²⟩
            msd[k, 3] = np.mean(diff[:, 3]**2)                                                              # ⟨Δc²⟩

            diff_xyz = trajectory_xyz[j:, :] - trajectory_xyz[:-j, :]
            msd_xyz[k, 0] = np.mean(diff_xyz[:, 0])                                                         # Δt average
            msd_xyz[k, 1] = np.mean(diff_xyz[:, 1]**2)                                                      # ⟨Δx²⟩
            msd_xyz[k, 2] = np.mean(diff_xyz[:, 2]**2)                                                      # ⟨Δy²⟩
            msd_xyz[k, 3] = np.mean(diff_xyz[:, 3]**2)                                                      # ⟨Δz²⟩

            msd_3D[k, 0] = np.mean(diff_xyz[:, 0])
            msd_3D[k, 1] = np.mean(diff_xyz[:, 1]**2+diff_xyz[:, 2]**2+diff_xyz[:, 3]**2)
        
        msd_all_runs.append(msd)
        msd_xyz_all_runs.append(msd_xyz)
        msd_3D_all_runs.append(msd_3D)

        #a b c direction:
        for direction in range(3):                                                                           

            cond1 = msd[:, direction + 1] > (data.grid["voxel_size"][direction]*data.grid["shape"][direction])**2           #Find start and end indices of the linear range of msd (between 1 and 4 unit cells)               
            cond2 = msd[:, direction + 1] > 4 * (data.grid["voxel_size"][direction]*data.grid["shape"][direction])**2
            idx_start = np.where(cond1)[0][0] if np.any(cond1) else None
            idx_end   = np.where(cond2)[0][0] if np.any(cond2) else None
            if idx_end == None:
                D[run, direction] = 0
            else:
                x_values = msd[idx_start:idx_end, 0]
                y_values = msd[idx_start:idx_end, direction+1]
                slope, intercept = np.polyfit(x_values, y_values, 1)
                D[run, direction] = 0.5 * slope * 10**(-16)

                """Save the MSD files for a b c directions"""
                directions = ["a2", "b2", "c2"]

        #x y z directions
        for direction in range(3):                                                                                                                               

            cond1 = msd_xyz[:, direction + 1] > (cell_size_xyz[direction])**2                               #Find start and end indices of the linear range of msd (between 1 and 4 unit cells)               
            cond2 = msd_xyz[:, direction + 1] > 4 * (cell_size_xyz[direction])**2
            idx_start = np.where(cond1)[0][0] if np.any(cond1) else None
            idx_end   = np.where(cond2)[0][0] if np.any(cond2) else None
            if idx_end == None:
                D_xyz[run, direction] = 0
            else:
                x_values = msd_xyz[idx_start:idx_end, 0]
                y_values = msd_xyz[idx_start:idx_end, direction+1]
                slope, intercept = np.polyfit(x_values, y_values, 1)
                D_xyz[run, direction] = 0.5 * slope * 10**(-16)

                """Save the MSD files for x y z directions"""
                directions = ["x2", "y2", "z2"]
        
        #Isotropic diffusion coefficient:
        cond1 = msd_3D[:, 1] > (cell_size_xyz[0]**2 + cell_size_xyz[1]**2 + cell_size_xyz[2]**2)            #Find start and end indices of the linear range of msd (between 1 and 4 unit cells)               
        cond2 = msd_3D[:, 1] > 4 * (cell_size_xyz[0]**2 + cell_size_xyz[1]**2 + cell_size_xyz[2]**2)
        idx_start = np.where(cond1)[0][0] if np.any(cond1) else None
        idx_end   = np.where(cond2)[0][0] if np.any(cond2) else None

        if idx_end == None:
            D_3D[run] = 0
        else:
            x_values = msd_3D[idx_start:idx_end, 0]
            y_values = msd_3D[idx_start:idx_end, 1]
            slope, intercept = np.polyfit(x_values, y_values, 1)
            D_3D[run] = 1.0/6.0 * slope * 10**(-16)

        logger.info(
            f"{run+1:8d} "
            f"{D[run, 0]:>15.3e} "
            f"{D[run, 1]:>15.3e} "
            f"{D[run, 2]:>15.3e} "
            f"{time.process_time()-timer_run:>10.2f} "
            f"{time.process_time()-timer_total:>10.2f} ")

    D_average = np.array([
    [np.mean(D[:, 0]), np.std(D[:, 0])],
    [np.mean(D[:, 1]), np.std(D[:, 1])],
    [np.mean(D[:, 2]), np.std(D[:, 2])]])

    D_xyz_average = np.array([
    [np.mean(D_xyz[:, 0]), np.std(D_xyz[:, 0])],
    [np.mean(D_xyz[:, 1]), np.std(D_xyz[:, 1])],
    [np.mean(D_xyz[:, 2]), np.std(D_xyz[:, 2])]])

    D_3D_average = np.array([np.mean(D_3D[:]), np.std(D_3D[:])])

    kMC_results = {}

    # a, b, c directions
    for i, direction in enumerate(("a", "b", "c")):
        kMC_results[direction] = {
            "D": float(D_average[i, 0]),
            "sd": float(D_average[i, 1]),
            "MSD": {
                f"run_{run+1}": {
                    "time": msd_all_runs[run][:, 0].tolist(),
                    "MSD": msd_all_runs[run][:, i+1].tolist(),
                }
                for run in range(N_kMC_runs)
            }
        }
    
    # x, y, z directions
    for i, direction in enumerate(("x", "y", "z")):
        kMC_results[direction] = {
            "D": float(D_xyz_average[i, 0]),
            "sd": float(D_xyz_average[i, 1]),
            "MSD": {
                f"run_{run+1}": {
                    "time": msd_xyz_all_runs[run][:, 0].tolist(),
                    "MSD": msd_xyz_all_runs[run][:, i+1].tolist(),
                }
                for run in range(N_kMC_runs)
            }
        }
    
    # Isotropic diffusion
    kMC_results["3D"] = {
        "D": float(D_3D_average[0]),
        "sd": float(D_3D_average[1]),
        "MSD": {
            f"run_{run+1}": {
                "time": msd_3D_all_runs[run][:, 0].tolist(),
                "MSD": msd_3D_all_runs[run][:, 1].tolist(),
            }
            for run in range(N_kMC_runs)
        }
    }

    return kMC_results