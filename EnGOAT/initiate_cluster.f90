subroutine find_new_clusters(level, N_clusters, i_0, j_0, k_0, &
                            level_matrix, Cluster_matrix_Levels, Cluster_matrix_IDs, &
                            cross_matrix_i, cross_matrix_j, cross_matrix_k, &
                            Energy_matrix, E_min)

    implicit none
    integer(8), intent(in) :: level                                                                     !Current level 
    integer(8), intent(in) :: i_0, j_0, k_0                                                             !Coordinates of the unexplored point
    integer(8), intent(in) :: level_matrix(:,:,:)                                                       !Level matrix
    real(8), intent(in) :: Energy_matrix(:,:,:)                                                         !Energy matrix (real numbers not levels)
    integer(8), intent(inout) :: Cluster_matrix_Levels(:,:,:), Cluster_matrix_IDs(:,:,:)                !Matrix to keep track of the level at which each point was first explored and the cluster ID of each point                      
    integer(8), intent(inout) :: cross_matrix_i(:,:,:), cross_matrix_j(:,:,:), cross_matrix_k(:,:,:)    !Matrix to keep track of crossing information in each dimension
    integer(8), intent(inout) :: N_clusters(1)                                                          !Total number of clusters (serves as cluster ID). A 1D array due to python compatibility issues
    real(8), intent(out) :: E_min                                                                       !Energy value of the point with the lowest energy in the cluster
 
    integer(8) :: grid_size(3)                                                                          !Grid size in each dimension
    integer(8) :: i, j, k                                                                               !coordinates of the current point                               
    integer(8) :: ip, im, jp, jm, kp, km, cross(6)                                                      !Periodic boundary condition neighbor point coordinates and crossing information
    integer(8) :: cross_i, cross_j, cross_k                                                             !Temporary scalars just for writing the crossing information into matrices and points
    integer(8) :: neighbor_ID(6)                                                                        !minID_Levels of the 6 neighboring points

    integer(8) :: N_points                                                                              !Number of points belonging to the current cluster
    integer(8) :: point_index                                                                           !Index of the point in the cluster being checked
    integer(8) :: boundary                                                                              !Flag to indicate if the current cluster touches the boundary                     

    integer(8), allocatable :: Cluster_points(:,:)                                                      !List of points in the current cluster. 9 arguments per point: 1:3 = point coordinates, 4:6 = cross vector information, 7 = boundary information, 8=level 9=TS information
    integer(8), allocatable :: temp(:,:)                                                                !Temporary matrix to store data
    integer(8) :: point                                                                                 !Just a counter for writing
    real(8), allocatable :: cluster_energies(:)                                                         !A matrix to store energies of the cluster points in order to find E_min

    grid_size(1) = size(level_matrix, 1)                                    
    grid_size(2) = size(level_matrix, 2)                                    
    grid_size(3) = size(level_matrix, 3)                                    

    i = i_0+1_8                                                                                         !Start from the unexplored point. 
    j = j_0+1_8                                                                                         !+1 is due to the difference in fortran and python labeling
    k = k_0+1_8                                                                                         !(python starts at 0, fortran at 1)

    call PBC3D(i, j, k, grid_size, ip, im, jp, jm, kp, km, cross)                                   

    neighbor_ID(1) = Cluster_matrix_IDs(ip, j, k)                                                    
    neighbor_ID(2) = Cluster_matrix_IDs(im, j, k)                                                    
    neighbor_ID(3) = Cluster_matrix_IDs(i, jp, k)                                                    
    neighbor_ID(4) = Cluster_matrix_IDs(i, jm, k)                                                    
    neighbor_ID(5) = Cluster_matrix_IDs(i, j, kp)                                                    
    neighbor_ID(6) = Cluster_matrix_IDs(i, j, km)                                                    

    if (sum(neighbor_ID)==0) then                                                                       !If all the neighbour points are yet unchecked
        N_points=1                                                                                      !Initialize new cluster                 
        N_clusters(1)=N_clusters(1)+1
        Cluster_matrix_Levels(i,j,k)=level
        Cluster_matrix_IDs(i,j,k)=N_clusters(1)
        point_index=1
        boundary=0         

        allocate(cluster_points(1,9))                                                                   !initialize list of points for this cluster and add the first point
        cluster_points(1,:)=[i, j, k, 0_8, 0_8, 0_8, boundary, level, 0_8]                              !1:3 = point coordinates, 4:6 = cross vector information, 7 = boundary information, 8=level, 9=TS info
        do while (point_index<=N_points)                                                                !Perform a neighbourhood search for points at the same level
            if (Cluster_matrix_Levels(ip, j, k)==0) then                  
                if (level_matrix(ip, j, k)==level) then                                                 !If yet unidentified neighbour point is found at a current level, add it to the cluster
                    Cluster_matrix_Levels(ip, j, k)=level                 
                    Cluster_matrix_IDs(ip, j, k)=N_clusters(1)                       

                    cross_i=cross_matrix_i(i, j, k)+cross(1)                    
                    cross_j=cross_matrix_j(i, j, k)                 
                    cross_k=cross_matrix_k(i, j, k)                                     
                    cross_matrix_i(ip, j, k)=cross_i                                                    !Cross_matrix entry of the neighbouring point is crossing info of the examined point + crossing vector
                    cross_matrix_j(ip, j, k)=cross_j                    
                    cross_matrix_k(ip, j, k)=cross_k                    

                    allocate(temp(N_points, 9))                                                         !store all old points in the temporary matrix
                    temp=cluster_points                 
                    deallocate(cluster_points)                  

                    allocate(Cluster_points(N_points+1, 9))                                             !add the new point (along with all the old ones) to the Cluster_points matrix
                    cluster_points(1:N_points, :)=temp(1:N_points, :)
                    cluster_points(N_points+1, :)=[ip, j, k, cross_i, cross_j, cross_k, 0_8, level, 0_8]
                    deallocate(temp)
                    N_points=N_points+1

                else
                    boundary=1                                                                          !If the point doesn't have neighbours on ALL sides, it is the boundary point
                end if                  
            end if                  
            if (Cluster_matrix_Levels(im, j, k)==0) then                                                !Repeat everything above for all 6 neighbours.
                if (level_matrix(im, j, k)==level) then                 
                    Cluster_matrix_Levels(im, j, k)=level
                    Cluster_matrix_IDs(im, j, k)=N_clusters(1)

                    cross_i=cross_matrix_i(i, j, k)+cross(2)
                    cross_j=cross_matrix_j(i, j, k)
                    cross_k=cross_matrix_k(i, j, k)                    
                    cross_matrix_i(im, j, k)=cross_i                                
                    cross_matrix_j(im, j, k)=cross_j
                    cross_matrix_k(im, j, k)=cross_k

                    allocate(temp(N_points, 9))                                    
                    temp=cluster_points
                    deallocate(cluster_points)
            
                    allocate(Cluster_points(N_points+1, 9))                        
                    cluster_points(1:N_points, :)=temp(1:N_points, :)
                    cluster_points(N_points+1, :)=[im, j, k, cross_i, cross_j, cross_k, 0_8, level, 0_8]
                    deallocate(temp)
                    N_points=N_points+1
                else
                    boundary=1
                end if
            end if
            if (Cluster_matrix_Levels(i, jp, k)==0) then
                if (level_matrix(i, jp, k)==level) then
                    Cluster_matrix_Levels(i, jp, k)=level
                    Cluster_matrix_IDs(i, jp, k)=N_clusters(1)

                    cross_i=cross_matrix_i(i, j, k)
                    cross_j=cross_matrix_j(i, j, k)+cross(3)
                    cross_k=cross_matrix_k(i, j, k)                    
                    cross_matrix_i(i, jp, k)=cross_i                              
                    cross_matrix_j(i, jp, k)=cross_j
                    cross_matrix_k(i, jp, k)=cross_k

                    allocate(temp(N_points, 9))                                   
                    temp=cluster_points
                    deallocate(cluster_points)
            
                    allocate(Cluster_points(N_points+1, 9))                        
                    cluster_points(1:N_points, :)=temp(1:N_points, :)
                    cluster_points(N_points+1, :)=[i, jp, k, cross_i, cross_j, cross_k, 0_8, level, 0_8]
                    deallocate(temp)
                    N_points=N_points+1
                else
                    boundary=1
                end if
            end if
            if (Cluster_matrix_Levels(i, jm, k)==0) then
                if (level_matrix(i, jm, k)==level) then 
                    Cluster_matrix_Levels(i, jm, k)=level
                    Cluster_matrix_IDs(i, jm, k)=N_clusters(1)

                    cross_i=cross_matrix_i(i, j, k)
                    cross_j=cross_matrix_j(i, j, k)+cross(4)
                    cross_k=cross_matrix_k(i, j, k)                    
                    cross_matrix_i(i, jm, k)=cross_i                               
                    cross_matrix_j(i, jm, k)=cross_j
                    cross_matrix_k(i, jm, k)=cross_k                    

                    allocate(temp(N_points, 9))                                 
                    temp=cluster_points
                    deallocate(cluster_points)
            
                    allocate(Cluster_points(N_points+1, 9))                      
                    cluster_points(1:N_points, :)=temp(1:N_points, :)
                    cluster_points(N_points+1, :)=[i, jm, k, cross_i, cross_j, cross_k, 0_8, level, 0_8]
                    deallocate(temp)
                    N_points=N_points+1
                else
                    boundary=1
                end if
            end if
            if (Cluster_matrix_Levels(i, j, kp)==0) then
                if (level_matrix(i, j, kp)==level) then
                    Cluster_matrix_Levels(i, j, kp)=level
                    Cluster_matrix_IDs(i, j, kp)=N_clusters(1)

                    cross_i=cross_matrix_i(i, j, k)
                    cross_j=cross_matrix_j(i, j, k)
                    cross_k=cross_matrix_k(i, j, k)+cross(5)                 
                    cross_matrix_i(i, j, kp)=cross_i                           
                    cross_matrix_j(i, j, kp)=cross_j
                    cross_matrix_k(i, j, kp)=cross_k

                    allocate(temp(N_points, 9))                                  
                    temp=cluster_points
                    deallocate(cluster_points)
            
                    allocate(Cluster_points(N_points+1, 9))                       
                    cluster_points(1:N_points, :)=temp(1:N_points, :)
                    cluster_points(N_points+1, :)=[i, j, kp, cross_i, cross_j, cross_k, 0_8, level, 0_8]
                    deallocate(temp)
                    N_points=N_points+1
                else
                    boundary=1
                end if
            end if
            if (Cluster_matrix_Levels(i, j, km)==0) then
                if (level_matrix(i, j, km)==level) then                    
                    Cluster_matrix_Levels(i, j, km)=level
                    Cluster_matrix_IDs(i, j, km)=N_clusters(1)

                    cross_i=cross_matrix_i(i, j, k)
                    cross_j=cross_matrix_j(i, j, k)
                    cross_k=cross_matrix_k(i, j, k)+cross(6)                 
                    cross_matrix_i(i, j, km)=cross_i                              
                    cross_matrix_j(i, j, km)=cross_j
                    cross_matrix_k(i, j, km)=cross_k

                    allocate(temp(N_points, 9))                                    
                    temp=cluster_points
                    deallocate(cluster_points)
            
                    allocate(Cluster_points(N_points+1, 9))                      
                    cluster_points(1:N_points, :)=temp(1:N_points, :)
                    cluster_points(N_points+1, :)=[i, j, km, cross_i, cross_j, cross_k, 0_8, level, 0_8]
                    deallocate(temp)
                    N_points=N_points+1
                else
                    boundary=1
                end if
            end if
            cluster_points(point_index,7)=boundary                                                      !Save the boundary information for the examined point
            if (point_index<N_points) then                                                              !Move to one of the neighbour points that were added to the cluster and examine it, repeat until all the points in the cluster are examined.
                point_index=point_index+1                   
                i=cluster_points(point_index,1)                 
                j=cluster_points(point_index,2)                 
                k=cluster_points(point_index,3)                 
                call PBC3D(i, j, k, grid_size, ip, im, jp, jm, kp, km, cross)                   
            else                    
                exit                    
            end if                  
                    
        end do                  
                    
        allocate(cluster_energies(N_points))                                                            !Find the energies of all points in the cluster
        open(10, file="temp1.dat")
        open(11, file="temp2.dat")
        do point = 1, N_points
            write(10, '(9I4)') cluster_points(point, :)
            cluster_energies(point)=Energy_matrix(cluster_points(point, 1), cluster_points(point, 2), cluster_points(point, 3))
            write(11, "(F15.8)") cluster_energies(point)
        end do
        close(10)
        close(11)
        E_min=minval(cluster_energies(:))                                                               !Find the lowest energy of the cluster
    end if

end subroutine find_new_clusters

Subroutine PBC3D(i, j, k, grid_size, ip, im, jp, jm, kp, km, cross)
    implicit none
    integer(8), intent(in):: i, j, k                                                                    !Current coordinates
    integer(8), intent(in):: grid_size(3)                                                               !Grid size in each dimension
    integer(8), intent(out):: ip, im, jp, jm, kp, km                                                    !Periodic boundary condition neighbor point coordinates
    integer(8), intent(out):: cross(6)                                                                  !Crossing information: 1 if crossing positive boundary, -1 if crossing negative boundary, 0 otherwise
    cross=0                                                                                             !cross=(positive x crossing, negative x crossing, positive y crossing, negative y crossing, positive z crossing, negative z crossing)
    if (i==grid_size(1)) then
        im=i-1
        ip=1
        cross(1)=1
    elseif (i==1) then
        im=grid_size(1)
        ip=i+1
        cross(2)=-1
    else
        im=i-1
        ip=i+1
    end if
    if (j==grid_size(2)) then
        jm=j-1
        jp=1
        cross(3)=1
    elseif (j==1) then
        jm=grid_size(2)
        jp=j+1
        cross(4)=-1
    else
        jm=j-1
        jp=j+1
    end if
    if (k==grid_size(3)) then
        km=k-1
        kp=1
        cross(5)=1
    elseif (k==1) then
        km=grid_size(3)
        kp=k+1
        cross(6)=-1
    else
        km=k-1
        kp=k+1
    end if
end subroutine PBC3D

!Fortran routine has to be compiled before running the code! Write the following line in the terminal:
!
!python3 -m numpy.f2py -c -m initiate_cluster initiate_cluster.f90