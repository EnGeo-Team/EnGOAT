Subroutine get_process_cross_vector(start_point, end_point, C1_C2_TS, minIDmatrix_clusters, process_cross_vector)
    implicit none
    integer(8), intent(in):: start_point(3), end_point(3)                                                              !Start and end point of the process cross vector search
    integer(8), intent(in):: C1_C2_TS(3)                                                                               !IDs of cluster1, cluster2, and transition state
    integer(8), intent(in) :: minIDmatrix_Clusters(:,:,:)                                                              !Matrices containing 1) the ID of the cluster a given point is asigned to (0 if the point belongs to no clusters and -1 if a point is a transition state) and 2) the ID of transition state a given point is asigned to (0 if the point is not a transition state)
    integer(8), intent(out):: process_cross_vector(3)                                                                  !PBC neighbours of the current point
    integer(8), allocatable:: cross_matrix_i(:,:,:), cross_matrix_j(:,:,:), cross_matrix_k(:,:,:)                      !Matrices holding crossing information of visited points
    integer(8), allocatable:: check_matrix(:,:,:)                                                                      !Matrix holding information on whether a point has already been checked (added to the point list) or not
    integer(8):: i_start, j_start, k_start, i_end, j_end, k_end                                                        !Start and end point coordinates
    integer(8):: i, j, k                                                                                               !Current point coordinates
    integer(8):: ip, im, jp, jm, kp, km                                                                                !Coordinates of neighbouring points.
    integer(8):: neighbours(6, 6)                                                                                      !Holds neighbour coordinates (1-3) and crossing information of neighbour points (4-6) relative to starting point
    integer(8), allocatable:: point_list(:,:), temp(:,:)                                                               !A list holding all visited points
    integer(8):: point_index, N_points                                                                                 !Index of the current point being checked. Total number of points visited (size of the point_list matrix)
    integer(8):: C1, C2, TS                                                                                            !IDs of cluster 1, cluster 2, and the transition state.
    integer(8):: cross(6)                                                                                              !Crossing information: 1 if crossing positive boundary, -1 if crossing negative boundary, 0 otherwise
    integer(8):: grid_size(3)                                                                                          !Size (in number of points) of the matrices
    integer(8):: neighbour, neighbour_ID, included                                                                     !Index for cycling through neighbours. minIDmatrix_clusters entry of the neighbour point. A flag indicating whether a neighbour point is included in C1 C2 TS group (1) or not (0).
            
    grid_size(1) = size(minIDmatrix_clusters, 1)             
    grid_size(2) = size(minIDmatrix_clusters, 2)             
    grid_size(3) = size(minIDmatrix_clusters, 3)
    
    allocate(cross_matrix_i(grid_size(1), grid_size(2), grid_size(3)))                                                 !Allocate the cross matrices and set their value to 0
    allocate(cross_matrix_j(grid_size(1), grid_size(2), grid_size(3)))
    allocate(cross_matrix_k(grid_size(1), grid_size(2), grid_size(3)))
    cross_matrix_i=0_8
    cross_matrix_j=0_8
    cross_matrix_k=0_8

    allocate(check_matrix(grid_size(1), grid_size(2), grid_size(3)))                                                   !Allocate the check matrix and set its value to 0
    check_matrix=0_8

    C1 = C1_C2_TS(1)                                                                                                   !Set the starting, ending cluster, and transition state IDs
    C2 = C1_C2_TS(2)
    TS = C1_C2_TS(3)

    i_start=start_point(1)+1_8                                                                                         !Set start and end point coordinates
    j_start=start_point(2)+1_8                                                                                         !+1 is due to the difference in fortran and python labeling
    k_start=start_point(3)+1_8                                                                                         !(python starts at 0, fortran at 1)  
    i_end=end_point(1)+1_8                                                                                             
    j_end=end_point(2)+1_8                                                                                             
    k_end=end_point(3)+1_8                                                                                                                                                                  

    N_points=1
    allocate(point_list(N_points, 6))
    point_list(N_points, :) = [i_start, j_start, k_start, 0_8, 0_8, 0_8]                                               !Add the starting point to the point list
    point_index=1_8                                                                                                    !Set the point_index to the first point in the point list (starting point)
    check_matrix(i_start, j_start, k_start)=1_8                                                                        !Set the starting point to "checked"

    do while (point_index<=N_points)

        i=point_list(point_index, 1)
        j=point_list(point_index, 2)
        k=point_list(point_index, 3)

        call PBC3D(i, j, k, grid_size, ip, im, jp, jm, kp, km, cross)

        neighbours(1, :) = [ip, j, k, &
            cross_matrix_i(i, j, k) + cross(1), cross_matrix_j(i, j, k), cross_matrix_k(i, j, k)]
        neighbours(2, :) = [im, j, k, &
            cross_matrix_i(i, j, k) + cross(2), cross_matrix_j(i, j, k), cross_matrix_k(i, j, k)]
        neighbours(3, :) = [i, jp, k, &
            cross_matrix_i(i, j, k), cross_matrix_j(i, j, k) + cross(3), cross_matrix_k(i, j, k)]
        neighbours(4, :) = [i, jm, k, &
            cross_matrix_i(i, j, k), cross_matrix_j(i, j, k) + cross(4), cross_matrix_k(i, j, k)]
        neighbours(5, :) = [i, j, kp, &
            cross_matrix_i(i, j, k), cross_matrix_j(i, j, k), cross_matrix_k(i, j, k) + cross(5)]
        neighbours(6, :) = [i, j, km, &
            cross_matrix_i(i, j, k), cross_matrix_j(i, j, k), cross_matrix_k(i, j, k) + cross(6)]

        do neighbour = 1, 6

            included=0_8                                                                                              !Check whether the neighbour is in the C1 C2 TS group. if the point in minIDmatrix is labeled as a transition state, only include the points of the TS transition state!
            neighbour_ID=minIDmatrix_clusters(neighbours(neighbour, 1), neighbours(neighbour, 2), neighbours(neighbour, 3))
            if ((neighbour_ID==C1).or.(neighbour_ID==C2)) then
                included=1_8
            end if

            if (included==1_8) then
                if (check_matrix(neighbours(neighbour, 1), neighbours(neighbour, 2), neighbours(neighbour, 3))==0) then     !Check whether the neighbour point has already been checked (added to the point list)
                
                    cross_matrix_i(neighbours(neighbour, 1), neighbours(neighbour, 2), neighbours(neighbour, 3))=&
                        neighbours(neighbour, 4)                                                                       !Save the crossing information of the neighbour point
                    cross_matrix_j(neighbours(neighbour, 1), neighbours(neighbour, 2), neighbours(neighbour, 3))=&
                        neighbours(neighbour, 5)
                    cross_matrix_k(neighbours(neighbour, 1), neighbours(neighbour, 2), neighbours(neighbour, 3))=&
                        neighbours(neighbour, 6)

                    allocate(temp(N_points, 6))                                                                        !store all old points in the temporary matrix
                    temp=point_list
                    deallocate(point_list)

                    allocate(point_list(N_points+1, 6))                                                                !add the new point (along with all the old ones) to the Cluster_points matrix
                    point_list(1:N_points, :)=temp(1:N_points, :)
                    point_list(N_points+1, :)=neighbours(neighbour, :)
                    deallocate(temp)
                    N_points=N_points+1
                    check_matrix(neighbours(neighbour, 1), neighbours(neighbour, 2), neighbours(neighbour, 3))=1_8     !Set the added point to "checked"

                    if ((neighbours(neighbour, 1)==i_end).and.& 
                        (neighbours(neighbour, 2)==j_end).and.&
                        (neighbours(neighbour, 3)==k_end)) then                                                        !Check if the end point was reached

                        process_cross_vector(1)=neighbours(neighbour, 4)
                        process_cross_vector(2)=neighbours(neighbour, 5)
                        process_cross_vector(3)=neighbours(neighbour, 6)
                        return
                    end if

                end if
            end if

        end do

        point_index=point_index+1_8

        if (point_index>N_points) then
            write(*,*) "Error! Endpoint was not found!"
            process_cross_vector=[100_8, 100_8, 100_8]            
            exit
        end if

    end do

end subroutine get_process_cross_vector


Subroutine PBC3D(i, j, k, grid_size, ip, im, jp, jm, kp, km, cross)
    implicit none
    integer(8), intent(in):: i, j, k                                                                                    !Current coordinates
    integer(8), intent(in):: grid_size(3)                                                                               !Grid size in each dimension
    integer(8), intent(out):: ip, im, jp, jm, kp, km                                                                    !Periodic boundary condition neighbour point coordinates
    integer(8), intent(out):: cross(6)                                                                                  !Crossing information: 1 if crossing positive boundary, -1 if crossing negative boundary, 0 otherwise
    cross=0                                                                                                             !cross=(positive x crossing, negative x crossing, positive y crossing, negative y crossing, positive z crossing, negative z crossing)
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
!python3 -m numpy.f2py -c -m get_process_cross_vector get_process_cross_vector.f90