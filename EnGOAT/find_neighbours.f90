Subroutine find_neighbours(level, point, neighbours, cross_matrix_i, cross_matrix_j, cross_matrix_k)
    implicit none
    integer(8), intent(in):: level
    integer(8), intent(in):: point(9)                                                                                  !Current point
    integer(8), intent(inout) :: cross_matrix_i(:,:,:), cross_matrix_j(:,:,:), cross_matrix_k(:,:,:)                   !Matrix to keep track of crossing information in each dimension
    integer(8), intent(out):: neighbours(6, 9)                                                                         !PBC neighbours of the current point
    integer(8):: i, j, k                                                                                               !Current point coordinates
    integer(8):: ip, im, jp, jm, kp, km                                                                                !Coordinates of neighbouring points. 9 arguments per point: 1:3 = point coordinates, 4:6 = cross vector information, 7 = boundary information, 8=level, 9=TS info
    integer(8):: cross(6)                                                                                              !Crossing information: 1 if crossing positive boundary, -1 if crossing negative boundary, 0 otherwise
    integer(8):: grid_size(3)                                                                                          !Size (in number of points) of the level matrix
            
    grid_size(1) = size(cross_matrix_i, 1)             
    grid_size(2) = size(cross_matrix_i, 2)             
    grid_size(3) = size(cross_matrix_i, 3)             
                
    i=point(1)+1_8                                                                                                     !+1 is due to the difference in fortran and python labeling
    j=point(2)+1_8                                                                                                     !(python starts at 0, fortran at 1)
    k=point(3)+1_8

    call PBC3D(i, j, k, grid_size, ip, im, jp, jm, kp, km, cross)

    neighbours(1,:) = [ip-1_8, j-1_8, k-1_8, &
        cross_matrix_i(i, j, k) + cross(1), cross_matrix_j(i, j, k), cross_matrix_k(i, j, k), 1_8, level, 0_8]         !Boundary info = 1, since new points are on the edge
    neighbours(2,:) = [im-1_8, j-1_8, k-1_8, &
        cross_matrix_i(i, j, k) + cross(2), cross_matrix_j(i, j, k), cross_matrix_k(i, j, k), 1_8, level, 0_8]
    neighbours(3,:) = [i-1_8, jp-1_8, k-1_8, &
        cross_matrix_i(i, j, k), cross_matrix_j(i, j, k) + cross(3), cross_matrix_k(i, j, k), 1_8, level, 0_8]
    neighbours(4,:) = [i-1_8, jm-1_8, k-1_8, &
        cross_matrix_i(i, j, k), cross_matrix_j(i, j, k) + cross(4), cross_matrix_k(i, j, k), 1_8, level, 0_8]
    neighbours(5,:) = [i-1_8, j-1_8, kp-1_8, &
        cross_matrix_i(i, j, k), cross_matrix_j(i, j, k), cross_matrix_k(i, j, k) + cross(5), 1_8, level, 0_8]
    neighbours(6,:) = [i-1_8, j-1_8, km-1_8, &
        cross_matrix_i(i, j, k), cross_matrix_j(i, j, k), cross_matrix_k(i, j, k) + cross(6), 1_8, level, 0_8]

end subroutine find_neighbours


Subroutine PBC3D(i, j, k, grid_size, ip, im, jp, jm, kp, km, cross)
    implicit none
    integer(8), intent(in):: i, j, k                                                                                    !Current coordinates
    integer(8), intent(in):: grid_size(3)                                                                               !Grid size in each dimension
    integer(8), intent(out):: ip, im, jp, jm, kp, km                                                                    !Periodic boundary condition neighbor point coordinates
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
!python3 -m numpy.f2py -c -m find_neighbours find_neighbours.f90