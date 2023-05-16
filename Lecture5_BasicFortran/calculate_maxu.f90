    subroutine calculate_maxu
    
    include 'com.txt'
    
    maxu = 0
    
    do i = 0,Nx
        do j = 0,Ny
            
            if (abs(uh(i,j)) > maxu) then
                maxu = abs(uh(i,j))
            end if
            
        end do
    end do
    
    end subroutine calculate_maxu