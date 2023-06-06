    subroutine save_solution
    
    include 'com.txt'
    
    tRK = tend
    call set_bc
    
    open(unit = 1,file = 'Q1.txt')
    open(unit = 2,file = 'Q2.txt')
    open(unit = 3,file = 'Q3.txt')
    open(unit = 4,file = 'Q4.txt')
    
    
    do d = 1,dimPk
        do j = 1,Ny
            do i = 1,Nx
                write(1,*) uh(i,j,d,1)
                write(2,*) uh(i,j,d,2)
                write(3,*) uh(i,j,d,3)
                write(4,*) uh(i,j,d,4)
            end do
        end do
    end do
    
    close(1)
    close(2)
    close(3)
    close(4)
    
    end subroutine save_solution