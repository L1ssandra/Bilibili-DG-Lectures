    subroutine save_solution
    
    include 'com.txt'
    
    open(unit = 1,file = 'Q1.txt')
    open(unit = 2,file = 'Q2.txt')
    open(unit = 3,file = 'Q3.txt')
    open(unit = 4,file = 'Q4.txt')
    open(unit = 5,file = 'Q5.txt')
    open(unit = 6,file = 'Q6.txt')
    open(unit = 7,file = 'Q7.txt')
    open(unit = 8,file = 'Q8.txt')
    
    
    do d = 1,dimPk
        do j = 1,Ny
            do i = 1,Nx
                write(1,*) uh(i,j,d,1)
                write(2,*) uh(i,j,d,2)
                write(3,*) uh(i,j,d,3)
                write(4,*) uh(i,j,d,4)
                write(5,*) uh(i,j,d,5)
                write(6,*) uh(i,j,d,6)
                write(7,*) uh(i,j,d,7)
                write(8,*) uh(i,j,d,8)
            end do
        end do
    end do
    
    close(1)
    close(2)
    close(3)
    close(4)
    close(5)
    close(6)
    close(7)
    close(8)
    
    end subroutine save_solution