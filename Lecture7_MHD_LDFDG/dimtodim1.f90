    subroutine dimtodim1(U,U1)
    
    real U(6,1),U1(8,1)
    
    U1 = 0
    
    U1(1:3,1) = U(1:3,1)
    U1(5:7,1) = U(4:6,1)
    
    end subroutine dimtodim1
    
    
    