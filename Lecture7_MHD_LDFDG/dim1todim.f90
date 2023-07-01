    subroutine dim1todim(U1,U)
    
    real U1(8,1),U(6,1)
    
    U(1:3,1) = U1(1:3,1)
    U(4:6,1) = U1(5:7,1)
    
    end subroutine dim1todim