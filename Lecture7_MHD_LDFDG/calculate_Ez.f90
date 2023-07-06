    subroutine calculate_Ez(Ez,u1,u2,B1,B2)
    
    real Ez,u1,u2,B1,B2
    
    Ez = u2*B1 - u1*B2
    
    end subroutine calculate_Ez