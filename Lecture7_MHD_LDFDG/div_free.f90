    subroutine div_free
    
    include 'com.txt'
    
    real BxR(k + 1),BxL(k + 1),ByU(k + 1),ByD(k + 1)
    real Bxint(dimPk),Byint(dimPk)
    real a0R,a1R,a2R,a0L,a1L,a2L
    real b0U,b1U,b2U,b0D,b1D,b2D
    real c00,c10,c01,c20,c11,c02
    real a00,a10,a01,a20,a11,a02,a30,a21,a12
    real b00,b10,b01,b20,b11,b02,b21,b12,b03
    real c5,c7,c8
    real a00new,a10new,a01new,a20new,a11new,a02new
    real b00new,b10new,b01new,b20new,b11new,b02new
    real rxy,ryx,ri,S1,S2
    
    do i = 1,Nx
        do j = 1,Ny
            
            a00 = uh(i,j,1,6)
            a10 = uh(i,j,2,6)
            a01 = uh(i,j,3,6)
            a20 = uh(i,j,4,6)
            a11 = uh(i,j,5,6)
            a02 = uh(i,j,6,6)
                
            b00 = uh(i,j,1,7)
            b10 = uh(i,j,2,7)
            b01 = uh(i,j,3,7)
            b20 = uh(i,j,4,7)
            b11 = uh(i,j,5,7)
            b02 = uh(i,j,6,7)
    
            ! The reconstruction of B = (Bx,By) from the interface
            a00new = a00
            a01new = a01
            a02new = a02
                
            b00new = b00
            b10new = b10
            b20new = b20
                
            c5 = (a10*hx - b01*hy)/(hx**2 + hy**2)
            c7 = (2*a20*hx - 5*b11*hy)/(2*hx**2 + 10*hy**2) 
            c8 = (5*a11*hx - 2*b02*hy)/(10*hx**2 + 2*hy**2)
                
            a10new = hx*c5
            a20new = hx*c7
            a11new = 2*hx*c8
                
            b01new = -hy*c5
            b11new = -2*hy*c7
            b02new = -hy*c8
                
            uh(i,j,1,6) = a00new
            uh(i,j,2,6) = a10new
            uh(i,j,3,6) = a01new
            uh(i,j,4,6) = a20new
            uh(i,j,5,6) = a11new
            uh(i,j,6,6) = a02new
            
            uh(i,j,1,7) = b00new
            uh(i,j,2,7) = b10new
            uh(i,j,3,7) = b01new
            uh(i,j,4,7) = b20new
            uh(i,j,5,7) = b11new
            uh(i,j,6,7) = b02new
        end do
    end do
    
    
    end subroutine div_free