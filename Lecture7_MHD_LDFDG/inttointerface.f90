    subroutine inttointerface
    
    include 'com.txt'
    
    real BxR(k + 1),BxL(k + 1),ByU(k + 1),ByD(k + 1)
    real Bxint(dimPk),Byint(dimPk)
    real a0R,a1R,a2R,a0L,a1L,a2L
    real b0U,b1U,b2U,b0D,b1D,b2D
    real a00,a10,a01,a20,a11,a02,a30,a12
    real b00,b10,b01,b20,b11,b02,b21,b03
    real rxy,ryx
    
    do i = 1,Nx
        do j = 1,Ny
            
            a00 = uh(i,j,1,5)
            a10 = uh(i,j,2,5)
            a01 = uh(i,j,3,5)
            a20 = uh(i,j,4,5)
            a11 = uh(i,j,5,5)
            a02 = uh(i,j,6,5)
            a30 = uh(i,j,7,5)
            a12 = uh(i,j,9,5)
            
            b00 = uh(i,j,1,6)
            b10 = uh(i,j,2,6)
            b01 = uh(i,j,3,6)
            b20 = uh(i,j,4,6)
            b11 = uh(i,j,5,6)
            b02 = uh(i,j,6,6)
            b21 = uh(i,j,8,6)
            b03 = uh(i,j,10,6)
            
            a0R = a00 + a10 + (2d0/3d0)*a20 + (2d0/5d0)*a30
            a1R = a01 + a11
            a2R = a02 + a12
            
            a0L = a00 - a10 + (2d0/3d0)*a20 - (2d0/5d0)*a30
            a1L = a01 - a11
            a2L = a02 - a12
            
            b0U = b00 + b01 + (2d0/3d0)*b02 + (2d0/5d0)*b03
            b1U = b10 + b11
            b2U = b20 + b21
            
            b0D = b00 - b01 + (2d0/3d0)*b02 - (2d0/5d0)*b03
            b1D = b10 - b11
            b2D = b20 - b21
            
            BxR(1) = a0R
            BxR(2) = a1R
            BxR(3) = a2R
    
            BxL(1) = a0L
            BxL(2) = a1L
            BxL(3) = a2L 
    
            ByU(1) = b0U
            ByU(2) = b1U
            ByU(3) = b2U
    
            ByD(1) = b0D
            ByD(2) = b1D
            ByD(3) = b2D
            
            Bx(i,j,:) = BxR
            Bx(i - 1,j,:) = BxL
            By(i,j,:) = ByU
            By(i,j - 1,:) = ByD
            
        end do
    end do
    
    
    end subroutine inttointerface