    subroutine get_basis
    
    include 'com.txt'
    
    phiRU(1) = 1
    phiRU(2) = 1
    phiRU(3) = 1
    phiRU(4) = 2d0/3d0
    phiRU(5) = 1
    phiRU(6) = 2d0/3d0
    
    phiLU(1) = 1
    phiLU(2) = -1
    phiLU(3) = 1
    phiLU(4) = 2d0/3d0
    phiLU(5) = -1
    phiLU(6) = 2d0/3d0
    
    phiRD(1) = 1
    phiRD(2) = 1
    phiRD(3) = -1
    phiRD(4) = 2d0/3d0
    phiRD(5) = -1
    phiRD(6) = 2d0/3d0
    
    phiLD(1) = 1
    phiLD(2) = -1
    phiLD(3) = -1
    phiLD(4) = 2d0/3d0
    phiLD(5) = 1
    phiLD(6) = 2d0/3d0
    
    do i = 1,NumGLP
        do j = 1,NumGLP
            phiG(i,j,1) = 1
            phiGLL(i,j,1,1) = 1
            phiGLL(i,j,1,2) = 1
            phiGR(j,1) = 1
            phiGL(j,1) = 1
            phiGU(i,1) = 1
            phiGD(i,1) = 1
            phixG(i,j,1) = 0
            phiyG(i,j,1) = 0
            mm(1) = 1
            
            phiG(i,j,2) = lambda(i)
            phiGLL(i,j,2,1) = lambdaL(i)
            phiGLL(i,j,2,2) = lambda(i)
            phiGR(j,2) = 1
            phiGL(j,2) = -1
            phiGU(i,2) = lambda(i)
            phiGD(i,2) = lambda(i)
            phixG(i,j,2) = 1d0/hx1
            phiyG(i,j,2) = 0
            mm(2) = 1d0/3d0
                
            phiG(i,j,3) = lambda(j)
            phiGLL(i,j,3,1) = lambda(j)
            phiGLL(i,j,3,2) = lambdaL(j)
            phiGR(j,3) = lambda(j)
            phiGL(j,3) = lambda(j)
            phiGU(i,3) = 1
            phiGD(i,3) = -1
            phixG(i,j,3) = 0
            phiyG(i,j,3) = 1d0/hy1
            mm(3) = 1d0/3d0
                
            phiG(i,j,4) = lambda(i)**2 - 1d0/3d0
            phiGLL(i,j,4,1) = lambdaL(i)**2 - 1d0/3d0
            phiGLL(i,j,4,2) = lambda(i)**2 - 1d0/3d0
            phiGR(j,4) = 2d0/3d0
            phiGL(j,4) = 2d0/3d0
            phiGU(i,4) = lambda(i)**2 - 1d0/3d0
            phiGD(i,4) = lambda(i)**2 - 1d0/3d0
            phixG(i,j,4) = 2d0*lambda(i)/hx1
            phiyG(i,j,4) = 0
            mm(4) = 4d0/45d0
                    
            phiG(i,j,5) = lambda(i)*lambda(j)
            phiGLL(i,j,5,1) = lambdaL(i)*lambda(j)
            phiGLL(i,j,5,2) = lambda(i)*lambdaL(j)
            phiGR(j,5) = lambda(j)
            phiGL(j,5) = -lambda(j)
            phiGU(i,5) = lambda(i)
            phiGD(i,5) = -lambda(i)
            phixG(i,j,5) = lambda(j)/hx1
            phiyG(i,j,5) = lambda(i)/hy1
            mm(5) = 1d0/9d0
                    
            phiG(i,j,6) = lambda(j)**2 - 1d0/3d0
            phiGLL(i,j,6,1) = lambda(j)**2 - 1d0/3d0
            phiGLL(i,j,6,2) = lambdaL(j)**2 - 1d0/3d0
            phiGR(j,6) = lambda(j)**2 - 1d0/3d0
            phiGL(j,6) = lambda(j)**2 - 1d0/3d0
            phiGU(i,6) = 2d0/3d0
            phiGD(i,6) = 2d0/3d0
            phixG(i,j,6) = 0
            phiyG(i,j,6) = 2d0*lambda(j)/hy1
            mm(6) = 4d0/45d0
        end do
        
    end do
    
    phiGLL(:,:,:,3) = phiG
    
    end subroutine get_basis