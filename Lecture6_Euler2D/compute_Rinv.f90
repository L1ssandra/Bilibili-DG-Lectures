    subroutine compute_Rinv(Rmat,Rinv,rho,u,v,E,n1,n2)
    
    real nf(2),n1,n2
    real Rmat(4,4)
    real Rinv(4,4)
    real c,rho,rhou,rhov,u,v,pr,eH,ek,unf,magnorm
    
    gam = 1.4
    
    nf(1) = n1
    nf(2) = n2

    magnorm = dsqrt(nf(1)**2+nf(2)**2)
    nf(1) = nf(1)/magnorm
    nf(2) = nf(2)/magnorm
    
    pr  = (E-0.5d0*rho*(u**2+v**2))*(gam-1)
    c   = dsqrt(gam*pr/rho) ! speed of sound
    eH  = (E+pr)/rho ! specific enthalpy
    unf = u*nf(1)+v*nf(2)
    ek = 0.5d0*(u**2+v**2)
    
    Rmat(1,1) = 1.d0
    Rmat(2,1) = u-c*nf(1)
    Rmat(3,1) = v-c*nf(2)
    Rmat(4,1) = eH-c*unf
    
    Rmat(1,2) = 1.d0
    Rmat(2,2) = u
    Rmat(3,2) = v
    Rmat(4,2) = ek

    Rmat(1,3) = 1.d0
    Rmat(2,3) = u+c*nf(1)
    Rmat(3,3) = v+c*nf(2)
    Rmat(4,3) = eH+c*unf

    Rmat(1,4) = 0.d0
    Rmat(2,4) = nf(2)
    Rmat(3,4) = -nf(1)
    Rmat(4,4) = u*nf(2)-v*nf(1)
    
    ! ===========================
    Rinv(1,1) = ((gam-1)*ek+c*unf)*0.5d0/(c**2.d0)
    Rinv(2,1) = (c**2-(gam-1)*ek)/(c**2.d0)
    Rinv(3,1) = ((gam-1)*ek-c*unf)*0.5d0/(c**2.d0)
    Rinv(4,1) = v*nf(1)-u*nf(2)
    
    Rinv(1,2) = ((1-gam)*u-c*nf(1))*0.5d0/(c**2.d0)
    Rinv(2,2) = (gam-1)*u/(c**2.d0)
    Rinv(3,2) = ((1-gam)*u+c*nf(1))*0.5d0/(c**2.d0)
    Rinv(4,2) = nf(2)

    Rinv(1,3) = ((1-gam)*v-c*nf(2))*0.5d0/(c**2.d0)
    Rinv(2,3) = (gam-1)*v/(c**2.d0)
    Rinv(3,3) = ((1-gam)*v+c*nf(2))*0.5d0/(c**2.d0)
    Rinv(4,3) = -nf(1)

    Rinv(1,4) = (gam-1)*0.5d0/(c**2.d0)
    Rinv(2,4) = (1-gam)/(c**2.d0)
    Rinv(3,4) = (gam-1)*0.5d0/(c**2.d0)
    Rinv(4,4) = 0.d0
    
    end subroutine compute_Rinv    
    
    
    
    
    