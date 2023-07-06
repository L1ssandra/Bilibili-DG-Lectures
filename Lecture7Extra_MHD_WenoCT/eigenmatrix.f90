    subroutine eigenmatrix(evr,evl,rho,u1,u2,u3,E,B1,B2,B3,n1,n2,n3)
    
    
    real(kind = 8) evl(8,8),evr(8,8),rho,u1,u2,u3,B1,B2,B3,n1,n2,n3,t1,t2,t3,p,gamma
    real(kind = 8) ru(8,8),lu(8,8)
    real(kind = 8) rhosq,a2,a,d,cf,cs,beta1,g1,sqg2,sq12
    real(kind = 8) alphas,alphaf,E,H,sqpr,sqpor,sq1og, bst
    real(kind = 8) b1s, b2s, b3s
    real(kind = 8) twosq
    real(kind = 8) um2,nen,nen2,nen31,nen32,BN,nen3
    real(kind = 8) nen51,nen52,nen71,nen72,Term1,Term2,Term3,Term4
    real(kind = 8) Term5,Term6,Term7,Term8,Term9,Term10,Term11
    real(kind = 8) Term12,Term13,Term14,Term15,Term16,Term17,Term18
    real(kind = 8) psq,sqgam,gu1sq
    real(kind = 8) TxN1,TxN2,TxN3,BT,BNs
    real(kind = 8) Term51,Term52,Term53,Term54,Term55,Term56,Term57
    real(kind = 8) Term71,Term72,Term73,Term74,Term75,Term76,Term77
    integer :: m
    parameter(gamma = 5.d0/3.d0)
    
    p = (gamma - 1.d0)*(E - 0.5d0*rho*(u1**2 + u2**2 + u3**2) - 0.5d0*(B1**2 + B2**2 + B3**2))
    
    ! R
    
    call compute_t_vector(n1, n2, n3, B1, B2, B3, t1, t2, t3)
    
    rhosq = sqrt(abs(rho))
    g1=gamma-1.d0
    a2 = abs(gamma*p/rho)
    a = sqrt(a2)

    sqg2=sqrt(1.d0/(2.d0*gamma))

    sq12 = sqrt(0.5d0)

    sqpr = sqrt(abs(p))/rho

    sqpor = sqrt(abs(p/rho))

    sq1og = sqrt(1.d0/gamma)

    b1s=B1/rhosq
    b2s=B2/rhosq
    b3s=B3/rhosq

    bN = b1s*n1+b2s*n2+b3s*n3

    d = a2 + (b1s**2 + b2s**2 + b3s**2)
    cf = sqrt(0.5d0*abs(d+sqrt(abs(d**2-4.d0*a2*(bN)**2))))
    cs = sqrt(0.5d0*abs(d-sqrt(abs(d**2-4.d0*a2*(bN)**2))))

    beta1 = dsign(1.d0,(b1s*n1+b2s*n2+b3s*n3)*1.d0)

    if ( abs(cf*cf-cs*cs).le. 1.d-12) then
       alphaf = dsin(datan(1.d0))
       alphas = dcos(datan(1.d0))
    else
       alphaf = sqrt(abs(a2 - cs*cs))/sqrt(abs(cf*cf-cs*cs))
       alphas = sqrt(abs(cf*cf - a2))/sqrt(abs(cf*cf-cs*cs))
    endif

    TxN1=n3*t2-n2*t3
    TxN2=n1*t3-n3*t1
    TxN3=n2*t1-n1*t2
    !
    !    # 1 - right eigenvector Entropy Wave
    ru(1,1)   =sqrt(g1/gamma)*rhosq
    ru(2,1)   =0.d0
    ru(3,1)   =0.d0
    ru(4,1)   =0.d0
    ru(5,1)   =0.d0
    ru(6,1)   =0.d0
    ru(7,1)   =0.d0
    ru(8,1)   =0.d0
    !
    !    # 2 - right eigenvector Divergence Wave
    ru(1,2)   =0.d0
    ru(2,2)   =0.d0
    ru(3,2)   =0.d0
    ru(4,2)   =0.d0
    ru(5,2)   =0.d0
    ru(6,2)   =sq1og*a*n1
    ru(7,2)   =sq1og*a*n2
    ru(8,2)   =sq1og*a*n3
    !
    !    # 3 - right eigenvector Alfven Wave
    !    lambda = V \times n + b \times n
    ru(1,3)   =0.d0
    ru(2,3)   =-sq12*(sqpr*(TxN1))
    ru(3,3)   =-sq12*(sqpr*(TxN2))
    ru(4,3)   =-sq12*(sqpr*(TxN3))
    ru(5,3)   =0.d0
    ru(6,3)   =sq12*sqpor*(TxN1)
    ru(7,3)   =sq12*sqpor*(TxN2)
    ru(8,3)   =sq12*sqpor*(TxN3)

    !
    !    # 4 - right eigenvector Alfven Wave
    !    lambda = V \times n - b \times n
    ru(1,4)   =ru(1,3)
    ru(2,4)   =-ru(2,3)
    ru(3,4)   =-ru(3,3)
    ru(4,4)   =-ru(4,3)
    ru(5,4)   =ru(5,3)
    ru(6,4)   =ru(6,3)
    ru(7,4)   =ru(7,3)
    ru(8,4)   =ru(8,3)

    !
    !    # 5 - right eigenvector
    !    lambda = V \times n + C_f
    bst = (b1s*t1+b2s*t2+b3s*t3)

    ru(1,5)   =sqg2*alphaf*rhosq
    ru(2,5)   =sqg2*((alphaf*a2*n1+alphas*a*(                     &
         (bst)*n1-(bN)*t1)))/(rhosq*cf)

    ru(3,5)   =sqg2*((alphaf*a2*n2+alphas*a*(                     &
         (bst)*n2-(bN)*t2)))/(rhosq*cf)

    ru(4,5)   =sqg2*((alphaf*a2*n3+alphas*a*(                     &
         (bst)*n3-(bN)*t3)))/(rhosq*cf)
    ru(5,5)   =sqg2*alphaf*rhosq*a2
    ru(6,5)   =sqg2*alphas*a*t1
    ru(7,5)   =sqg2*alphas*a*t2
    ru(8,5)   =sqg2*alphas*a*t3

    !
    !    # 6 - right eigenvector
    !    lambda = V \times n - C_f
    ru(1,6) = ru(1,5)
    ru(2,6) = -ru(2,5)
    ru(3,6) = -ru(3,5)
    ru(4,6) = -ru(4,5)
    ru(5,6) = ru(5,5)
    ru(6,6) = ru(6,5)
    ru(7,6) = ru(7,5)
    ru(8,6) = ru(8,5)

    !
    !    # 7 - right eigenvector
    !    lambda = V \times n + C_s
    ru(1,7) =sqg2*alphas*rhosq
    ru(2,7) =beta1*sqg2*(alphaf*cf**2*t1+a*n1*alphas*(bN))/(rhosq*cf)
    ru(3,7) =beta1*sqg2*(alphaf*cf**2*t2+alphas*a*(bN)*n2)/(rhosq*cf)
    ru(4,7) =beta1*sqg2*(alphaf*cf**2*t3+alphas*a*(bN)*n3)/(rhosq*cf)
    ru(5,7) =a**2*sqg2*alphas*rhosq
    ru(6,7) =-sqg2*alphaf*a*t1
    ru(7,7) =-sqg2*alphaf*a*t2
    ru(8,7) =-sqg2*alphaf*a*t3
    !
    !    # 8 - right eigenvector
    !    lambda = V \times n - C_s
    ru(1,8) =ru(1,7)
    ru(2,8) =-ru(2,7)
    ru(3,8) =-ru(3,7)
    ru(4,8) =-ru(4,7)
    ru(5,8) =ru(5,7)
    ru(6,8) =ru(6,7)
    ru(7,8) =ru(7,7)
    ru(8,8) =ru(8,7)

    !-   ------------------------------------------------------------------------
    !    CONSERVATIVE E_VECTORS
    !-   ------------------------------------------------------------------------
    !
    do m=1,8

       evr(1,m)=ru(1,m)/g1
       evr(2,m)=(ru(1,m)*u1 + ru(2,m)*rho)/g1
       evr(3,m)=(ru(1,m)*u2 + ru(3,m)*rho)/g1
       evr(4,m)=(ru(1,m)*u3 + ru(4,m)*rho)/g1
       evr(5,m)=(ru(5,m)/g1+B1*ru(6,m)+B2*ru(7,m)+B3*ru(8,m)+0.5d0*ru(1,m) &
            *(u1**2+u2**2+u3**2) + ru(2,m)*u1*rho + ru(3,m)*u2*rho     &
            +ru(4,m)*u3*rho)/g1
       evr(6,m)=ru(6,m)/g1
       evr(7,m)=ru(7,m)/g1
       evr(8,m)=ru(8,m)/g1

    enddo
    
    
    ! L
    
    call compute_t_vector(n1, n2, n3, B1, B2, B3, t1, t2, t3)
    
    rhosq = sqrt(rho)
    g1=gamma-1.d0
    a2 = abs(gamma*p/rho)
    a = sqrt(a2)

    b1s=B1/rhosq
    b2s=B2/rhosq
    b3s=B3/rhosq

    BNs = b1s*n1+b2s*n2+b3s*n3

    BN = (B1*n1+B2*n2+B3*n3)

    d = a2 + (b1s**2 + b2s**2 + b3s**2)
    cf = sqrt(0.5d0*abs(d+sqrt(abs(d**2-4.d0*a2*(BNs)**2))))
    cs = sqrt(0.5d0*abs(d-sqrt(abs(d**2-4.d0*a2*(BNs)**2))))

    if ( abs(cf*cf-cs*cs).le. 1.d-12) then
       alphaf = dsin(datan(1.d0))
       alphas = dcos(datan(1.d0))
    else
       alphaf = sqrt(abs(a2 - cs*cs))/sqrt(abs(cf*cf-cs*cs))
       alphas = sqrt(abs(cf*cf - a2))/sqrt(abs(cf*cf-cs*cs))
    endif


    psq   = sqrt(abs(p))
    twosq = sqrt(2.d0)
    um2 = 0.5d0*(u1**2 + u2**2 + u3**2)

    beta1 = dsign(1.d0,(BNs)*1.d0)

    sqgam = sqrt(g1/gamma)

    TxN1=n3*t2-n2*t3
    TxN2=n1*t3-n3*t1
    TxN3=n2*t1-n1*t2

    BT=(B1*t1+B2*t2+B3*t3)

    gu1sq=sqrt(1.d0/gamma)


    !    # ---------------------------------------------------
    !    # PRIMITIVE E-VECTORS
    !    # ---------------------------------------------------
    !
    !    # 1 - left eigenvector
    lu(1,1)   = 1.d0/(sqgam*rhosq)
    lu(1,2)   = 0.d0
    lu(1,3)   = 0.d0
    lu(1,4)   = 0.d0
    lu(1,5)   = -1.d0/(a2*sqgam*rhosq)
    lu(1,6)   = 0.d0
    lu(1,7)   = 0.d0
    lu(1,8)   = 0.d0

    !
    !    # 2 - left eigenvector
    nen = (n3**2*(t1**2+t2**2)-2.d0*n1*n3*t1*t3-2.d0*n2*t2*       &
         (n1*t1+n3*t3)+n2**2*(t1**2+t3**2)+n1**2*(t2**2+t3**2))

    nen2=a*gu1sq*nen

    lu(2,1)   = 0.d0
    lu(2,2)   = 0.d0
    lu(2,3)   = 0.d0
    lu(2,4)   = 0.d0
    lu(2,5)   = 0.d0
    lu(2,6)   = (-n2*t1*t2-n3*t1*t3+n1*(t2**2+t3**2))/nen2
    lu(2,7)   = (-t2*(n1*t1+n3*t3)+n2*(t1**2+t3**2))/nen2
    lu(2,8)   = (n3*(t1**2+t2**2)-(n1*t1+n2*t2)*t3)/nen2
    !
    !    # 3 - left eigenvector

    nen3 = sqrt(2.d0)*sqrt(abs(p))*nen
    nen31 = nen3/sqrt(abs(rho))

    lu(3,1)   = 0.d0
    lu(3,2)   = rho*(-TxN1)/nen3
    lu(3,3)   = rho*(-TxN2)/nen3
    lu(3,4)   = rho*(-TxN3)/nen3
    lu(3,5)   = 0.d0
    lu(3,6)   = TxN1/nen31
    lu(3,7)   = TxN2/nen31
    lu(3,8)   = TxN3/nen31
    !
    !    # 4 - left eigenvector

    lu(4,1)   = lu(3,1)
    lu(4,2)   = -lu(3,2)
    lu(4,3)   = -lu(3,3)
    lu(4,4)   = -lu(3,4)
    lu(4,5)   = lu(3,5)
    lu(4,6)   = lu(3,6)
    lu(4,7)   = lu(3,7)
    lu(4,8)   = lu(3,8)
    !
    !    # 5 - left eigenvector

    Term51 = rho*cf*((n2*t1*t2+n3*t1*t3-n1*(t2**2+t3**2))*rhosq*cf**2* &
         alphaf-a*BN*alphas*(-n2**2*t1+n1*n2*t2+n3*TxN2))

    Term52 = rho*cf*((-t2*(n1*t1+n3*t3)+n2*(t1**2+t3**2))*rhosq*cf**2 &
         *alphaf-a*BN*alphas*(-n1*n2*t1+n1**2*t2+n3*TxN1))

    Term53 = rho*cf*((n3*(t1**2+t2**2)-(n1*t1+n2*t2)*t3)*rhosq*cf**2* &
         alphaf-a*(-n1*n3*t1+n1**2*t3+n2*(-n3*t2+n2*t3))*BN*alphas)

    Term54 = alphaf/(twosq*a**2*gu1sq*rhosq*(alphaf**2+alphas**2))

    Term55 = (n2**2*t1-n1*n2*t2+n3*(n3*t1-n1*t3))*alphas

    Term56 = alphas*(-n1*n2*t1+n1**2*t2+n3*TxN1)
    Term57 = (-n1*n3*t1+n1**2*t3+n2*(-n3*t2+n2*t3))*alphas

    nen51 = twosq*a*nen*gu1sq*(a*BN**2*alphas**2+rhosq*cf**2*alphaf*  &
         (a*rhosq*alphaf+BT*alphas))

    nen52 = twosq*a*gu1sq*(alphaf**2+alphas**2)*nen


    lu(5,1)  = 0.d0
    lu(5,2)  = -Term51/nen51
    lu(5,3)  = Term52/nen51
    lu(5,4)  = Term53/nen51
    lu(5,5)  = Term54
    lu(5,6)  = Term55/nen52
    lu(5,7)  = Term56/nen52
    lu(5,8)  = Term57/nen52
    !
    !    # 6 - left eigenvector
    !
    lu(6,1)  = lu(5,1)
    lu(6,2)  = -lu(5,2)
    lu(6,3)  = -lu(5,3)
    lu(6,4)  = -lu(5,4)
    lu(6,5)  = lu(5,5)
    lu(6,6)  = lu(5,6)
    lu(6,7)  = lu(5,7)
    lu(6,8)  = lu(5,8)

    !    # 7 - left eigenvector
    Term71 = rho*cf*(a*(n2**2*t1-n1*n2*t2+n3*(n3*t1-n1*t3))*rhosq*  &
         alphaf+((B3*n2*t1-B2*n3*t1-B3*n1*t2+B2*n1*t3)*(-n3*t2  &
         +n2*t3)+B1*(n2**2*t1**2+n3**2*t1**2-2.d0*n1*n2*t1*t2-2.d0 &
         *n1*n3*t1*t3+n1**2*(t2**2+t3**2)))*alphas)

    Term72 = rho*cf*(((n3*t1-n1*t3)*(B3*n2*t1-B3*n1*t2+B1*n3*t2-B1  &
         *n2*t3)+B2*((n1**2+n3**2)*t2**2-2.d0*n2*t2*(n1*t1+n3*t3)  &
         +n2**2*(t1**2+t3**2)))*alphas+a*rhosq*alphaf*(-n1*n2*t1   &
         +n1**2*t2+n3*TxN1))

    Term73 = rho*cf*(a*(-n1*n3*t1+n1**2*t3+n2*(-n3*t2+n2*t3))*rhosq &
         *alphaf+alphas*(B3*(n3**2*(t1**2+t2**2)-2.d0*n3*      &
         (n1*t1+n2*t2)*t3+(n1**2+n2**2)*t3**2)+B2*             &
         (n3*t1-n1*t3)*TxN3+B1*(-n3*t2+n2*t3)*TxN3))

    Term74 = alphas/(twosq*a**2*gu1sq*rhosq*(alphaf**2+alphas**2))

    Term75 = -alphaf*(n2**2*t1-n1*n2*t2+n3*(n3*t1-n1*t3))

    Term76 = -alphaf*(-n1*n2*t1+n1**2*t2+n3*TxN1)

    Term77 = -alphaf*(-n1*n3*t1+n1**2*t3+n2*(-n3*t2+n2*t3))


    nen71 = twosq*beta1*nen                                         &
         *gu1sq*(a*BN**2*alphas**2+rhosq*cf**2*alphaf*(a*rhosq*  &
         alphaf+BT*alphas))

    nen72 = nen52

    lu(7,1)= 0.d0
    lu(7,2)= Term71/nen71
    lu(7,3)= Term72/nen71
    lu(7,4)= Term73/nen71
    lu(7,5)= Term74
    lu(7,6)= Term75/nen72
    lu(7,7)= Term76/nen72
    lu(7,8)= Term77/nen72

    !
    !    # 8 - left eigenvector
    lu(8,1)=lu(7,1)
    lu(8,2)=-lu(7,2)
    lu(8,3)=-lu(7,3)
    lu(8,4)=-lu(7,4)
    lu(8,5)=lu(7,5)
    lu(8,6)=lu(7,6)
    lu(8,7)=lu(7,7)
    lu(8,8)=lu(7,8)

    !
    !    # ---------------------------------------------------
    !    # CONSERVATIVE E-VECTORS
    !    # ---------------------------------------------------
    do m=1,8
       evl(m,1) =lu(m,1)*g1-lu(m,2)*u1*g1/rho-lu(m,3)*u2*g1/rho-   &
            lu(m,4)*u3*g1/rho+lu(m,5)*g1**2*(u1**2+u2**2+u3**2)*.5d0
       evl(m,2) =-lu(m,5)*u1*g1**2+lu(m,2)*g1/rho
       evl(m,3) =-lu(m,5)*u2*g1**2+lu(m,3)*g1/rho
       evl(m,4) =-lu(m,5)*u3*g1**2+lu(m,4)*g1/rho
       evl(m,5) =lu(m,5)*g1**2
       evl(m,6) =lu(m,6)*g1-B1*lu(m,5)*g1**2
       evl(m,7) =lu(m,7)*g1-B2*lu(m,5)*g1**2
       evl(m,8) =lu(m,8)*g1-B3*lu(m,5)*g1**2
    enddo
    
    
    
    
    end subroutine eigenmatrix
    
    
    
    subroutine compute_t_vector(n1, n2, n3, B1, B2, B3, t1, t2, t3)
    ! compute a vector t such that |t|=1, t.n = 0, and t is in span{n, B}
    ! assumption: |n|=1
    real(kind = 8) :: n1, n2, n3, B1, B2, B3
    real(kind = 8) :: t1, t2, t3
    real(kind = 8) eps
    parameter(eps = 1e-7) ! a small positive number
    real(kind = 8) :: nB1, nB2, nB3
    real(kind = 8) :: tp1, tp2, tp3, tpnorm   ! t, before normalization

    call compute_cross(n1, n2, n3, B1, B2, B3, nB1, nB2, nB3)

    if(nB1**2 + nB2**2 + nB3**2 < eps) then
       call pick_orthogonal(n1, n2, n3, tp1, tp2, tp3)
    else
       call compute_cross(nB1, nB2, nB3, n1, n2, n3, tp1, tp2, tp3)
    endif

    tpnorm = sqrt(tp1**2 + tp2**2 + tp3**2)
    t1 = tp1/tpnorm
    t2 = tp2/tpnorm
    t3 = tp3/tpnorm
    end subroutine compute_t_vector
    
    subroutine compute_cross(a1, a2, a3, b1, b2, b3, c1, c2, c3)
    ! compute cross product of a and b, and store result in c
    real(kind = 8) :: a1, a2, a3, b1, b2, b3
    real(kind = 8) :: c1, c2, c3

    c1 = a2*b3 - a3*b2
    c2 = a3*b1 - a1*b3
    c3 = a1*b2 - a2*b1
    end subroutine compute_cross

    real function norm_sq(v1, v2, v3)
    real(kind = 8) :: v1, v2, v3

    norm_sq = v1*v1 + v2*v2 + v3*v3
    end function norm_sq

    subroutine pick_orthogonal(v1, v2, v3, r1, r2, r3)
    ! return some r such that v.r = 0
    ! assumption:  norm of v is not very small
    real(kind = 8) :: v1, v2, v3
    real(kind = 8) :: r1, r2, r3
    if(abs(v1) >= abs(v2) .and. abs(v1) >= abs(v3)) then
       r2 = 1.0
       r3 = 1.0
       r1 = (-v2-v3) / v1
    elseif(abs(v2) >= abs(v1) .and. abs(v2) >= abs(v3)) then
       r1 = 1.0
       r3 = 1.0
       r2 = (-v1-v3) / v2
    elseif(abs(v3) >= abs(v1) .and. abs(v3) >= abs(v1)) then
       r1 = 1.0
       r2 = 1.0
       r3 = (-v1-v2) / v3
    endif
    end subroutine pick_orthogonal
    
    
    
    
    
    