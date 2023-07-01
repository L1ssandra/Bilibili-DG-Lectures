    subroutine calculate_totaldiv
    
    include 'com.txt'
    
    real div0
    real,allocatable :: edgediv(:,:)
    real,allocatable :: ujumpRL(:,:,:)
    real,allocatable :: ujumpUD(:,:,:)
    
    allocate(edgediv(Nx,Ny))
    allocate(ujumpRL(0:Nx,Ny,NumGLP))
    allocate(ujumpUD(Nx,0:Ny,NumGLP))
    
    totaldiv = 0
    open(unit = 10,file = 'div.txt')
    
    call set_bc
    
    UR = 0
    UL = 0
    UU = 0
    UD = 0
    
    edgediv = 0
    
    do i = 0,Nx
        do j = 1,Ny
            do d = 1,dimPk
                UR(i,j,:,6) = UR(i,j,:,6) + uh(i,j,d,6)*phiGR(:,d)
                UL(i + 1,j,:,6) = UL(i + 1,j,:,6) + uh(i + 1,j,d,6)*phiGL(:,d)
            end do
        end do
    end do
    
    do i = 1,Nx
        do j = 0,Ny
            do d = 1,dimPk
                UU(i,j,:,7) = UU(i,j,:,7) + uh(i,j,d,7)*phiGU(:,d)
                UD(i,j + 1,:,7) = UD(i,j + 1,:,7) + uh(i,j + 1,d,7)*phiGD(:,d)
            end do
        end do
    end do
    
    ujumpRL = UL(1:Nx1,:,:,6) - UR(:,:,:,6)
    ujumpUD = UD(:,1:Ny1,:,7) - UU(:,:,:,7)
    
    do j = 1,Ny
        do i = 1,Nx
            do i1 = 1,NumGLP
                edgediv(i,j) = edgediv(i,j) + hx1*weight(i1)*abs(ujumpRL(i,j,i1))
                edgediv(i,j) = edgediv(i,j) + hx1*weight(i1)*abs(ujumpRL(i - 1,j,i1))
                edgediv(i,j) = edgediv(i,j) + hy1*weight(i1)*abs(ujumpUD(i,j,i1))
                edgediv(i,j) = edgediv(i,j) + hy1*weight(i1)*abs(ujumpUD(i,j - 1,i1))
            end do
        end do
    end do
    
    do j = 1,Ny
        do i = 1,Nx
            uGdiv = 0
            div0 = 0
            do d = 1,dimPk
                uGdiv = uGdiv + uh(i,j,d,6)*phixG(:,:,d) + uh(i,j,d,7)*phiyG(:,:,d)
            end do
            
            do i1 = 1,NumGLP
                do j1 = 1,NumGLP
                    totaldiv = totaldiv + hx1*hy1*weight(i1)*weight(j1)*abs(uGdiv(i1,j1))
                    div0 = div0 + hx1*hy1*weight(i1)*weight(j1)*abs(uGdiv(i1,j1))
                end do
            end do
            
            totaldiv = totaldiv + edgediv(i,j)
            div0 = div0 + edgediv(i,j)
            
            write(10,*) div0
            
        end do
    end do
    
    close(10)
    
    end subroutine calculate_totaldiv