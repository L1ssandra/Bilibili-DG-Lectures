    program main
    
    include 'com.txt'
    
    call get_GLP
    
    call init_data
    
    call RK3
    
    call save_XY
    
    if (tend <= 5) then
        t20 = tend
    else
        t20 = 0
    end if
    
    call save_solution
    
    call calculate_L2error

    end program main

