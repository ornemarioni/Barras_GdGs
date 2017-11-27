program boludez
    implicit none    
    REAL, dimension (3)   :: r, m,p
    REAL, PARAMETER     :: G = 4.299e-6
    REAL                :: aux, dist
    INTEGER             :: i, j, l,n       
    
    n= 3
    r(1)=18.26794243
    r(2)=17.41469765
    r(3)=17.70786858 
    m(1) = 3.02298827e-06
    m(2) = 3.02298827e-06
    m(3) = 3.02298827e-06
    DO i = 1, 3  
        p(i) = 0.
    END DO
    DO i = 1, 3
        DO j = 1, n
            dist = abs(r(i)-r(j))
            IF (i /= j .AND. dist /=0. ) THEN     
                   aux = G*m(i)*m(j)*1.e20/dist
                   p(i) = aux + p(i) 
            END IF
        END DO
    END DO
    print *,  p
end program 
