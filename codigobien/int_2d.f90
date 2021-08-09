subroutine inicializacion2d(npart,longy,alt,temp,tempz,r,v)

    !! En esta subrutina vamos a asignar posiciones y velocidades a todas
    !! las particulas del sistema de forma que consigamos las temperaturas
    !! temp y tempz que queramos y que la velocidad total sea nula
    !!  npart--->numero de particulas
    !! longx-->longitud en dimensiones y del confinamiento en el eje horizontal. 
    
    implicit none
    REAL(kind = 8), EXTERNAL :: dran_u
    INTEGER::npart
    INTEGER:: i,number
    REAL(kind=8)::longy,alt
    REAL(kind=8),DIMENSION(npart,2):: v
    REAL(kind=8),DIMENSION(npart,2)::r
    REAL(kind=8)::sumvsq,sumvsqz,temp,tempz
    REAL(kind=8),DIMENSION(2)::sumv
    REAL(kind=8):: x1,x2,y1,y2
    REAL(kind=8)::size
    REAL(kind=8), parameter :: pi = 4 * atan (1.0_8)
    !! inicialicemos los acumuladores 
         sumv(:)=0
         sumvsq=0
         sumvsqz=0
    
 

    !!!! introducimos las variables caracteristicas del sistema      
         
         number=npart +1 !te da aproximadamente la mitad de las particulas para distribuirlas en x e y
         size=longy/REAL(number,8) ! te divide el espacio de tal forma que las particulas se distribuyan uniformmemente
         
    !!!! establecemos las posiciones y velocidades de las particulas  
        DO  i=1,npart
            r(i,1)=i*size-longy/2.0d0
            r(i,2)=alt*0.5d00
            ! r(i,2)=0.d00
        END DO
           
        DO i=1,npart,2 
            ! CALL RANDOM_NUMBER(x1)
            ! CALL RANDOM_NUMBER(y1)
            ! CALL RANDOM_NUMBER(x2)
            ! CALL RANDOM_NUMBER(y2)
            x1=dran_u()
            y1=dran_u()
            x2=dran_u()
            y2=dran_u()

            !determina las velocidades a traves de una distribucion maxwelliana
            v(i,1)=sqrt(-2.*log(x1))*cos(2.*pi*y1)
            v(i,2)=sqrt(-2.*log(x2))*cos(2.*pi*y2)
            v(i+1,1)=sqrt(-2.*log(x1))*sin(2.*pi*y1)
            v(i+1,2)=sqrt(-2.*log(x2))*sin(2.*pi*y2)
        END DO
        DO i=1,2              
        sumv(i)=sum(v(:,i))
        END DO 
    !!!! Corregimos las velocidades para que el momento total del sistema sea cero
    
        DO i=1,npart
            v(i,1)=v(i,1)-sumv(1)/dble(npart)
            v(i,2)=v(i,2)-sumv(2)/dble(npart)
            sumvsq=sumvsq+v(i,1)**2
            sumvsqz=sumvsqz+v(i,2)**2
        END DO
     
    !!!! Escalamos las velocidades para que las temperaturas sean las deseadas
    
        DO  i=1,npart !!multiplicamos por la raiz cuadrada de la temperatura /temperatura promedio <v^2>/2=v^2/2n (para particulas con masa unidad)
            v(i,1)=sqrt(dble(npart)*temp/sumvsq)*v(i,1) 
            v(i,2)=sqrt(dble(npart)*tempz/sumvsqz)*v(i,2)
        END DO
        
END SUBROUTINE inicializacion2d
    !!!! A partir de aquí está la subrutina para generar números aleatorios
    !!!! Hay que asignar un valor Ijkl en el programa principal      

    Subroutine Ran_Init(Ijkl)
           Implicit None
           Integer Ijkl
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !       This Is The Initialization Routine For Ran_Uniform()                  !
    !       Note: The Seed Variable Should Be In The Range 0 <= 900 000 000       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
           Integer I,J,K,L,Ij,Kl
           Integer Ii,Jj,M
           Double Precision S,T
    
           Double Precision    U(97),C,Cd,Cm
           Integer I97,J97
           Logical Initialised
           Common  /Raset1/ U,C,Cd,Cm,I97,J97,Initialised
    
           Initialised = .False.
    
            IF( (Ijkl <0) .Or.  (Ijkl > 900000000)) THEN
               Write(6,*) 'The Random Number Seed Must Have A Value Between 0 And 900 000 000'
                 Stop
           
            END IF
    
           Ij = Ijkl / 30082
           Kl = Ijkl - 30082 * Ij
    
           I = Mod(Ij/177,177) + 2
           J = Mod(Ij    ,177) + 2
           K = Mod(Kl/169,178) + 1
           L = Mod(Kl,   169)
    
           Do Ii = 1,97
              S = 0.0d0
              T = 0.5d0
              Do Jj = 1,24
                 M = Mod(Mod(I*J,179)*K,179)
                 I = J
                 J = K
                 K = M
                 L = Mod(53*L+1,169)
                 If (Mod(L*M,64) .Ge. 32) S = S + T
                 T = 0.5d0 * T
              End Do
              U(Ii) = S
           End Do
    
           C           = 362436.0d0 / 16777216.0d0
           Cd          = 7654321.0d0 / 16777216.0d0
           Cm          = 16777213.0d0 /16777216.0d0
           I97         = 97
           J97         = 33
           Initialised = .True.
    
           Return
    
    END SUBROUTINE Ran_Init