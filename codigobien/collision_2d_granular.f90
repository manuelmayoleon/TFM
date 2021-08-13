!final_version.f90
PROGRAM final_version

    !!!!!!!!!!!!!!!!!!Programa para calcular las colisiones!!!!!!!!!!!
    !!  Calculo de colisiones en 2d con condiciones periodicas      !!
    !!                                                              !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !*****************************************************************************80
  !  Modified:
  !
  !    19 Jule 2021
  !
  !  Author:
  !
  !    Manuel Mayo León 
  !



implicit none



    INTEGER:: i,j,k,l,m
    REAL(kind=8),ALLOCATABLE,DIMENSION(:,:):: r,v !vector con posiciones y velocidades
    REAL(kind=8),ALLOCATABLE,DIMENSION(:,:,:)::sumv ! suma de velocidades para cada colision
    REAL(kind=8),ALLOCATABLE,DIMENSION(:,:):: tmp !temperaturas en las dos direcciones del espacio
    REAL(kind=8)::temp,tempz,H,longy,sigma,epsilon,rho !temperaturas en "y" y en "z", altura vertical y anchura, sigma==tamaño de la particula, rho=density
    REAL(kind=8)::alpha, vp  !! coeficiente de restitucion y velocidad que se introduce a traves de la pared 
    LOGICAL :: boolean,granular
    REAL(kind=8)::tcol,colt !tiempo de colision, tiempo para comparar y tiempo inicial
    INTEGER::rep,iter,n,iseed !numero de repeticiones que se realizan (tiempo) y numero de iteraciones  (numero de copias)
    REAL(kind=8),ALLOCATABLE,DIMENSION(:)::rab,vab !distancias y velocidades relativas
    INTEGER,DIMENSION(2)::ni !particulas que colisionan
    INTEGER,ALLOCATABLE,DIMENSION(:)::colisiones !numero de colisiones
    REAL(kind=8)::bij,qij,discr,t !bij=(ri-rj)*(vi-vj), discr es el discriminante de la solucion de segundo grado, t=tiempo de colision
    REAL(kind=8),ALLOCATABLE,DIMENSION(:)::tiempos,deltas !tiempos de colision
    REAL(kind=8), parameter :: pi = 4 * atan (1.0_8)
    !!! para deteminar el tiempo de cálculo
    REAL(kind=4):: start, finish
    character(len=10)::alfa,eps


    !notar que consideramos KT=1
    !inicializamos variables
    temp=1.0d00
    ! tempz=0.d001*temp
    tempz=5.0*temp
    ! temp=1.d00
    ! tempz=5.d00
    
    sigma=1.0d00
    H=1.5*sigma
    n=500
    ! rho=0.06d00

    rho=0.015d00
    epsilon=(H-sigma)/sigma
    longy=REAL(n,8)/(rho*(H-sigma))
    ! rep=550000
    rep=30000000
    iter=1

    alpha=0.95
    ! alpha=1.0
    ! vp=0.001*temp
    vp=0.0001
    ! vp=0.0d0
  


    ALLOCATE(r(n,2),v(n,2),sumv(iter,rep,2),tmp(rep,2),rab(2),vab(2),colisiones(iter),tiempos(rep),deltas(rep))

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '            MD 2D SIMULATION                 '
    write ( *, '(a)' ) '            FORTRAN90 version                '
    
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Temperature y axis = ', temp
    write ( *, '(a,g14.6)' ) '  Temperature z axis = ', tempz
    write ( *, '(a,i8)' ) '  number of steps = ', rep
    write ( *, '(a,i8)' ) &
      '  The number of iterations taken  = ', iter
      write ( *, '(a,i8)' ) '  N = ', n
    write ( *, '(a,g14.6)' ) '  diameter (sigma) = ', sigma
    write ( *, '(a,g14.6)' ) '  density (rho) = ', rho
    write ( *, '(a,g14.6)' ) '  epsilon = ', epsilon
    write ( *, '(a,g14.6)' ) '  alpha = ', alpha
    write ( *, '(a,g14.6)' ) '  v_p = ', vp
     
    write ( *, '(a)' ) ' '


    WRITE(alfa,'(F10.2)')   alpha
    WRITE(eps,'(F10.2)')   epsilon
    !!!! para guardar los valores de las posiciones y velocidades iniciales!!!!!!

    call save_initial_distribution()

  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! llamamos al generador de numeros aleatorios
    iseed=2312
    call dran_ini(iseed)
    !pongo el numero de colisiones a cero
    colisiones(:)=0.d00
    !inicializo el tiempo para calcular el tiempo de computo
    call cpu_time(start)
    !Abro los archivos en los que voy a gaurdar los valores de las temperaturas, velocidades, posiciones, etc...
    
    ! If granular eqv true, the particles lost energy on each collision (inelastic disks) 
    !else, the collision is elastic between particles 
    granular=.TRUE.
    ! granular=.FALSE.
    
    ! If boolean eqv false, the particle is confined between two rigid plates 
    !else, the lower plate is vibrating in a sawtooth way 
    boolean=.TRUE.
    ! boolean=.FALSE.

    DO i=1,iter
        !inicializo los tiempos 
        t=0.0
        deltas(1)=0.0
        CALL inicializacion2d(n,longy,H,temp,tempz,r,v)
        !PRINT*,'temperatura en y', 0.5d00*sum(v(:,1)**2)/n, 'temperatura en z' , 0.5d00*sum(v(:,2)**2)/n
        DO j=1,rep
            ni=0
            colt = HUGE(1.0)
            DO k=1,n

                CALL calcular_tiempo_de_colision(k)


            END DO
         
            !para saber si los tiempos son negativos de nuevo
            IF (colt<0) THEN 
                WRITE(*,*) 'TIEMPOS NEGATIVOS'
                PRINT*,  'iteracion', j
                STOP 
            END IF 
            !!!Hacer avanzar las posiciones hasta el momento de colisión
            t=t+colt
            tiempos(j)=t

            CALL evolve_positions()
        

            !colision entre particulas  
            IF (ni(2)<=n .AND. ni(2)>0) THEN
             
            CALL collide(ni(1),ni(2))
            colisiones(i)=colisiones(i)+1
            END IF

            !colision entre particula a y muro
            IF (ni(2)>n) THEN

                CALL wall_collide(ni(1),ni(2))            

                   

            !colisiones(i)=colisiones(i)+1
            END IF

            ! OBTENEMOS LOS VALORES DE LAS TEMPERATURA Y TIEMPO MEDIO ENTRE COLISIONES
               DO l=1,2
                !    sumv(i,j,l)=0.5d00*sum(v(:,l)**2)/n
                sumv(i,j,l)=sum(v(:,l)**2)/n

               END DO
               IF (j>1) THEN
               deltas(j)=(0.06*sigma*epsilon*(sumv(i,j,1)-sumv(i,1,1))*(tiempos(j)))/(sqrt(pi*sumv(i,j,1)))
               END IF
        END DO

    
        CALL superpuesto()
       

        ! Guardamos los valores de las velocidades para todas las trayectorias de forma consecutiva. Asi aumentamos la estadística

        OPEN(11,FILE='velocidad_' // trim(adjustl(alfa)) // '.txt',  FORM ='FORMATTED',STATUS='UNKNOWN',POSITION='APPEND'&
        ,ACTION='READWRITE')   
                DO l=1,n
                    WRITE(11,*) v(l,1), v(l,2)
                END DO
            CLOSE(11)
       


     

        !PRINT*, "numero de colisiones en la iteración ",i ,":", colisiones(i)
    END DO

    DO l=1,rep
        DO m=1,2
        tmp(l,m)=sum(sumv(:,l,m))/iter
        ! tmp(l,m)=2*tmp(l,m)/(temp+tempz)
        END DO
    END DO 

 
   
    
    

    OPEN(9,FILE='temperaturas_' // trim(adjustl(alfa)) // '_' // trim(adjustl(eps)) // '.txt',STATUS='unknown')
    DO l=1,rep
        WRITE(9,*) tmp(l,1), tmp(l,2)
    END DO 
    CLOSE(9)   

        

    OPEN(10,FILE='pos_' // trim(adjustl(alfa)) // '.txt',STATUS='unknown')   
        DO l=1,n
         WRITE(10,*) r(l,1), r(l,2)
        END DO
    CLOSE(10) 

    
    OPEN(12,FILE='tiemposdecol_' // trim(adjustl(alfa)) // '.txt',STATUS='unknown') 
    DO l=1,rep
        WRITE(12,*) tiempos(l)
    END DO
    CLOSE(12) 
    
    OPEN(13,FILE='sparam_' // trim(adjustl(alfa)) // '.txt',STATUS='unknown')
    DO l=1,rep
        WRITE(13,*) deltas(l)
    END DO
    CLOSE(13)
   
    call save_data_file()
 
    !final del programa
    call cpu_time(finish)
    !calcula el tiempo de computo
    WRITE(*,*) '(Tiempo = ', finish-start , 'segundos.)'

    

   


    DEALLOCATE(r,v,sumv,tmp,rab,vab,colisiones,tiempos,deltas)



    CONTAINS
        !!!! para guardar los valores de las posiciones y velocidades iniciales!!!!!!
        subroutine save_initial_distribution()
          
            OPEN(7,FILE='velocidad_init.txt',STATUS='unknown')                       
            OPEN(8,FILE='posiciones_init.txt',STATUS='unknown')                      
            CALL inicializacion2d(n,longy,H,temp,tempz,r,v)                         
            DO i=1,n                                                                
                WRITE(7,*) v(i,1), v(i,2)
                WRITE(8,*)  r(i,1), r(i,2)
            END DO
            CLOSE(7)
            CLOSE(8)

        end subroutine save_initial_distribution


        subroutine save_data_file()
            implicit none 
            OPEN(UNIT=35,FILE='data.txt', FORM ='FORMATTED',STATUS='UNKNOWN',POSITION='APPEND'&
            ,ACTION='READWRITE')
    
            write(35,* ) 'N', n, 'T_y', temp, 'T_z', tempz, 'epsilon', epsilon,  ' alpha    '&
            , alpha, ' vp ', vp  , ' rho ', rho 
            write(35,* ) ' '
          
        end subroutine save_data_file
   

        SUBROUTINE collide ( a, b)
            IMPLICIT NONE
            INTEGER, INTENT(in)  :: a, b   ! Colliding atom indices
            ! LOGICAL :: granular
            ! This routine implements collision dynamics, updating the velocities
            ! The colliding pair (i,j) is assumed to be in contact already

            REAL(kind=8), DIMENSION(2) :: rij, vij
            REAL(kind=8)               :: factor

            rij(:) = r(a,:) - r(b,:)
            rij(1) = rij(1) - longy*ANINT ( rij(1)/(longy) ) ! Separation vector
            vij(:) = v(a,:) - v(b,:)           ! Relative velocity

            factor = DOT_PRODUCT ( rij, vij )
            if(granular .EQV. .TRUE.) THEN 
                vij    = -((1.0d0+alpha)*factor * rij)/(2.0d00)
            ELSE
                vij    = -factor * rij
            END IF

            v(a,:) = v(a,:) + vij
            v(b,:) = v(b,:) - vij
                ! PRINT*, "velocidad particula 1",v(i,:)
                ! PRINT*, "velocidad particula 2,",v(j,:)
        END SUBROUTINE collide
    
        subroutine calcular_tiempo_de_colision(c) 
            IMPLICIT NONE

            INTEGER, INTENT(in)  :: c
            
            

            IF(c/=n) THEN    
                rab(:)=r(c,:)-r(c+1,:) ! calculamos posiciones relativas
                rab(1)=rab(1)-longy*ANINT(rab(1)/(longy)) ! condiciones periodicas
                vab(:)=v(c,:)-v(c+1,:)   !calculamos velocidades relativas
                bij    = DOT_PRODUCT ( rab, vab )   ! obtenemos el producto escalar (ri-rj)*(vi-vj)

                !! FIRST WAY TO COMPUTE 

                    IF (bij<0 ) THEN
                    discr=bij**2-(SUM(rab**2)-sigma**2)*SUM(vab**2)
                    IF( discr>0.0) THEN ! si colisiona con la sucesiva particula
                        ! tcol = ( -bij - SQRT ( discr ) ) / ( SUM ( vab**2 ) )
                    !! ALTERNATIVE WAY 
                        
                        qij=-(bij+sign(1.0d00,bij)*dsqrt(discr))    
                        !  
                        tcol=MIN(qij/abs(dsqrt(sum(vab**2))),(sum(rab**2)-sigma**2)/qij )
                !comprobar que los tiempos no son negativos
                        IF (tcol<0) THEN 
                        PRINT*, 'colisión:',c,c+1,'tiempo',tcol
                        END IF 
                        
                        IF (tcol<colt ) THEN
                            ! PRINT*, 'colisión:',k,k+1,'tiempo',tcol
                            colt=tcol
                            ni(1)=c
                            ni(2)=c+1
                        END IF
                    END IF
                    END IF
      
            END IF
            
            
            !!!!!!!!!!!!!!!! Si consideramos la partícula 1, vemos si esta colisiona con la última!!!!!!!!!!!!!!    
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                    
            IF (c==1) THEN
                rab(:)=r(c,:)-r(n,:) ! calculamos posiciones relativas
                rab(1)=rab(1)-longy*ANINT(rab(1)/(longy)) ! condiciones periodicas
                vab(:)=v(c,:)-v(n,:)   !calculamos velocidades relativas
                bij    = DOT_PRODUCT ( rab, vab )   ! obtenemos el producto escalar (ri-rj)*(vi-vj)
                IF (bij<0 ) THEN
                discr=bij**2-(SUM(rab**2)-sigma**2)*SUM(vab**2)
                IF( discr>0.0) THEN ! si colisiona con la sucesiva particula
                    tcol = ( -bij - SQRT ( discr ) ) / ( SUM ( vab**2 ) )
                    
                    !comprobar que los tiempos no son negativos

                    IF (tcol<0) THEN 
                        PRINT*, 'colisión:',c,' con',n,'. Tiempo',tcol
                    END IF  
                    
                    IF (tcol<colt ) THEN
                        ! PRINT*, 'colisión:',k,' con',n,'. Tiempo',tcol
                        colt=tcol
                        ni(1)=c
                        ni(2)=n
                    END IF
                END IF
                END IF
    
    
    
            
        
        
            END IF

                
                !!!!!!!!!!!!!!!!!!!!!!!!!!! Colisión con las paredes !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



                IF (v(c,2)>0 ) THEN
                    tcol=(H-sigma*0.5-r(c,2))/v(c,2)

                    !comprobar que los tiempos no son negativos
                    IF (tcol<0) THEN 
                        PRINT*, 'colisión:',c,' con pared. Tiempo',tcol
                    END IF 
                
                        IF (tcol<colt ) THEN
                            colt=tcol
                            ni(1)=c
                            ni(2)=n+1
                        END IF
                END IF
                IF (v(c,2)<0 ) THEN
                    IF(boolean .EQV. .TRUE.) THEN 
                        tcol=(sigma*0.5-r(c,2))/(v(c,2)-vp)
                    ELSE
                        tcol=(sigma*0.5-r(c,2))/v(c,2)
                    END IF
                        !comprobar que los tiempos no son negativos
                        IF (tcol<0 ) THEN 
                            PRINT*, 'colisión:',k,' con pared. Tiempo',tcol
                        END IF 

                        IF (tcol<colt) THEN
                        colt=tcol
                        ni(1)=c
                        ni(2)=n+2
                        END IF
                END IF
                



        end subroutine calcular_tiempo_de_colision

        subroutine  evolve_positions()
            IMPLICIT NONE

            r(:,1)     = r(:,1) + colt * v(:,1) 
            r(:,2)     = r(:,2) + colt* v(:,2)   ! Advance all positions by t (box=1 units)
            r(:,1)     = r(:,1) - longy*ANINT(r(:,1)/(longy)) ! Apply periodic boundaries
        
        end subroutine evolve_positions
        !!!!! COLISIONES ENTRE UNA PARTICULA Y LA PARED. 
        ! granular=.FALSE. en el caso de colisiones con dos paredes rígidas
        ! granular =.TRUE. en el caso de que la pared de abajo sea de tipo diente de sierra 
        subroutine wall_collide(p,q)
            IMPLICIT NONE
            INTEGER :: p,q 
            ! LOGICAL :: boolean
            
            if (boolean .EQV. .FALSE.) THEN 
                v(p,2)=-v(p,2)
            ELSE 
                IF(q==n+1) THEN 
                    v(p,2)=-v(p,2)
                ELSE 
                    v(p,2)=2.0d00*vp-v(p,2)
                    ! print*, "collision with bottom wall"
                END IF 
            END IF
        end subroutine wall_collide

        SUBROUTINE superpuesto() 
        LOGICAL::super
        INTEGER::q

        super=.FALSE. !! 1 si es falso, 0 si es verdadero
        iloop:DO q=1,n
            IF(ABS(r(q+1,1)-r(q,1)) <( 0.95) .AND. ABS(r(q+1,2)-r(q,2)) <( 0.95)) THEN
                !PRINT*, 'particula ',q,'con',q+1,'superpuestas',ABS(r(q+1,1)-r(q,1))
            super=.TRUE.
            !PRINT*, 'particula',q,'distansia', ABS(r(q+1,1)-r(q,1))
                EXIT iloop !salir del loop 
            END IF
        END DO iloop

        IF (super.EQV. .TRUE.) THEN 
        PRINT*, "LAS PARTICULAS ESTAN SUPERPUESTAS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! :,("
        ELSE
            PRINT*, "LAS PARTICULAS ESTAN BIEN :)"
        END IF

        END SUBROUTINE superpuesto



END PROGRAM final_version