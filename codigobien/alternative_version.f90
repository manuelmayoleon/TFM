!alternative_version.f90
PROGRAM alternative_version

    !!!!!!!!!!!!!!!!!!Programa para calcular las colisiones!!!!!!!!!!!
    !!  Calculo de colisiones en 2d con condiciones periodicas      !!
    !!                                                              !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none



    INTEGER:: i,j,k,l
    REAL(kind=8),ALLOCATABLE,DIMENSION(:,:):: r,v !vector con posiciones y velocidades
    REAL(kind=8),ALLOCATABLE,DIMENSION(:,:,:)::sumv ! suma de velocidades para cada colision
    REAL(kind=8),ALLOCATABLE,DIMENSION(:,:):: tmp !temperaturas en las dos direcciones del espacio
    REAL(kind=8)::temp,tempz,H,longy,sigma,epsilon,densidad !temperaturas en "y" y en "z", altura vertical y anchura, sigma==tamaño de la particula
    REAL(kind=8)::alpha !coeficiente de restitucion
    REAL(kind=8)::tcol,colt !tiempo de colision, tiempo para comparar y tiempo inicial
    INTEGER::rep,iter,n !numero de repeticiones que se realizan (tiempo) y numero de iteraciones  (numero de copias)
    REAL(kind=8),ALLOCATABLE,DIMENSION(:)::rab,vab !distancias y velocidades relativas
    INTEGER,DIMENSION(2)::ni !particulas que colisionan
    INTEGER,ALLOCATABLE,DIMENSION(:)::colisiones !numero de colisiones
    REAL(kind=8)::bij=0.d00,discr=0.d00,prodvij=0.d00,prodrij=0.d00,qij=0.d00,t !bij=(ri-rj)*(vi-vj), discr es el discriminante de la solucion de segundo grado, t=tiempo de colision
    REAL(kind=8),ALLOCATABLE,DIMENSION(:)::tiempos,deltas !tiempos de colision
    REAL(kind=8), parameter :: pi = 4 * atan (1.0_8)
    !!! para deteminar el tiempo de cálculo
    REAL(kind=4):: start, finish
   


    !notar que consideramos KT=1
    !inicializamos variables
    temp=1.0d00
    tempz=5*temp
    sigma=1d00
    H=1.5*sigma
    n=500
    densidad=0.06d0
    epsilon=(H-sigma)/sigma
    longy=REAL(n,8)/(densidad*H)
    rep=1000000
    iter=1

    alpha =0.2

    ALLOCATE(r(n,2),v(n,2),sumv(iter,rep,2),tmp(rep,2),rab(2),vab(2),colisiones(iter),tiempos(rep),deltas(rep))

    prodrij=0.0d0
    prodvij=0.0d0
    
    sumv=0.d00
    tmp=0.d00
    
    rab=0.d00
    vab=0.d00
    tiempos=0.d00
    deltas=0.d00
    
    !!!! para guardar los valores de las posiciones y velocidades iniciales!!!!!!
    ! OPEN(9,FILE='velocidad_init.txt',STATUS='unknown')                       
    ! OPEN(10,FILE='posiciones_init.txt',STATUS='unknown')                      
    ! CALL inicializacion2d(n,longy,H,temp,tempz,r,v)                         
    ! DO i=1,n                                                                
    !     WRITE(9,*) v(i,1), v(i,2)
    !     WRITE(10,*)  r(i,1), r(i,2)
    ! END DO
    ! CLOSE(9)
    ! CLOSE(10)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !pongo el numero de colisiones a cero
    colisiones(:)=0
    !inicializo el tiempo para calcular el tiempo de computo
    call cpu_time(start)
    
  
    DO i=1,iter
        !inicializo los tiempos 
        t=0.0
        deltas(1)=0.0
        CALL inicializacion2d(n,longy,H,temp,tempz,r,v)
        !PRINT*,'temperatura en y', 0.5d00*sum(v(:,1)**2)/n, 'temperatura en z' , 0.5d00*sum(v(:,2)**2)/n
        
        DO j=1,rep
            ni=0
            colt = HUGE(1.0)
            tcol=0.d00
            bij=0.d00
            discr=0.d00
            prodvij=0.d00
            prodrij=0.d00
            qij=0.d00

            DO k=1,n
                DO l=k+1,n
                        IF(k/=n) THEN 
                        rab(:)=r(k,:)-r(l,:) ! calculamos posiciones relativas
                        rab(1)=rab(1)-longy*ANINT(rab(1)/(longy)) ! condiciones periodicas
                        vab(:)=v(k,:)-v(l,:)   !calculamos velocidades relativas
                        bij    = producto_escalar ( rab, vab )   ! obtenemos el producto escalar (ri-rj)*(vi-vj)
                            IF (bij<0 ) THEN
                                prodvij=producto_escalar(vab,vab)
                                prodrij=producto_escalar(rab,rab)
                                ! PRINT*, prodrij
                                ! PRINT*, sum(rab**2)
                                discr=bij**2-(prodrij-sigma**2)*(prodvij)
                                ! discr=bij**2-(SUM(rab**2)-sigma**2)*SUM(vab**2)
                                IF( discr>0.0) THEN ! si colisiona con la sucesiva particula
                                
                                    ! tcol = ( -bij - sqrt( discr )) / ( prodvij )
                                
                                    !! ALTERNATIVE WAY 
                                
                                    qij=-(bij+sign(1.0d00,bij)*dsqrt(discr))    
                                
                                    tcol=MIN(qij/prodvij,(prodrij-sigma**2)/qij )
                                    
                                    !comprobar que los tiempos no son negativos
                                    IF (tcol<0) THEN 
                                        PRINT*, 'colisión:',k,k+1,'tiempo',tcol
                                    END IF 
                                
                                    IF (tcol<colt .AND. tcol> 0 ) THEN
                                        colt=tcol
                                        ni(1)=k
                                        ni(2)=k+1
                                    END IF
                                END IF
                            END IF
                        END IF
                        
                            

                    
                    
                    
                            
                END DO
                !!!!!!!!!!!!!!!!!!!!!!!!!!! Colisión con las paredes !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



                IF (v(k,2)>0 ) THEN
                    ! tcol=(H-sigma*0.5-r(k,2))/v(k,2)
                    tcol=((H-sigma)*0.5-r(k,2))/v(k,2)
                    !comprobar que los tiempos no son negativos
                    IF (tcol<0) THEN 
                    PRINT*, 'colisión:',k,' con pared. Tiempo',tcol
                    END IF 
                
                        IF (tcol<colt ) THEN
                            colt=tcol
                            ni(1)=k
                            ni(2)=n+1
                        END IF
                END IF
                IF (v(k,2)<0 ) THEN
                    ! tcol=(sigma*0.5-r(k,2))/v(k,2)
                    tcol=(-(H-sigma)*0.5-r(k,2))/v(k,2)
                        !comprobar que los tiempos no son negativos
                        IF (tcol<0 ) THEN 
                        PRINT*, 'colisión:',k,' con pared. Tiempo',tcol
                        END IF 

                        IF (tcol<colt) THEN
                        colt=tcol
                        ni(1)=k
                        ni(2)=n+2
                        END IF
                END IF


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
            r(:,1)     = r(:,1) + colt * v(:,1) 
            r(:,2)     = r(:,2) + colt* v(:,2)   ! Advance all positions by t (box=1 units)
            r(:,1)     = r(:,1) - longy*ANINT(r(:,1)/(longy)) ! Apply periodic boundaries

            !colision entre particulas  
            IF (ni(2)<=n .AND. ni(2)>0) THEN
            CALL collide(ni(1),ni(2))
            colisiones(i)=colisiones(i)+1
            END IF

            !colision entre particula a y muro
            IF (ni(2)>n) THEN
            v(ni(1),2)=-v(ni(1),2)
            !colisiones(i)=colisiones(i)+1
            END IF
               DO l=1,2
                   sumv(i,j,l)=0.5d00*sum(v(:,l)**2)/n

               END DO
            !    IF (j>1) THEN
            !    deltas(j)=(0.06*sigma*epsilon*(sumv(i,j,1)-sumv(i,1,1))*(tiempos(j)))/(sqrt(pi*sumv(i,j,1)))
            !    END IF
        END DO
        
        CALL superpuesto()
        PRINT*, "numero de colisiones en la iteración ",i ,":", colisiones(i)
    END DO

    !final del programa
    call cpu_time(finish)
    !calcula el tiempo de computo
    WRITE(*,*) '(Tiempo = ', finish-start , 'segundos.)'

        DO j=1,rep
                DO l=1,2
                tmp(j,l)=sum(sumv(:,j,l))/iter
                tmp(j,l)=2*tmp(j,l)/(temp+tempz)
                END DO
        END DO 

        OPEN(9,FILE='temperaturas.txt',STATUS='unknown')
        DO j=1,rep
        WRITE(9,*) tmp(j,1), tmp(j,2)
        END DO 
        CLOSE(9)



        OPEN(10,FILE='pos.txt',STATUS='unknown')
        DO i=1,n
         WRITE(10,*) r(i,1), r(i,2)
        END DO
        CLOSE(10)
        OPEN(11,FILE='velocidad.txt',STATUS='unknown')
        DO i=1,n
            WRITE(11,*) v(i,1), v(i,2)
        END DO
        CLOSE(11)


        OPEN(12,FILE='tiemposdecol.txt',STATUS='unknown')
        DO i=1,rep
            WRITE(12,*) tiempos(i)
        END DO
        CLOSE(12)
 
        ! OPEN(13,FILE='sparam.txt',STATUS='unknown')
        ! DO i=1,rep
        !     WRITE(13,*) deltas(i)
        ! END DO
        ! CLOSE(13)



    DEALLOCATE(r,v,sumv,tmp,rab,vab,colisiones,tiempos,deltas)



    CONTAINS

        FUNCTION producto_escalar(vector1,vector2) result(producto)
            implicit none
            REAL(kind=8), DIMENSION(2) :: vector1, vector2
            REAL(kind=8) :: producto
            INTEGER :: w
            producto=0.d00
            DO w=1,2
            producto=producto+vector1(w)*vector2(w)
            END DO 
            RETURN 
        END FUNCTION producto_escalar
   

      SUBROUTINE collide ( a, b)
        IMPLICIT NONE
        INTEGER, INTENT(in)  :: a, b   ! Colliding atom indices

        ! This routine implements collision dynamics, updating the velocities
        ! The colliding pair (i,j) is assumed to be in contact already

        REAL(kind=8), DIMENSION(2) :: rij, vij
        REAL(kind=8)               :: factor

        rij(:) = r(a,:) - r(b,:)
        rij(1) = rij(1) - longy*ANINT ( rij(1)/(longy) ) ! Separation vector
        vij(:) = v(a,:) - v(b,:)           ! Relative velocity
        rij(:)=rij(:)/dsqrt(producto_escalar(rij,rij))
        factor = producto_escalar ( rij, vij ) /(dsqrt(producto_escalar(rij,rij)))
        vij    = -(factor * rij)

        v(a,:) = v(a,:) + vij
        v(b,:) = v(b,:) - vij
            ! PRINT*, "velocidad particula 1",v(i,:)
            ! PRINT*, "velocidad particula 2,",v(j,:)
      END SUBROUTINE collide
      SUBROUTINE collide_granular ( a, b,alfa)
        IMPLICIT NONE
        INTEGER, INTENT(in)  :: a, b   ! Colliding atom indices
        ! This routine implements collision dynamics, updating the velocities
        ! The colliding pair (i,j) is assumed to be in contact already

        REAL(kind=8), DIMENSION(2) :: rij, vij
        REAL(kind=8)               :: factor,alfa

        rij(:) = r(a,:) - r(b,:)
        rij(1) = rij(1) - longy*ANINT ( rij(1)/(longy) ) ! Separation vector
        vij(:) = v(a,:) - v(b,:)           ! Relative velocity

        factor = DOT_PRODUCT ( rij, vij )
        vij    = -factor * rij

        v(a,:) = v(a,:) + vij*(1+alfa)/2d00
        v(b,:) = v(b,:) - vij*(1+alfa)/2d00
            ! PRINT*, "velocidad particula 1",v(i,:)
            ! PRINT*, "velocidad particula 2,",v(j,:)
      END SUBROUTINE collide_granular



     SUBROUTINE superpuesto() 
        LOGICAL::super
        INTEGER::q

        super=.FALSE. !! 1 si es falso, 0 si es verdadero
        iloop:DO q=1,n
            IF(ABS(r(q+1,1)-r(q,1)) <( 0.95) .AND. ABS(r(q+1,2)-r(q,2)) <( 0.95)  ) THEN
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

END PROGRAM alternative_version