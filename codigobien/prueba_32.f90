!prueba_32.f90
PROGRAM prueba_32

    !!!!!!!!!!!!!!!!!!Programa para calcular las colisiones!!!!!!!!!!!
    !!  Calculo de colisiones en 2d con condiciones periodicas      !!
    !!                                                              !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none



    INTEGER:: i,j,k,l
    REAL(kind=8),ALLOCATABLE,DIMENSION(:,:):: r,v !vector con posiciones y velocidades
    REAL(kind=8),ALLOCATABLE,DIMENSION(:,:,:)::sumv ! suma de velocidades para cada colision
    REAL(kind=8),ALLOCATABLE,DIMENSION(:,:):: tmp !temperaturas en las dos direcciones del espacio
    REAL(kind=8)::temp,tempz,H,longy,sigma !temperaturas en "y" y en "z", altura vertical y anchura, sigma==tamaño de la particula
    REAL(kind=8)::tcol,colt !tiempo de colision, tiempo para comparar y tiempo inicial
    INTEGER::rep,iter,n !numero de repeticiones que se realizan (tiempo) y numero de iteraciones  (numero de copias)
    REAL(kind=8),ALLOCATABLE,DIMENSION(:)::rab,vab
    INTEGER,DIMENSION(2)::ni
    INTEGER,ALLOCATABLE,DIMENSION(:)::colisiones
    REAL(kind=8)::bij,discr,t 

    !!! para deteminar el tiempo de cálculo
    !REAL(kind=4):: start, finish
    ! para ver si solapan en la configuracion
    !LOGICAL::ov
    !notar que consideramos KT=1
    !inicializamos variables
    temp=1.0d00
    tempz=5*temp
    sigma=1d00
    H=1.5*sigma
    n=500
    longy=REAL(n,8)/(0.06*(H-sigma))
    rep=550000
    iter=10
    ALLOCATE(r(n,2),v(n,2),sumv(iter,rep,2),tmp(rep,2),rab(2),vab(2),colisiones(iter))




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


    colisiones(:)=0
    i=1
     DO i=1,iter
    t=0.0
   
        CALL inicializacion2d(n,longy,H,temp,tempz,r,v)
        !PRINT*,'temperatura en y', 0.5d00*sum(v(:,1)**2)/n, 'temperatura en z' , 0.5d00*sum(v(:,2)**2)/n
        !  r(:,1) = r(:,1) / longy
        !  r(:,2)=r(:,2)/(H)                                    ! convertir en unidades de la caja
        !  r(:,1) = r(:,1) - ANINT ( r(:,1) )                        ! condiciones periodicas
        DO j=1,rep
            ni=0
            colt = HUGE(1.0)
            DO k=1,(n-1)
                rab(:)=r(k,:)-r(k+1,:) ! calculamos posiciones relativas
                rab(1)=rab(1)-longy*ANINT(rab(1)/(longy)) ! condiciones periodicas
                !  rab(1) = rab(1) * longy ! unidades correctas
                !  rab(2)=rab(2)*(H)
                vab(:)=v(k,:)-v(k+1,:)   !calculamos velocidades relativas
                bij    = DOT_PRODUCT ( rab, vab )   ! obtenemos el producto escalar (ri-rj)*(vi-vj)
                    IF (bij<0 ) THEN
                    discr=bij**2-(SUM(rab**2)-sigma**2)*SUM(vab**2)
                    IF( discr>0.0) THEN ! si colisiona con la sucesiva particula
                        tcol = ( -bij - SQRT ( discr ) ) / ( SUM ( vab**2 ) )
                            
                        IF (tcol<0) THEN 
                           PRINT*, 'colisión:',k,k+1,'tiempo',tcol
                        END IF 
                        
                        IF (tcol<colt ) THEN
                            colt=tcol
                            ni(1)=k
                            ni(2)=k+1
                            !PRINT*, 'parejas que colisionan',ni
                        END IF
                    END IF
                    END IF

                    
                    !!!!!!!!!!!!!!!!!!!!!!!!!!! Colisión con las paredes !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



                    IF (v(k,2)>0 ) THEN
                        tcol=(H-sigma*0.5-r(k,2))/v(k,2)

                        IF (tcol<0) THEN 
                        PRINT*, 'colisión:',k,' con pared. Tiempo',tcol
                        END IF 
                        !PRINT*,'tiempo de colisión', tcol
                            IF (tcol<colt ) THEN
                                colt=tcol
                                ni(1)=k
                                ni(2)=n+1
                                !PRINT*, 'parejas que colisionan',ni(1),'pared de arriba'
                            END IF
                    END IF
                    IF (v(k,2)<0 ) THEN
                        tcol=(sigma*0.5-r(k,2))/v(k,2)

                            IF (tcol<0 ) THEN 
                            PRINT*, 'colisión:',k,' con pared. Tiempo',tcol
                            END IF 

                            IF (tcol<colt) THEN
                            colt=tcol
                            ni(1)=k
                            ni(2)=n+2
                            !PRINT*, 'parejas que colisionan',ni(1),'pared de abajo'
                            END IF
                    END IF
                    !!!!!!!!!!!!!!!! Si consideramos la partícula 1, vemos si esta colisiona con la última!!!!!!!!!!!!!!    
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    IF (k==1) THEN
                        rab(:)=r(k,:)-r(n,:) ! calculamos posiciones relativas
                        rab(1)=rab(1)-longy*ANINT(rab(1)/(longy)) ! condiciones periodicas
                        ! rab(1) = rab(1) * longy ! unidades correctas
                        ! rab(2)=rab(2)*(H )
                        vab(:)=v(k,:)-v(n,:)   !calculamos velocidades relativas
                        bij    = DOT_PRODUCT ( rab, vab )   ! obtenemos el producto escalar (ri-rj)*(vi-vj)
                        IF (bij<0 ) THEN
                        discr=bij**2-(SUM(rab**2)-sigma**2)*SUM(vab**2)
                        IF( discr>0.0) THEN ! si colisiona con la sucesiva particula
                            tcol = ( -bij - SQRT ( discr ) ) / ( SUM ( vab**2 ) )

                            IF (tcol<0) THEN 
                                PRINT*, 'colisión:',k,' con',n,'. Tiempo',tcol
                            END IF 
                            
                            IF (tcol<colt ) THEN
                                colt=tcol
                                ni(1)=k
                                ni(2)=n
                                !PRINT*, 'parejas que colisionan',ni
                            END IF
                        END IF
                        END IF
            
            
            
                    
                
                
                    END IF


            END DO

            IF (colt<0) THEN 
                WRITE(*,*) 'TIEMPOS NEGATIVOS'
                PRINT*,  'iteracion', j
                STOP 
            END IF 
                !!!Hacer avanzar las posiciones hasta el momento de colisión
            t=t+colt
            r(:,1)     = r(:,1) + colt * v(:,1) 
            r(:,2)     = r(:,2) + colt* v(:,2)   ! Advance all positions by t (box=1 units)
            r(:,1)     = r(:,1) - longy*ANINT(r(:,1)/(longy)) ! Apply periodic boundaries

            
            !PRINT*,'particula',ni(1), 'posicion en y',r(ni(1),1),'particula',ni(2), 'posicion en y',r(ni(2),1)     ! Apply periodic boundaries

            !PRINT*, 'parejas que colisionan finalmente ',ni
            IF (ni(2)<=n .AND. ni(2)>0) THEN
            CALL collide(ni(1),ni(2))
            colisiones(i)=colisiones(i)+1
            END IF

            IF (ni(2)>n) THEN
            v(ni(1),2)=-v(ni(1),2)
            !colisiones(i)=colisiones(i)+1
            END IF
               DO l=1,2
                   sumv(i,j,l)=0.5d00*sum(v(:,l)**2)/n

               END DO



        END DO
        
        CALL superpuesto()
        !PRINT*, "numero de colisiones en la iteración ",i ,":", colisiones(i)
    END DO

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

 


    DEALLOCATE(r,v,sumv,tmp,rab,vab,colisiones)



    CONTAINS


   

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

        factor = DOT_PRODUCT ( rij, vij )
        vij    = -factor * rij

        v(a,:) = v(a,:) + vij
        v(b,:) = v(b,:) - vij
            ! PRINT*, "velocidad particula 1",v(i,:)
            ! PRINT*, "velocidad particula 2,",v(j,:)
      END SUBROUTINE collide

     SUBROUTINE superpuesto() 
        LOGICAL::super
        INTEGER::q

        super=.FALSE. !! 1 si es falso, 0 si es verdadero
        iloop:DO q=1,n
            IF(ABS(r(q+1,1)-r(q,1)) <( 0.95)) THEN
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

END PROGRAM prueba_32