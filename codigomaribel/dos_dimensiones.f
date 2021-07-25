ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc  Monocapa vibrada   
cccc
cccccccccccccccccccccccccccccccccc
ccc
ccc   Hemos elininado el calculo de la funcion de correlacion de la energia
ccc   Hemos introducido la modificacion para el tiempo real
ccc   Hemos definido ipvt,ipvt2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 

      program dinamicamolecular

      implicit none

      double precision temp,tempz,
     &densidad,timbig,intervalof,intervalog,ai,
     &w0,t0

      double precision rparedx

      integer npart,i,j,ranf,npasos,nceldas,maxf,maxg,ijkl,ptoinicial,
     &ipvt,ipvt2

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc  ptoinicial: ptoinicial*5=numero de colisiones por particula a partir 
cccc     de la que se empieza a promediar la funcion de correlacion espacial
cccc     (hay 6 elasticas)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      integer ndist

      parameter (ndist=100)

ccc   DENSIDAD SERIA LA DE LOS DISCOS EN 2D. ESTO HAY QUE CAMBIARLO

      parameter(npart=500,densidad=0.02d0,temp=1.0d0,tempz=5.0d0,
     &     timbig=100000000000000.0d0,ranf=50,npasos=2000000,
     &maxf=50,maxg=400,intervalof=0.1d0,ai=1.0d0,w0=0.0d0,
     &nceldas=4,ptoinicial=17,
     &ipvt=50,ipvt2=100)
     
cccc  Defino las longitudes del sistema

      double precision longx,longz

      parameter (longz=1.2d0)     

c      parameter(rparedx=1.0d0)

cccc  Defino la magnitud de la vibracion de la pared inferior

      double precision vibra
      parameter (vibra=0.0d0)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
ccc   ipvt=nos dice cuando empezamos a calcular la funcion de correlacion
ccc        de la velocidad transversal
ccc   ipvt2=nos dice cuando empezamos a calcular la funcion de correlacion
ccc        de la velocidad transversal por el segundo metodo (esperando a
ccc        que llegue al estacionario)
ccc 
ccc   Por definicion ipvt<ipvt2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision box,rx(npart),rz(npart),
     &vx(npart),vz(npart),
     &tcol(npart+2),tcolmin,vsq,sumvsq,tempi,g(maxg),
     &sumg(maxg),t,t1,dt,virial,sumvirial,gmas,sumgmas,sumgmas1,
     &gsigmamas(ranf),sumvx,sumvy,sumvx2,sumvy2,sumvz2,
     &sumvx4,sumvy4,sumv4,
     &sumvxvy,sumvx2vy2,sumvxmenosvy2,sumv2vx,
     &datos(npasos/npart+1,18),sumdatos(npasos/npart+1,18)


      double precision fvx(-maxf:maxf),fvy(-maxf:maxf),
     &fvz(-maxf:maxf),
     &fvxvy(-maxf:maxf,-maxf:maxf),
     &sumfvx(-maxf:maxf),mfvx(npasos/(5*ndist),-maxf:maxf+2),
     &summfvx(npasos/(5*ndist),-maxf:maxf+2),
     &sumfvy(-maxf:maxf),mfvy(npasos/(5*ndist),-maxf:maxf+2),
     &summfvy(npasos/(5*ndist),-maxf:maxf+2),
     &sumfvz(-maxf:maxf),mfvz(npasos/(5*ndist),-maxf:maxf+2),
     &summfvz(npasos/(5*ndist),-maxf:maxf+2),
     &sumfvxvy(-maxf:maxf,-maxf:maxf),
     &mfvxvy(npasos/(5*ndist),-maxf:maxf,-maxf:maxf),
     &summfvxvy(npasos/(5*ndist),-maxf:maxf,-maxf:maxf),

     &csuave(0:500),ssuave(0:500),
     &csuave2(0:500),ssuave2(0:500),

     &mgdr(10,maxg),summgdr(10,maxg),ecz,
     &datospg(npasos/(5*npart),3),sumdatospg(npasos/(5*npart),3),a,w,nf,
     &ec,flu(npasos/(5*npart),2),sumflu(npasos/(5*npart),2),t0menos1,
     &thcs,thcsexp1,thcsexp2,th(npasos/npart),pi,a2,
     &tiempo(npasos/npart),javier(ranf,npasos/npart),
     &javierxy(ranf,npasos/npart),javierv2x(ranf,npasos/npart),

     &sumjavier2(npasos/npart),
     &sumjavier3(npasos/npart),sumjavier4(npasos/npart),
     &sumjavier5(npasos/npart),sumjavier6(npasos/npart),

     &umediaceldax(nceldas,nceldas,npasos/(5*npart)),
     &umediacelday(nceldas,nceldas,npasos/(5*npart)),
     &nmediacelda(nceldas,nceldas,npasos/(5*npart)),
     &dEcmedia(nceldas,nceldas,npasos/(5*npart)),
     &Tempmedia(nceldas,nceldas,npasos/(5*npart)),
     &dTempmedia(nceldas,nceldas,npasos/(5*npart)),
     &corr(0:nceldas-1),correlacion(npasos/(npart),0:nceldas-1),
     &sumcorrelacion(npasos/(npart),0:nceldas-1),
     &sumtcorrelacion(0:nceldas-1),colp,corrymedio(nceldas),
     &correlacionymedio(npasos/(npart),nceldas),
     &sumcorrelacionymedio(npasos/(npart),nceldas),ux(nceldas),
     &uxpas(ranf,npasos/(npart),nceldas),
     &sumux(npasos/(npart),nceldas),densi(nceldas),
     &densipas(ranf,npasos/(npart),nceldas),
     &sumdensi(npasos/(npart),nceldas),temperature(npart),
     &temperaturepas(ranf,npasos/(npart),nceldas),
     &sumtemperature(npasos/(npart),nceldas),corry2(nceldas),
     &correlaciony2(npasos/(npart),nceldas),
     &sumcorrelaciony2(npasos/(npart),nceldas),
     &upromedio(nceldas),deltaudeltau(nceldas),
     &deltaudeltau2(npasos/(npart),nceldas),
     &densipromedio(nceldas),deltandeltan(nceldas),
     &deltandeltan2(npasos/(npart),nceldas),tempetotal



      double precision ycuadrado,vxpory,ecmod,sumv6

      integer partner(npart+2),ii,jj,p,q,k,l,m,c1,c2,c3,c4,c5,c6,
     &  contabilizador,contador2(ranf),c2max,contnuevo

      integer colisionpart,contadorpeli
      
      integer kk,ll

      double precision dencosx,dencosy,densinx,densiny
      double precision vxcosx,vxcosy,vxsinx,vxsiny
      double precision vycosx,vycosy,vysinx,vysiny,v2
      
      double precision acumpresion1,presion1,acumpresion2,presion2,
     &acumpresion3,presion3

      
      longx=npart/(densidad*longz)
     
cccccc inicializamos subrutina numeros aleatorios
cccccc Ijkl es la semilla que introducimos

       Ijkl=154
       call Ran_Init(Ijkl)

ccccc inicializamos la matriz tiempo()

      do p=1,npasos/npart
         tiempo(p)=timbig*timbig*timbig*timbig
      enddo    

      do p=1,npasos/npart
         sumjavier2(p)=0
         sumjavier3(p)=0
         sumjavier4(p)=0
         sumjavier5(p)=0
         sumjavier6(p)=0
      enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccc inicializamos funciones que acumulan entre los bucles cccc
ccccc de pasos y de trayectorias                            cccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


ccccc Matriz sumdatos

      do p=1,npasos/npart
         do q=1,18

            sumdatos(p,q)=0.0d0

         enddo
      enddo

ccccc  Funcion de distribucion

      
      do p=1,npasos/(5*ndist)
         do q=-maxf,maxf+2

            summfvx(p,q)=0.0d0
            summfvy(p,q)=0.0d0
            summfvz(p,q)=0.0d0

         enddo
      enddo

ccccc  Funcion de correlacion espacial g(r)
 
      do p=1,10
         do q=1,maxg

            summgdr(p,q)=0.0d0

         enddo
      enddo
       do p=1,ranf

          contador2(p)=0
c          print*,contador2(p)
       enddo

       
ccccc  Control de inestabilidades (manera antigua) y funcion de
ccccc  correlacion de las velocidades de una capa del fluido con otra

      do p=1,npasos/(5*npart)
         sumflu(p,1)=0.0d0
         sumflu(p,2)=0.0d0

         do k=0,nceldas-1
            sumcorrelacion(p,k)=0.0d0
         enddo
      do k=1,nceldas
            sumcorrelacionymedio(p,k)=0.0d0
            sumcorrelaciony2(p,k)=0.0d0
            sumux(p,k)=0.0d0
            sumdensi(p,k)=0.0d0
            sumtemperature(p,k)=0.0d0
         enddo
      enddo
      do p=1,ranf
         do q=1,npasos/(npart)
            do k=1,nceldas
               uxpas(p,q,k)=0.0d0
               densipas(p,q,k)=0.0d0
               temperaturepas(p,q,k)=0.0d0
            enddo
         enddo
      enddo





  
ccccc calculamos el tama�o de la caja y t0menos1

CCCC CUIDADO !!!!!!!!!!!!!!!!!!!

      box=(dble(npart)/densidad)**(1.0d0/2.0d0)

      intervalog=box/(2.0d0*dble(maxg))

      pi=3.141592654
 
      a2=16*(1-ai)*(1-2*ai*ai)/(57-25*ai+30*ai*ai-30*ai*ai*ai)

      t0menos1=(1-ai*ai)*densidad*sqrt(pi)/2*
     &(1+densidad*pi*(25-4*densidad*pi)/(4*(4-densidad*pi)**2))*
     &(1+3/16*a2)

    

      open(10,file='momentos0.02ranf100al90u01N2000.dat',
     &status='unknown')
      open(20,file='fvx0.02ranf100al90u01N2000.dat',status='unknown')
      open(25,file='gdr0.02ranf100al90u01N2000.dat',status='unknown')
      open(30,file='fluctuaciones0.02ranf100al90u01N2000.dat',
     &status='unknown')
     
      open(50,file='hist0.02ranf100al90u01N2000.dat',status='unknown')
      open(55,file='tempesranf100al90u01N2000.dat',status='unknown')
      open(44,file='tempromediaranf100al90u01N2000.dat',
     &status='unknown')
cc      open(13,file='autosuave0.02ranf100al90u01N2000.dat',
cc     &status='unknown')
      open(94,file='tmodos0.02ranf100al90u01N2000.dat',
     &status='unknown')
      open(95,file='cmodos0.02ranf100al90u01N2000.dat',
     &status='unknown')
      open(66,file='corres0.02ranf100al90u01N2000.dat',
     &status='unknown')
      open(63,file='pruebacorres0.02ranf100al90u01N2000.dat',
     &     status='unknown')
      open(64,file='correspromediada0.02ranf100al90u01N2000.dat',
     &     status='unknown')
      open(67,file='corresymedio0.02ranf100al90u01N2000.dat',
     &status='unknown')
      open(68,file='ux-densi-T-media0.02ranf100al90u01N2000.dat',
     &status='unknown')
      open(76,file='promediot0.02ranf100al90u01N2000.dat',
     &status='unknown')        
      open(78,file='vt_hist0.02ranf100al90u01N2000.dat',
     &status='unknown')
      open(79,file='vt_vt0.02ranf100al90u01N2000.dat',
     &status='unknown')
      open(80,file='vt_vt20.02ranf100al90u01N2000.dat',
     &status='unknown')

      open(93,file='energiainternaranf100al90u01N2000.dat',
     &status='unknown')


      c5=1





      do 2000 i=1,ranf
      
      
      print*,i

c      call flush()
CC     Gaussiana
c      call inicializacion(npart,longxy,temp,rx,ry,rz,vx,vy,vz,pi)
CC     Gaussiana tambien en Z
      call inicializacion0(npart,longx,temp,tempz,rx,rz,vx,vz,pi)
CC     Cuatro picos
c     call inicializacion2(npart,longxy,temp,rx,ry,rz,vx,vy,vz)
CC     Dos picos
c      call inicializacion3(npart,longxy,temp,rx,ry,rz,vx,vy,vz)
CC     Plana
c      call inicializacion4(npart,longxy,temp,rx,ry,rz,vx,vy,vz)
CC     Asimetrica VERIFICAR!!
c       call inicializacion5(npart,longxy,temp,rx,ry,rz,vx,vy,vz)
      
      do p=1,npart
c         print*,npart,p,rx(p),rz(p)
      write(89,*)npart,p,rx(p),rz(p)

      enddo
c       print*,npart,p,rx(500),rz(500)
       call tcolinicial(npart,longx,longz,timbig,rx,rz,vx,vz,
     &partner,tcol)

 
      a=1
      w=0
      t0=1
      c1=1
      c2=1
      c3=1
      c4=1
      c6=1
      contnuevo=0
      colisionpart=0
      contadorpeli=0
      presion1=0.0d0
      presion2=0.0d0
      presion3=0.0d0

cccc  inicializamos el tiempo en la primera colision inelastica
      

         t=0
      sumvsq=0

      do 250 j=1,npasos

         if (i.eq.1) then

            if ((j/(4*npart)).eq.contadorpeli) then
               print*,j,contadorpeli
c
c               call flush()
               do p=1,npart

c      write(90,*)p,rx(p),ry(p),vx(p),vy(p)
          write(100+contadorpeli,888)rx(p),rz(p),vx(p),vz(p)
c                close(1000+contadorpeli)
        enddo

        close(100+contadorpeli)
 888    format(4e20.10)
        
c             dencosx=0
c             dencosy=0
c             densinx=0
c             densiny=0
c             vxcosx=0
c             vxcosy=0
c             vxsinx=0
c             vxsiny=0
c             vycosx=0
c             vycosy=0
c             vysinx=0
c             vysiny=0
c             v2=0

c             do p=1,npart
c                dencosx=dencosx+cos(2*pi*rx(p)/longxy)
c                dencosy=dencosy+cos(2*pi*ry(p)/longxy)
c                densinx=densinx+sin(2*pi*rx(p)/longxy)
c                densiny=densiny+sin(2*pi*ry(p)/longxy)

c                vxcosx=vxcosx+vx(p)*cos(2*pi*rx(p)/longxy)
c                vxcosy=vxcosy+vx(p)*cos(2*pi*ry(p)/longxy)
c                vxsinx=vxsinx+vx(p)*sin(2*pi*rx(p)/longxy)
c                vxsiny=vxsiny+vx(p)*sin(2*pi*ry(p)/longxy)

c                vycosx=vycosx+vy(p)*cos(2*pi*rx(p)/longxy)
c                vycosy=vycosy+vy(p)*cos(2*pi*ry(p)/longxy)
c                vysinx=vysinx+vy(p)*sin(2*pi*rx(p)/longxy)
c                vysiny=vysiny+vy(p)*sin(2*pi*ry(p)/longxy)

c                v2=v2+(vx(p)*vx(p)+vy(p)*vy(p))/2.0d0
c             enddo
c             v2=v2/npart

c             write(56,765)j,dencosx,dencosy,densinx,densiny,
c     &vxcosx,vxcosy,vxsinx,vxsiny,
c     &vycosx,vycosy,vysinx,vysiny,v2,colisionpart
c 765         format(i20,13e20.10,i20)

             contadorpeli=contadorpeli+1
            endif

         endif
c         print*,j

c         call flush()

cccc  el j=30 tiene todavia energia cinetica uno y el J=31 menor que uno

c      if (j.gt.1*npart) then

          a=ai
c          w=w0
c      endif

cccc calculamos vx,vy,...y todas las magnitudes en tiempo real

      if (i.eq.1.and.(j/npart).eq.c4-1) then
c         print*,j/npart,c4-1,t
             tiempo(c4)=t
c      print*,j/npart,c4,tiempo(c4)
     
             c4=c4+1

      
      endif

cc      if (j/npart.ge.1.and.t.ge.tiempo(c1).and.c1.le.npasos/npart)
cc     & then

ccccccccccccccccc
cc     Meto la modificacion de tiempo real, a ver que pasa
cc
        if ((i.ne.1.and.j/npart.ge.0.and.t.ge.tiempo(c1).and. 
     &c1.le.npasos/npart).or.i.eq.1.and.t.eq.tiempo(c1)) 
     & then 
        
          sumvx=0
          sumvy=0
          sumvx2=0
          sumvy2=0
          sumvz2=0
          sumvx4=0
          sumvy4=0
          sumv4=0
          sumv6=0
          sumvxvy=0
          sumv2vx=0
          sumvx2vy2=0
          sumvxmenosvy2=0

c          sumvtcosx=0.0d0
c          sumvtsinx=0.0d0
c          sumvtcosy=0.0d0
c          sumvtsiny=0.0d0

          do p=1,npart
             sumvx=sumvx+vx(p)      
c             sumvy=sumvy+vy(p)
             sumvx2=sumvx2+vx(p)*vx(p)
c             sumvy2=sumvy2+vy(p)*vy(p)
             sumvz2=sumvz2+vz(p)*vz(p)
             sumvx4=sumvx4+vx(p)*vx(p)*vx(p)*vx(p)
c             sumvy4=sumvy4+vy(p)*vy(p)*vy(p)*vy(p)
             sumv4=sumv4+vx(p)*vx(p)*vx(p)*vx(p)
c             sumv6=sumv6+(vx(p)*vx(p)+vy(p)*vy(p))**3
c             sumvxvy=sumvxvy+vx(p)*vy(p)

c             sumv2vx=sumv2vx+(vx(p)*vx(p)+vy(p)*vy(p))*vx(p)

c             sumvx2vy2=sumvx2vy2+vx(p)*vx(p)*vy(p)*vy(p)
c             sumvxmenosvy2=sumvxmenosvy2+(vx(p)*vx(p)-vy(p)*vy(p))*
c     &       (vx(p)*vx(p)-vy(p)*vy(p))




          enddo

          sumvx=sumvx/npart
c          sumvy=sumvy/npart
          sumvx2=sumvx2/npart
c          sumvy2=sumvy2/npart
          sumvz2=sumvz2/npart

          sumvx4=sumvx4/npart
c          sumvy4=sumvy4/npart
          sumv4=sumv4/npart
          sumv6=sumv6/npart

cc          print*,datos(c1,1)
c          if (w0.ne.0)  datos(c1,1)=1.0d0/w0*dlog(t/t0)
          datos(c1,2)=t
          datos(c1,3)=j
          datos(c1,4)=sumvx
          datos(c1,5)=sumvy
          datos(c1,6)=sumvx2
          datos(c1,7)=sumvy2
          datos(c1,8)=(sumvx2+sumvy2)/2
          datos(c1,9)=sumvxvy*sumvxvy
          datos(c1,10)=sumvx2vy2
          datos(c1,11)=sumv4
          datos(c1,12)=(sumvx2+sumvy2)*(sumvx2+sumvy2)/4
          datos(c1,13)=presion3
          datos(c1,14)=sumv6
          datos(c1,15)=sumvz2
          datos(c1,16)=presion1
          datos(c1,17)=presion2
          datos(c1,18)=colisionpart
          
c          print*,j,(sumvx2+sumvy2)/2,sumvz2
c          call flush()
      
          c1=c1+1

c          colisionpart=0
c      call control(npart,npasos,densidad,box,rx,ry,vx,vy,i,j,
c     &flu,c3,umediaceldax,umediacelday,nmediacelda,nceldas,dEcmedia,
c     &Tempmedia,dTempmedia)

c      do kk=1,nceldas

c         do ll=1,nceldas

c            write(88,*)c3,kk,ll,nmediacelda(kk,ll,c3),
c     &umediaceldax(kk,ll,c3),
c     &umediacelday(kk,ll,c3)

c         enddo

c      enddo

c      c3=c3+1

          endif
     
      call tcolminimo(npart,tcol,partner,timbig,tcolmin,ii,jj)

c      print*,ii,jj
      if (ii.le.npart.and.jj.le.npart) then
      
      colisionpart=colisionpart+1

      endif
      if(tcolmin.lt.0) then
         print*,tcolmin,ii,jj,rz(ii),rz(jj),rx(ii),rx(jj)!,rx(ii),ry(ii),vx(ii),vy(ii)
c         call flush()
      endif

ccccc       DETERMINACION COLAPSO INELASTICO
ccccc      Lo que vamos a ver es si la distancia entre las
ccccc      particulas que van a colisionar es menor que 10e-15

      if(ii.le.npart.and.jj.le.npart) then
       colp=(rx(ii)-rx(jj))*(rx(ii)-rx(jj))+
     &(rz(ii)-rz(jj))*(rz(ii)-rz(jj))
  
      If (colp.le.(1+1.0d-12)) then
         print*,'Cuidado,puede haber colapso inelastico',colp,j,ii,jj,
     &tcolmin,longx,rx(ii),
     &rx(jj),rz(ii),rz(jj)
      endif

      endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
cccc  calculamos el tiempo, tomando como origen de tau la ultima collision elas

c      if (j.gt.1*npart) then
              t=t+tcolmin 
c      endif

      call movimiento(npart,longx,longz,
     &tcol,tcolmin,rx,rz,vx,vz,ii,jj,virial,
     &           gmas,a,t,vibra,acumpresion1,acumpresion2,acumpresion3)
c      ,
c     &rparedx)
      
      presion1=presion1+acumpresion1
      presion2=presion2+acumpresion2
      presion3=presion3+acumpresion3

c      print*,ii,jj
c      call flush()

      call nuevostcol(npart,longx,longz,
     &timbig,rx,rz,vx,vz,partner,tcol,
     &ii,jj)

cccc  Restamos (momento total)/(npart), a todas las velocidades,
cccc  para que el momento total sea siempre cero
cc      if (j.gt.1*npart+3) then
cc      sumvx=0
cc      sumvy=0

cc      do p=1,npart
cc         sumvx=sumvx+vx(p)
cc         sumvy=sumvy+vy(p)
cc      enddo

cc      do p=1,npart
cc         vx(p)=vx(p)-sumvx/npart
cc         vy(p)=vy(p)-sumvy/npart
cc      enddo
cc      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccc  calculamos la fvx y gdr 10 veces en todo el programa

      if (j/(5*ndist).eq.c6) then

c          call fdistribucion(npart,maxf,intervalof,vx,vy,vz,fvx,fvy,fvz)
     
c          do p=-maxf,maxf
c             mfvx(c6,p)=fvx(p)
c             mfvy(c6,p)=fvy(p)
c             mfvz(c6,p)=fvz(p)
c          enddo
          
c          mfvx(c6,maxf+1)=t
c          mfvy(c6,maxf+1)=t
c          mfvz(c6,maxf+1)=t
c          ec=0
c          ecz=0

c          do p=1,npart
c             ec=ec+(vx(p)*vx(p)+vy(p)*vy(p))/2
c             ecz=ecz+(vz(p)*vz(p))/2
c          enddo

c          ec=ec/npart
   
c          mfvx(c6,maxf+2)=sqrt(2*ec)   
c          mfvy(c6,maxf+2)=sqrt(2*ec)
c          mfvz(c6,maxf+2)=sqrt(2*ecz)   
c          c6=c6+1

          contador2(i)=contador2(i)+1

      endif          
     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccc  calculamos flu  cada 5 colisiones  por particula

      if (j/(5*npart).eq.c3) then
         print*,j,j/npart
c         call flush() 

ccccccc         call control(npart,npasos,densidad,longxy,rx,ry,vx,vy,i,j,
ccccccc     &flu,c3,umediaceldax,umediacelday,nmediacelda,nceldas,dEcmedia,
ccccccc     &Tempmedia,dTempmedia)


c         if ((c3/5)*5.eq.c3) then
c         print*,c3
c      do p=1,npart

c      write(90,*)p,rx(p),ry(p),rz(p),vx(p),vy(p),vz(p)
c      write(1000+c3/5,*)rx(p),ry(p)
c      enddo
c      endif
      
c         print*,c3
cccccc         do kk=1,nceldas

cccccc         do ll=1,nceldas

cccccc            write(88,*)c3,kk,ll,nmediacelda(kk,ll,c3),
cccccc     &umediaceldax(kk,ll,c3),
cccccc     &umediacelday(kk,ll,c3)

cccccc         enddo

cccccc      enddo

         c3=c3+1

c      print*,i,contador2(i)
      endif

      if (i.eq.1) then

      if (j/npart.eq.100*contnuevo) then
         do p=1,npart
            
            write(97,*)j,p,rx(p),rz(p),vx(p),vz(p)
            
         enddo
         contnuevo=contnuevo+1
      endif

      endif
cccccccccccccccccccccccccccccccccccccccccccccccccc
     
 250  continue
 
      

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc      Escribimos las velocidades en cada celda cada 5 pasos
ccc  por particula para la primera trayectoria a fin de ver si el
ccc  sistema es estable
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


 
cccc  promediamos en trayectorias

      do p=1,npasos/npart
         do q=1,18
            sumdatos(p,q)=sumdatos(p,q)+datos(p,q)

         enddo
      enddo

      write(*,*)datos(1,6),datos(1,7),datos(1,15)
      write(*,*)sumdatos(1,6),sumdatos(1,7),sumdatos(1,15)

c      do p=1,npasos/(5*ndist)
c         write(*,*)p,mfvx(p,maxf+2)
c      enddo
      
ccc      do p=1,npasos/(5*ndist)
ccc         do q=-maxf,maxf+2
ccc            summfvx(p,q)=summfvx(p,q)+mfvx(p,q)
ccc            summfvy(p,q)=summfvy(p,q)+mfvy(p,q)
ccc            summfvz(p,q)=summfvz(p,q)+mfvz(p,q)
ccc         enddo
c         write(*,*)p,summfvx(p,maxf+2)
ccc      enddo
      


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc  GUARDAMOS LOS MOMENTOS Y FLUC cada 5 trayectorias  ccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc




 2000 continue

cccc  Escribimos en los ficheros

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      escribimos (E(t)-<E(t)>)/<E(t)> para ranf trayectorias
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do p=1,npasos/npart
       
          write(10,190)sumdatos(p,1)/ranf,
     &               sumdatos(p,2)/ranf,
     &               sumdatos(p,3)/ranf,
     &               sumdatos(p,4)/ranf,
     &               sumdatos(p,5)/ranf,
     &               sumdatos(p,6)/ranf,
     &               sumdatos(p,7)/ranf,
     &               sumdatos(p,8)/ranf,
     &               sumdatos(p,9)/ranf,
     &               sumdatos(p,10)/ranf,
     &              (sumdatos(p,11)/ranf)/(2*(2*sumdatos(p,8)/ranf)**2),
     &              (sumdatos(p,12)/ranf)/(sumdatos(p,8)/ranf)**2-1,
     &               sumdatos(p,13),
     &               (sumdatos(p,14)/ranf)/(2*sumdatos(p,8)/ranf)**3,
     &     (sumdatos(p,15)/ranf),
     &     (sumdatos(p,16)/ranf),    
     &     (sumdatos(p,17)/ranf), 
     &     (sumdatos(p,18)/ranf),
     &     (sumjavier2(p)/ranf)/(sumdatos(p,8)/ranf)**2,
     &     (sumjavier3(p)/ranf)/(sumdatos(p,8)/ranf)**3,
     &     (sumjavier4(p)/ranf)/(sumdatos(p,8)/ranf)**4,    
     &     (sumjavier5(p)/ranf)/(sumdatos(p,8)/ranf)**5, 
     &     (sumjavier6(p)/ranf)/(sumdatos(p,8)/ranf)**6,
     &      sumjavier2(p)/ranf,
     &      sumjavier3(p)/ranf,
     &      sumjavier4(p)/ranf,
     &      sumjavier5(p)/ranf,
     &      sumjavier6(p)/ranf
    
      enddo
   
 190  format(28e20.10)


       c2max=contador2(1)
       do p=1,ranf
          print*,contador2(p)
c          call flush()
          c2max=min(contador2(p),c2max)

       enddo

       print*,c2max
c       call flush()

ccccc      do p=1,c2max
c        write(*,*)summfvx(p,maxf+2)
c        write(20,*)
cccc         do q=-maxf,maxf
cccc            write(20,555)(q-0.5)*intervalof,summfvx(p,q)/
cccc     &                 (npart*ranf*intervalof),summfvx(p,maxf+1)/ranf,
cccc     &summfvy(p,q)/
cccc     &                 (npart*ranf*intervalof),summfvy(p,maxf+1)/ranf,
cccc     &summfvz(p,q)/
cccc     &           (npart*ranf*intervalof),summfvz(p,maxf+1)/ranf

cccc            write(58,*)(q-0.5)*intervalof/(summfvx(p,maxf+2)/ranf),
cccc     &summfvx(p,q)*summfvx(p,maxf+2)/(npart*ranf*intervalof*ranf),
cccc     &summfvx(p,maxf+1)/ranf

cccc            write(59,*)(q-0.5)*intervalof/(summfvy(p,maxf+2)/ranf),
cccc     &summfvy(p,q)*summfvy(p,maxf+2)/(npart*ranf*intervalof*ranf),
cccc     &summfvy(p,maxf+1)/ranf

cccc            write(60,*)(q-0.5)*intervalof/(summfvz(p,maxf+2)/ranf),
cccc     &summfvz(p,q)*summfvz(p,maxf+2)/(npart*ranf*intervalof*ranf),
cccc     &summfvz(p,maxf+1)/ranf
cccc         enddo
cccc      enddo

 555  format(7e20.10)


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      


      end







      subroutine correlacionespacial(npart,box,rx,ry,vx,
     &vy,nceldas,corr,corrymedio,ux,densi,temperature,corry2,tempetotal)

      implicit none
      integer nceldas,npart,celdilla(npart)

      double precision ux(nceldas),uy(nceldas),rx(npart),ry(npart),
     &vx(npart),vy(npart),box,corr(0:nceldas-1),corrymedio(nceldas),
     &densi(nceldas),corry2(nceldas),tempetotal

      integer k,p,l,cont(0:nceldas-1),m

      double precision temperature(nceldas)


CCCCCCCC   Inicializamos la matriz ux, que es la velocidad total de una capa
      do l=0,nceldas-1

         corr(l)=0.0d0
         cont(l)=0
        
      enddo 


      do k=1,nceldas
         densi(k)=0.0d0
         ux(k)=0.0d0
         uy(k)=0.0d0
         temperature(k)=0.0d0
         corrymedio(k)=0.0d0
         corry2(k)=0.0d0
      enddo


ccccc      Asignamos las posiciones de las particulas al indice k de cada capa


         do p=1,npart

ccccc      Distinguimos si el numero de celdas es par o impar cccccccc

ccccccccccccccccc  PAR  ccccccccccccccccccccccccccccccccccc

         if (nceldas/2*2.eq.nceldas) then
c            if (rx(p).gt.0) then
c               k=int(rx(p)*nceldas/box)+nceldas/2+1
c            else
c               k=int(rx(p)*nceldas/box)+nceldas/2
c            endif


            if (ry(p).eq.box/2.0d0) then
               k=nceldas
c               print*,'limite superior'
c               call flush()

            else if (ry(p).eq.-box/2.0d0) then
               k=1
c               print*,'limite inferior'
c               call flush()
               
            else if (ry(p).gt.0) then
               k=int(ry(p)*nceldas/box)+nceldas/2+1
            else
               k=int(ry(p)*nceldas/box)+nceldas/2
            endif 

         else  


cccccccccccccccc  IMPAR  cccccccccccccccccccccccccccccccccccccc

c            if (rx(p).gt.-box/nceldas) then
c               k=int(rx(p)*nceldas/box+0.5)+int(nceldas/2)+1
c            else
c               k=int(rx(p)*nceldas/box+0.5)+int(nceldas/2)
c            endif
       
           
            if (ry(p).eq.box/2.0d0) then
               k=nceldas

c               print*,'limite superior'
c               call flush()

            else if (ry(p).eq.-box/2.0d0) then
               k=1
c               print*,'limite inferior'
c               call flush()


            else if (ry(p).gt.-box/(2*nceldas)) then
               k=int(ry(p)*nceldas/box+0.5)+int(nceldas/2)+1
            else
               k=int(ry(p)*nceldas/box+0.5)+int(nceldas/2)
            endif 
         endif
cc         print*,k,vx(p)
cc         call flush()
         
c         k=int(k)
         if (k.gt.npart.or.k.lt.1) then

            print*,'error al determinar k',k,ry(p),box/2,ry(p)+box/2

c            call flush()

         endif


            ux(k)=ux(k)+vx(p)
            uy(k)=uy(k)+vy(p)
            densi(k)=densi(k)+1
            celdilla(p)=k

         enddo
         do p=1,nceldas

            ux(p)=ux(p)/densi(p)
            uy(p)=uy(p)/densi(p)
            
         enddo

         do p=1,npart
            m=celdilla(p)
            temperature(m)=temperature(m)+(vx(p)-ux(m))*
     &           (vx(p)-ux(m))/2+
     &           (vy(p)-uy(m))*(vy(p)-uy(m))/2

         enddo
         do p=1,nceldas
            temperature(p)=temperature(p)/densi(p)
         enddo
cccccccccc       
         tempetotal=0.0d0

      do k=1,nceldas

       tempetotal=tempetotal+temperature(k)

      enddo

      tempetotal=tempetotal/nceldas




cccccccccc       Calculamos la funcion de correlacion

      do k=1,nceldas
         do l=0,nceldas-1

            corr(l)=corr(l)+ux(k)*ux(k+l)
            cont(l)=cont(l)+1

         enddo

      enddo


      do k=1,nceldas

         corrymedio(k)=corrymedio(k)+ux(nceldas/2)*ux(k)
         corry2(k)=corry2(k)+ux(2)*ux(k)

      enddo

      return



      end








      subroutine control(npart,npasos,densidad,longxy,rx,ry,vx,vy,i,j,
     &flu,c3,umediaceldax,umediacelday,nmediacelda,nceldas,dEcmedia,
     &Tempmedia,dTempmedia)

      implicit none
      integer npart,nceldas,npasos,i,j,p,k,l,ll,kk,c3
c      parameter(nceldas=6)
      double precision rx(npart),ry(npart),vx(npart),vy(npart),longxy,
     &n(nceldas,nceldas),nf,densidad,vceldax(nceldas,nceldas),
     &vcelday(nceldas,nceldas),dEc(nceldas,nceldas),
     &sumvcelda2,sumvsq,flu(npasos/(5*npart),2),
     &umediaceldax(nceldas,nceldas,npasos/(5*npart)),
     &umediacelday(nceldas,nceldas,npasos/(5*npart)),
     &nmediacelda(nceldas,nceldas,npasos/(5*npart)),
     &dEcmedia(nceldas,nceldas,npasos/(5*npart)),Ec,
     &Tempmedia(nceldas,nceldas,npasos/(5*npart)),
     &dTempmedia(nceldas,nceldas,npasos/(5*npart))

            Ec=0.0d0
      do kk=1,nceldas
         do ll=1,nceldas
            n(ll,kk)=0
            vceldax(ll,kk)=0
            vcelday(ll,kk)=0
            dEc(ll,kk)=0
         enddo
      enddo

     

      do p=1,npart

cccc  Distinguimos si el numero de celdas es par o impar cccccccc
ccccccccccccccccc  PAR  ccccccccccccccccccccccccccccccccccc

         if (nceldas/2*2.eq.nceldas) then
            if (rx(p).gt.0) then
               k=int(rx(p)*nceldas/longxy)+nceldas/2+1
            else
               k=int(rx(p)*nceldas/longxy)+nceldas/2
            endif
         
            if (ry(p).gt.0) then
               l=int(ry(p)*nceldas/longxy)+nceldas/2+1
            else
               l=int(ry(p)*nceldas/longxy)+nceldas/2
            endif 

         else  


cccccccccccccccc  IMPAR  cccccccccccccccccccccccccccccccccccccc

            if (rx(p).gt.-longxy/(2*nceldas)) then
               k=int(rx(p)*nceldas/longxy+0.5)+int(nceldas/2)+1
            else
               k=int(rx(p)*nceldas/longxy+0.5)+int(nceldas/2)
            endif
         
            if (ry(p).gt.-longxy/(2*nceldas)) then
               l=int(ry(p)*nceldas/longxy+0.5)+int(nceldas/2)+1
            else
               l=int(ry(p)*nceldas/longxy+0.5)+int(nceldas/2)
            endif 

          endif

         n(l,k)=n(l,k)+1
         vceldax(l,k)=vceldax(l,k)+vx(p)
         vcelday(l,k)=vcelday(l,k)+vy(p)
         dEc(l,k)=dEc(l,k)+1.0d0/2.0d0*(vx(p)*vx(p)+vy(p)*vy(p))
      enddo
      
      nf=0
      sumvcelda2=0
      do k=1,nceldas
         do l=1,nceldas
            if (n(l,k).eq.0) then
             print*,'mierda,n es cero',l,k,n(l,k)
            endif
            if (n(l,k).ne.0) then

cc  MAL             nf=nf+(n(l,k)/(longxy/nceldas)**2-densidad)**2
cc  MAL           vceldax(l,k)=vceldax(l,k)/n(l,k)
cc  MAL            vcelday(l,k)=vcelday(l,k)/n(l,k)
cc  MAL             dEc(l,k)=dEc(l,k)/n(l,k)
ccccccccccccccc MODIFICACION VORTICES 


            if (i.eq.1) then

             do p=1,npart
                Ec=Ec+1.0d0/2.0d0*(vx(p)*vx(p)+vy(p)*vy(p))
             enddo
                Ec=Ec/npart
             umediaceldax(l,k,c3)=vceldax(l,k)         
             umediacelday(l,k,c3)=vcelday(l,k)
             nmediacelda(l,k,c3)=n(l,k)  
             dEcmedia(l,k,c3)=dEc(l,k)-Ec
             Tempmedia(l,k,c3)=dEc(l,k)-
     &                 1.0d0/2.0d0*(vceldax(l,k)*vceldax(l,k)+
     &                                       vcelday(l,k)*vcelday(l,k))
             dTempmedia(l,k,c3)=Tempmedia(l,k,c3)-Ec

            endif     

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            endif

            sumvcelda2=sumvcelda2+vceldax(l,k)**2+vcelday(l,k)**2
         enddo
      enddo

      nf=nf/(nceldas*nceldas)
      sumvcelda2=sumvcelda2/(nceldas*nceldas)
      
      sumvsq=0
      do p=1,npart
         sumvsq=sumvsq+vx(p)*vx(p)+vy(p)*vy(p)
      enddo

         sumvsq=sumvsq/(2*npart)

      flu(c3,1)=nf
      flu(c3,2)=sumvcelda2/sumvsq

      return
      end

      
      subroutine f2distribution(npart,maxf,intervalof,f2vx,vx)
      implicit none
      integer npart, maxf,a,b,m,n
      double precision f2vx(-maxf:maxf,-maxf:maxf),vx(npart),intervalof

      do a=-maxf,maxf
         do b=-maxf,maxf
            f2vx(a,b)=0.0d0
         enddo
      enddo
      
      do 1000 a=1,npart-1
         do 2000 b=a+1,npart

         if (vx(a).ge.0) then
            m=int(vx(a)/intervalof)+1
         else 
            m=int(vx(a)/intervalof)
         endif
         if (vx(b).ge.0) then
            n=int(vx(b)/intervalof)+1
         else 
            n=int(vx(b)/intervalof)
         endif
         if ((m.ge.-maxf).and.(m.le.maxf).and.(n.ge.-maxf).and.
     &       (n.le.maxf)) then
         f2vx(m,n)=f2vx(m,n)+1
         f2vx(n,m)=f2vx(n,m)+1
         endif

         
 2000    continue
 1000 continue   

      




      return
      end


      subroutine gdr(npart,box,densidad,maxg,intervalog,rx,ry,g)

      implicit none
      integer npart,r,i,j,maxg
      double precision box,densidad,intervalog,rx(npart),ry(npart),rxij,
     &ryij,rijsq,rij,g(maxg),rinf,rsup,nideal

      do i=1,maxg
         g(i)=0
      enddo
      do 100 i=1,npart-1
  
         do 99 j=i+1,npart
            rxij=rx(i)-rx(j)
            ryij=ry(i)-ry(j)
            rxij=rxij-box*anint(rxij/box)
c            ryij=ryij-box*anint(ryij/box)
       
            rijsq=rxij**2+ryij**2
          
            rij=sqrt(rijsq)
cccccccccccccccccccccccccccccccccccc
            if (rij.lt.0.9) then
               print*,'mierda',rij,i,j
c               call flush()
            endif
cccccccccccccccccccccccccccccccccccc
          
            r=int(rij/intervalog)+1
            
            if (r.le.maxg) then
            
               g(r)=g(r)+2
       
            endif
 99         continue
 100     continue
cccccccccccccccccccccccccccccccccccccccc
ccccc dividamos por nideal y por npart
cccccccccccccccccccccccccccccccccccccccc
      do i=1,maxg
         rinf=dble(i-1)*intervalog
         rsup=dble(i)*intervalog
         nideal=3.141592*(rsup*rsup-rinf*rinf)*densidad
         g(i)=g(i)/(npart*nideal)
      enddo         
      return
      end   

      subroutine fdistribucion(npart,maxf,intervalof,
     &vx,vy,vz,fvx,fvy,fvz)
      
      implicit none
      integer npart,maxf,j,k,m,p,l,r
      double precision intervalof, vx(npart),vy(npart),vz(npart),
     &fvx(-maxf:maxf),fvy(-maxf:maxf),fvz(-maxf:maxf)
     

         do l=-maxf,maxf
         fvx(l)=0
         fvy(l)=0
         fvz(l)=0
         enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         
         do 500 j=1,npart
    
         if (vx(j).gt.0) then
         k=int(vx(j)/intervalof)+1
         else
         k=int(vx(j)/intervalof)
         endif
         if (vy(j).gt.0) then
         m=int(vy(j)/intervalof)+1
         else
         m=int(vy(j)/intervalof)
         endif
         if (vz(j).gt.0) then
         r=int(vz(j)/intervalof)+1
         else
         r=int(vz(j)/intervalof)
         endif
         if ((k.ge.-maxf).and.(k.le.maxf)) then
         fvx(k)=fvx(k)+1
         endif
         if ((m.ge.-maxf).and.(m.le.maxf)) then
         fvy(m)=fvy(m)+1
         endif
         if ((r.ge.-maxf).and.(r.le.maxf)) then
         fvz(r)=fvz(r)+1
         endif
 500     continue

         return
         end





      subroutine nuevostcol(npart,longx,longz,
     &timbig,rx,rz,vx,vz,partner,tcol,
     & ii,jj)

ccccc calcula los tiempos de colision y con que particula colisiona cada una 
ccccc ????????

      integer i,j,npart,ii,jj,partner(npart+2),p
      double precision timbig,rx(npart),rz(npart),
     &vx(npart),vz(npart),
     &rxij,rzij,
     &rijsq,vxij,vzij,vijsq,bij,disc,tij,tcol(npart+2)

      double precision longx,longz,rparedx

      
     
      do i=npart+1,npart+2

               if ((i.eq.ii).or.(i.eq.jj).or.(partner(i).eq.ii).or.
     &(partner(i).eq.jj)) then

                  tcol(i)=timbig

                  endif
               enddo

      do 100 i=1,npart
         if ((i.eq.ii).or.(i.eq.jj).or.(partner(i).eq.ii).or.
     &(partner(i).eq.jj)) then
         tcol(i)=timbig

         do 200 j=1,npart
            if (i.ne.j)then
            rxij=rx(i)-rx(j) 
            
            rzij=rz(i)-rz(j)
            rxij=rxij-longx*anint(rxij/longx)
           

            rijsq=rxij**2.+rzij**2.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccc     DETECCION DE COLAPSO INELASTICO   ccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c            if ((i.eq.ii.and.j.ne.jj).or.(i.eq.jj.and.j.ne.ii)
c     &              .and.rijsq.le.(1.0d0+2.0d-15))then
c                print*,'PARECE QUE HAY COLAPSO INELASTICO MARIBEL'
c                
c            endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


            vxij=vx(i)-vx(j)
            
            vzij=vz(i)-vz(j)
            vijsq=vxij**2.+vyij**2.+vzij**2.
            bij=rxij*vxij+rzij*vzij

ccccc si es posible la colision entre las particulas (i,j) (bij<0 y disc>0) se ccccc calcula el tiempo de colision y se actualiza la matriz tcol()

            if (bij.lt.0) then
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccc     DETECCION DE COLAPSO INELASTICO   ccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c            if ((i.eq.ii).or.(i.eq.jj).and.rijsq.le.(1.0d0+2.0d-15))then
c                print*,'PARECE QUE HAY COLAPSO INELASTICO PABLO'
c            endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
ccccc calculamos el discriminabte

                disc=bij**2-vijsq*(rijsq-1) 
                if (disc.gt.0) then
                   tij=(-bij-disc**(1.0d0/2.0d0))/vijsq
                      if(tij.lt.tcol(i)) then
                        tcol(i)=tij
                        partner(i)=j
                      endif
                      if (tij.lt.tcol(j))then
                         tcol(j)=tij
                         partner(j)=i
                      endif
                endif
               endif
              endif 
 200         continue
       

      if (tcol(i).lt.0) then

         write(77,*)i,tcol(i)

         endif

      endif

 100         continue

ccc   Tiempos de colision con las paredes
 
      do i=npart+1,npart+2
c
                  
         do p=1,npart

cc            if (i.eq.npart+3) then

cc               if (rx(p).lt.rparedx.and.vx(p).gt.0) then

cc                  tij=(-rx(p)+rparedx)/vx(p)
c                 print*,p,tij
                
cc                  if (tij.lt.tcol(p)) then
                     
cc                     tcol(p)=tij
cc                     partner(p)=i

cc                  endif
cc                  if (tcol(i).gt.tij) then
cc                     tcol(i)=tij
cc                     partner(i)=p
cc                  endif

cc            endif

cc            if (rx(p).gt.rparedx.and.vx(p).lt.0) then

cc                  tij=(-rx(p)+rparedx)/vx(p)
c                 print*,p,tij
                
cc                  if (tij.lt.tcol(p)) then
                     
cc                     tcol(p)=tij
cc                     partner(p)=i

cc                  endif
cc                  if (tcol(i).gt.tij) then
cc                     tcol(i)=tij
cc                     partner(i)=p
cc                  endif

cc            endif


cc         endif
            
            if(i.eq.npart+2) then

               if (vz(p).lt.0) then

                  tij=(-(longz-1.0d0)/2.0d0-rz(p))/vz(p)
                
                  if (tij.lt.tcol(p)) then

                     tcol(p)=tij
                     partner(p)=i

               
                  endif
                  if (tcol(i).gt.tij) then
                     tcol(i)=tij
                     partner(i)=p
                  endif
               endif
            endif


            if(i.eq.npart+1) then
               if (vz(p).gt.0) then
                  
                  tij=((longz-1.0d0)/2.0d0-rz(p))/vz(p)

                  if (tij.lt.tcol(p)) then

                     tcol(p)=tij
                     partner(p)=i

                  endif
  
                  if (tcol(i).gt.tij) then
                     tcol(i)=tij
                     partner(i)=p
                  endif
               endif
            endif
            
         enddo 
                  
      enddo  
c      print*,tcol(254)
       return
       end  



      
      subroutine movimiento(npart,longx,longz,
     &tcol,tcolmin,rx,rz,vx,vz,ii,jj,
     &     virial,gmas,a,t,vibra,acumpresion1,acumpresion2,acumpresion3)
C,rparedx      )

ccccc  esta subrutina mueve todas las particulas tcolmin, aplica las
ccccc reglas de colision ,actualiza tcol() y calcula el virial
      implicit none
      integer i,npart,ii,jj,p
      double precision longx,longz,rx(npart),rz(npart)
      double precision vx(npart),vz(npart),
     &tcolmin,tcol(npart+2),rxiijj,rziijj,
     &vxiijj,vziijj,biijj,virial,
     &gmas,a,t,colp,vibra,acumpresion1,acumpresion2,acumpresion3
     inicializacion2d
c      double precision rparedx

      acumpresion1=0.0d0
      acumpresion2=0.0d0
      acumpresion3=0.0d0


ccccc movemos las particulas tcolmin,aplicamos condiciones de contorno 
ccccc y actualizamos tcol()
      
      Do 100 i=1,npart

         rx(i)=rx(i)+vx(i)*tcolmin
         
         rz(i)=rz(i)+vz(i)*tcolmin
         rx(i)=rx(i)-longx*anint(rx(i)/longx)
         
         
ccc Modificacion USF

c         ry(i)=ry(i)-box*anint(ry(i)/box)

ccc  Fin modificacion USF

         tcol(i)=tcol(i)-tcolmin

 100  continue

cc        if (ii.eq.npart+3) then

cc           rx(jj)=rparedx
c            print*,jj,rx(jj)


cc         else if (jj.eq.npart+3) then

cc            rx(ii)=rparedx
c            print*,ii,rx(ii)

c         endif

      do i=npart+1,npart+2
         tcol(i)=tcol(i)-tcolmin

      enddo

ccc  Modificacion USF
      if (ii.le.npart.and.jj.le.npart) then

     
ccc Fin USF

      rxiijj=rx(ii)-rx(jj)      
      rziijj=rz(ii)-rz(jj)
      rxiijj=rxiijj-longx*anint(rxiijj/longx)     
      vxiijj=vx(ii)-vx(jj)    
      vziijj=vz(ii)-vz(jj)
      biijj=rxiijj*vxiijj+rziijj*vziijj
ccccc ahora aplicamos las reglas de colision
      vx(ii)=vx(ii)-(1+a)/2*biijj*rxiijj 
      vz(ii)=vz(ii)-(1+a)/2*biijj*rziijj
      vx(jj)=vx(jj)+(1+a)/2*biijj*rxiijj 
      vz(jj)=vz(jj)+(1+a)/2*biijj*rziijj
ccccc calculamos el virial y la gsigma+
      virial=-(1+a)/2*biijj
      gmas=-1/biijj
ccc  Modificacion USF
      endif

ccc  Colisiones con las paredes 

      if (ii.eq.npart+1) then
        
         acumpresion1=2*vz(jj)
c         print*,acumpresion1

         vz(jj)=-vz(jj)
         
      endif

      if (ii.eq.npart+2) then
         
         acumpresion2=-2*vz(jj)+2*vibra !Sera asi?
         vz(jj)=-vz(jj)+2*vibra
       
      endif

cc      if (ii.eq.npart+3) then
         
cc         acumpresion3=Abs(vx(jj)) 
c         print*,acumpresion3
       
c      endif

      if (jj.eq.npart+1) then

         acumpresion1=2*vz(ii)
c         print*,acumpresion1

         vz(ii)=-vz(ii)
         
      endif
      

      if (jj.eq.npart+2) then

         acumpresion2=-2*vz(ii)+2*vibra !sera asi?
         vz(ii)=-vz(ii)+2*vibra
       
      endif

cc      if (jj.eq.npart+3) then

cc         acumpresion3=Abs(vx(ii)) 
c         print*,ii,acumpresion3,rx(ii)
         
cc      endif
      return
      end

      subroutine inicializacion(npart,longxy,temp,rx,ry,rz,vx,vy,vz,pi)

ccccc  En esta subrutina vamos a asignar posiciones y velocidades a todas
ccccc las particulas del sistema de forma que consigamos la temperatura T
ccccc que queramos y que la velocidad total sea nula


      implicit none
      integer i,j,k,npart,number,nplace,Ijkl,ii
      double precision rx(npart),ry(npart),rz(npart),
     & size,vx(npart),vy(npart),vz(npart),
     & sumvx,sumvy,sumvz,sumvsq,temp,Ran_Uniform,alea

      double precision longxy

      double precision x1,x2,y1,y2,pi
cccccc inicialicemos los acumuladoress 
    
      nplace=0
      sumvx=0
      sumvy=0
      sumvz=0
      sumvsq=0

            
cccccc introducimos las variables caracteristicas del sistema      
      
c     npart=100
c      box=11.0d0
      number=idint(dble(npart)**(1.0d0/2.0d0)+1)
      size=longxy/dble(number)
c      temp=1
      
      
cccccc establecemos las posiciones y velocidades de las particulas  
          
       do 300 i=1,number
        do 200 j=1,number
           
           nplace=nplace+1
           if(nplace.le.npart) then
            rx(nplace)=(i-0.5)*size-longxy/2.0d0
            ry(nplace)=(j-0.5)*size-longxy/2.0d0
            rz(nplace)=0.0d0

           endif
 
  200    continue
  300  continue
 
	do ii=1,npart,2
 2012      continue
		x1=Ran_Uniform()
		y1=Ran_Uniform()
                x2=Ran_Uniform()
                y2=Ran_Uniform()

                if((x1*y1).lt.1.0e-15) goto 2012
                if((x2*y2).lt.1.0e-15) goto 2012

		vx(ii)=sqrt(-2.*log(x1))*cos(2.*pi*y1)
                vy(ii)=sqrt(-2.*log(x2))*cos(2.*pi*y2)
		vx(ii+1)=sqrt(-2.*log(x1))*sin(2.*pi*y1)
                vy(ii+1)=sqrt(-2.*log(x2))*sin(2.*pi*y2)
                vz(ii)=Ran_Uniform()-0.5d0
                vz(ii+1)=Ran_Uniform()-0.5d0
                sumvx=sumvx+vx(ii)+vx(ii+1)
                sumvy=sumvy+vy(ii)+vy(ii+1)
                sumvz=sumvz+vz(ii)+vz(ii+1)
	enddo       
          

         
cccccc corregimos las velocidades para que el momento total del sistema sea cero

       do 100 i=1,npart
         vx(i)=vx(i)-sumvx/dble(npart)
         vy(i)=vy(i)-sumvy/dble(npart)
         vz(i)=vz(i)-sumvz/dble(npart)
         sumvsq=sumvsq+vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i)
  100  continue
  
cccccc escalamos las velocidades para que la temperatura sea la deseada

       do 400 i=1,npart
         vx(i)=(3.0d0*dble(npart)*temp/sumvsq)**(1.0d0/2.0d0)*vx(i) 
         vy(i)=(3.0d0*dble(npart)*temp/sumvsq)**(1.0d0/2.0d0)*vy(i)
         vz(i)=(3.0d0*dble(npart)*temp/sumvsq)**(1.0d0/2.0d0)*vz(i)
     
         
 400   continue

       return      
       end

      subroutine inicializacion0(npart,longx,temp,tempz,
     &  rx,rz,vx,vz,pi)

ccccc  En esta subrutina vamos a asignar posiciones y velocidades a todas
ccccc las particulas del sistema de forma que consigamos la temperatura T
ccccc que queramos y que la velocidad total sea nula


      implicit none
      integer i,j,k,npart,number,nplace,Ijkl,ii
      double precision rx(npart),rz(npart),
     & size,vx(npart),vz(npart),
     & sumvx,sumvz,sumvsq,sumvsqz,temp,tempz,Ran_Uniform,alea

      double precision longx

      double precision x1,x2,x3,y1,y2,y3,pi
cccccc inicialicemos los acumuladoress 
    
      nplace=0
      sumvx=0
      sumvz=0
      sumvsq=0
      sumvsqz=0
            
cccccc introducimos las variables caracteristicas del sistema      
      
c     npart=100
c      box=11.0d0
      number=npart
      size=longx/dble(number)
c      temp=1
      
      
cccccc establecemos las posiciones y velocidades de las particulas  
          
       do 300 i=1,number
             
           nplace=nplace+1
           if(nplace.le.npart) then
            rx(nplace)=(i-0.5)*size-longx/2.0d0
            rz(nplace)=0.0d0

           endif
 

  300  continue
 
	do ii=1,npart,2
 2012      continue
		x1=Ran_Uniform()
		y1=Ran_Uniform()
                
               
                x3=Ran_Uniform()
                y3=Ran_Uniform()

                if((x1*y1).lt.1.0e-15) goto 2012

                if((x3*y3).lt.1.0e-15) goto 2012

		vx(ii)=sqrt(-2.*log(x1))*cos(2.*pi*y1)
                
                vZ(ii)=sqrt(-2.*log(x3))*cos(2.*pi*y3)
		vx(ii+1)=sqrt(-2.*log(x1))*sin(2.*pi*y1)
                
                vZ(ii+1)=sqrt(-2.*log(x3))*sin(2.*pi*y3)

                sumvx=sumvx+vx(ii)+vx(ii+1)
               
                sumvz=sumvz+vz(ii)+vz(ii+1)
	enddo       
          

         
cccccc corregimos las velocidades para que el momento total del sistema sea cero

       do 100 i=1,npart
         vx(i)=vx(i)-sumvx/dble(npart)
         
         vz(i)=vz(i)-sumvz/dble(npart)
         sumvsq=sumvsq+vx(i)*vx(i)
         sumvsqz=sumvsqz+vz(i)*vz(i)
  100  continue
  
cccccc escalamos las velocidades para que la temperatura sea la deseada

       do 400 i=1,npart
         vx(i)=(1.0d0*dble(npart)*temp/sumvsq)**(1.0d0/2.0d0)*vx(i) 
         
         vz(i)=(1.0d0*dble(npart)*tempz/sumvsqz)**(1.0d0/2.0d0)*vz(i)
     
         
 400   continue

       return      
       end

   
      subroutine inicializacion4(npart,longxy,temp,rx,ry,rz,vx,vy,vz)

ccccc  En esta subrutina vamos a asignar posiciones y velocidades a todas
ccccc las particulas del sistema de forma que consigamos la temperatura T
ccccc que queramos y que la velocidad total sea nula


      implicit none
      integer i,j,k,npart,number,nplace,Ijkl
      double precision rx(npart),ry(npart),rz(npart),
     & size,vx(npart),vy(npart),vz(npart),
     & sumvx,sumvy,sumvz,sumvsq ,temp,Ran_Uniform

      double precision longxy
cccccc inicialicemos los acumuladoress 
    
      nplace=0
      sumvx=0
      sumvy=0
      sumvz=0
      sumvsq=0

            
cccccc introducimos las variables caracteristicas del sistema      
      
c     npart=100
c      box=11.0d0
      number=idint(dble(npart)**(1.0d0/2.0d0)+1)
      size=longxy/dble(number)
c      temp=1
      
      
cccccc establecemos las posiciones y velocidades de las particulas  
          
       do 300 i=1,number
        do 200 j=1,number
           
           nplace=nplace+1
           if(nplace.le.npart) then
            rx(nplace)=(i-0.5)*size-longxy/2.0d0
            ry(nplace)=(j-0.5)*size-longxy/2.0d0
            rz(nplace)=0.0d0

            vx(nplace)=Ran_Uniform()-0.5d0
            vy(nplace)=Ran_Uniform()-0.5d0
            vz(nplace)=Ran_Uniform()-0.5d0
            sumvx=sumvx+vx(nplace)
            sumvy=sumvy+vy(nplace)
            sumvz=sumvz+vz(nplace)

           endif
 
  200    continue
  300  continue
  
cccccc corregimos las velocidades para que el momento total del sistema sealongx cero

       do 100 i=1,npart
         vx(i)=vx(i)-sumvx/dble(npart)
         vy(i)=vy(i)-sumvy/dble(npart)
         vz(i)=vz(i)-sumvz/dble(npart)
         sumvsq=sumvsq+vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i)
  100  continue
  
cccccc escalamos las velocidades para que la temperatura sea la deseada

       do 400 i=1,npart
         vx(i)=(3.0d0*dble(npart)*temp/sumvsq)**(1.0d0/2.0d0)*vx(i) 
         vy(i)=(3.0d0*dble(npart)*temp/sumvsq)**(1.0d0/2.0d0)*vy(i)
         vz(i)=(3.0d0*dble(npart)*temp/sumvsq)**(1.0d0/2.0d0)*vz(i)
     
         
 400   continue

       return      
       end
   
       


       subroutine inicializacion2(npart,longxy,temp,rx,ry,rz,vx,vy,vz)

ccccc  En esta subrutina vamos a asignar posiciones y velocidades a todas
ccccc las particulas del sistema de forma que consigamos la temperatura T
ccccc que queramos y que la velocidad total sea nula


      implicit none
      integer i,j,k,npart,number,nplace,Ijkl,ii
      double precision rx(npart),ry(npart),rz(npart),
     & size,vx(npart),vy(npart),vz(npart),
     & sumvx,sumvy,sumvz,sumvsq,temp,Ran_Uniform,alea

      double precision longxy
cccccc inicialicemos los acumuladoress 
    
      nplace=0
      sumvx=0
      sumvy=0
      sumvz=0
      sumvsq=0

            
cccccc introducimos las variables caracteristicas del sistema      
      
c     npart=100
c      box=11.0d0
      number=idint(dble(npart)**(1.0d0/2.0d0)+1)
      size=longxy/dble(number)
c      temp=1
      
      
cccccc establecemos las posiciones y velocidades de las particulas  
          
       do 300 i=1,number
        do 200 j=1,number
           
           nplace=nplace+1
           if(nplace.le.npart) then
            rx(nplace)=(i-0.5)*size-longxy/2.0d0
            ry(nplace)=(j-0.5)*size-longxy/2.0d0
            rz(nplace)=0.0d0

           endif
 
  200    continue
  300  continue
 
       do ii=1,npart

          alea=Ran_Uniform()
          if(alea.lt.0.25) then
             vx(ii)=1.0d0+Ran_Uniform()/10.0d0
             vy(ii)=1.0d0+Ran_Uniform()/10.0d0
          else if(alea.ge.0.25.and.alea.lt.0.50) then
             vx(ii)=-1.0d0+Ran_Uniform()/10.0d0
             vy(ii)=-1.0d0+Ran_Uniform()/10.0d0
          else if(alea.ge.0.50.and.alea.lt.0.75) then
             vx(ii)=3.0d0+Ran_Uniform()/10.0d0
             vy(ii)=3.0d0+Ran_Uniform()/10.0d0
          else
             vx(ii)=-3.0d0+Ran_Uniform()/10.0d0
             vy(ii)=-3.0d0+Ran_Uniform()/10.0d0
          endif
          vz(ii)=Ran_Uniform()-0.5d0
          sumvx=sumvx+vx(ii)
          sumvy=sumvy+vy(ii)
          sumvz=sumvz+vz(ii)

          
       enddo
cccccc corregimos las velocidades para que el momento total del sistema sea cero

       do 100 i=1,npart
         vx(i)=vx(i)-sumvx/dble(npart)
         vy(i)=vy(i)-sumvy/dble(npart)
         vz(i)=vz(i)-sumvz/dble(npart)
         sumvsq=sumvsq+vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i)
  100  continue
  
cccccc escalamos las velocidades para que la temperatura sea la deseada

       do 400 i=1,npart
         vx(i)=(3.0d0*dble(npart)*temp/sumvsq)**(1.0d0/2.0d0)*vx(i) 
         vy(i)=(3.0d0*dble(npart)*temp/sumvsq)**(1.0d0/2.0d0)*vy(i)
         vz(i)=(3.0d0*dble(npart)*temp/sumvsq)**(1.0d0/2.0d0)*vz(i)
     
         
 400   continue

       return      
       end
   
       
       subroutine inicializacion3(npart,longxy,temp,rx,ry,rz,vx,vy,vz)

ccccc  En esta subrutina vamos a asignar posiciones y velocidades a todas
ccccc las particulas del sistema de forma que consigamos la temperatura T
ccccc que queramos y que la velocidad total sea nula


      implicit none
      integer i,j,k,npart,number,nplace,Ijkl,ii
      double precision rx(npart),ry(npart),rz(npart),
     & size,vx(npart),vy(npart),vz(npart),
     & sumvx,sumvy,sumvz,sumvsq,temp,Ran_Uniform,alea

      double precision longxy
cccccc inicialicemos los acumuladoress 
    
      nplace=0
      sumvx=0
      sumvy=0
      sumvz=0
      sumvsq=0

            
cccccc introducimos las variables caracteristicas del sistema      
      
c     npart=100
c      box=11.0d0
      number=idint(dble(npart)**(1.0d0/2.0d0)+1)
      size=longxy/dble(number)
c      temp=1
      
      
cccccc establecemos las posiciones y velocidades de las particulas  
          
       do 300 i=1,number
        do 200 j=1,number
           
           nplace=nplace+1
           if(nplace.le.npart) then
            rx(nplace)=(i-0.5)*size-longxy/2.0d0
            ry(nplace)=(j-0.5)*size-longxy/2.0d0
            rz(nplace)=0.0d0

           endif
 
  200    continue
  300  continue
 
       do ii=1,npart

          alea=Ran_Uniform()
          if(alea.lt.0.5) then
          vx(ii)=1.0d0+Ran_Uniform()/10.0d0
          vy(ii)=1.0d0+Ran_Uniform()/10.0d0
          else
          vx(ii)=-1.0d0+Ran_Uniform()/10.0d0
          vy(ii)=-1.0d0+Ran_Uniform()/10.0d0
          endif
          vz(ii)=Ran_Uniform()-0.5d0
          sumvx=sumvx+vx(ii)
          sumvy=sumvy+vy(ii)
          sumvz=sumvz+vz(ii)

          
       enddo
cccccc corregimos las velocidades para que el momento total del sistema sea cero

       do 100 i=1,npart
         vx(i)=vx(i)-sumvx/dble(npart)
         vy(i)=vy(i)-sumvy/dble(npart)
         vz(i)=vz(i)-sumvz/dble(npart)
         sumvsq=sumvsq+vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i)
  100  continue
  
cccccc escalamos las velocidades para que la temperatura sea la deseada

       do 400 i=1,npart
         vx(i)=(3.0d0*dble(npart)*temp/sumvsq)**(1.0d0/2.0d0)*vx(i) 
         vy(i)=(3.0d0*dble(npart)*temp/sumvsq)**(1.0d0/2.0d0)*vy(i)
         vz(i)=(3.0d0*dble(npart)*temp/sumvsq)**(1.0d0/2.0d0)*vz(i)
     
         
 400   continue

       return      
       end
   
       subroutine inicializacion5(npart,longxy,temp,rx,ry,rz,vx,vy,vz)

ccccc  En esta subrutina vamos a asignar posiciones y velocidades a todas
ccccc las particulas del sistema de forma que consigamos la temperatura T
ccccc que queramos y que la velocidad total sea nula


      implicit none
      integer i,j,k,npart,number,nplace,Ijkl,ii
      double precision rx(npart),ry(npart),rz(npart),
     & size,vx(npart),vy(npart),vz(npart),
     & sumvx,sumvy,sumvz,sumvsq,temp,Ran_Uniform,alea

      double precision longxy

cccccc inicialicemos los acumuladoress 
    
      nplace=0
      sumvx=0
      sumvy=0
      sumvz=0
      sumvsq=0

            
cccccc introducimos las variables caracteristicas del sistema      
      
c     npart=100
c      box=11.0d0
      number=idint(dble(npart)**(1.0d0/2.0d0)+1)
      size=longxy/dble(number)
c      temp=1
      
      
cccccc establecemos las posiciones y velocidades de las particulas  
          
       do 300 i=1,number
        do 200 j=1,number
           
           nplace=nplace+1
           if(nplace.le.npart) then
            rx(nplace)=(i-0.5)*size-longxy/2.0d0
            ry(nplace)=(j-0.5)*size-longxy/2.0d0
            rz(nplace)=0.0d0

           endif
 
  200    continue
  300  continue
 
       do ii=1,npart

          alea=Ran_Uniform()
          if(alea.lt.1.0d0/3.0d0) then
          vx(ii)=-3.0d0+Ran_Uniform()/10.0d0
          vy(ii)=-3.0d0+Ran_Uniform()/10.0d0
          else if(alea.ge.1.0d0/3.0d0.and.alea.lt.2.0d0/3.0d0) then
          vx(ii)=1.0d0+Ran_Uniform()/10.0d0
          vy(ii)=1.0d0+Ran_Uniform()/10.0d0
          else 
          vx(ii)=5.0d0+Ran_Uniform()/10.0d0
          vy(ii)=5.0d0+Ran_Uniform()/10.0d0
          endif
          vz(ii)=Ran_Uniform()-0.5d0
                
          sumvx=sumvx+vx(ii)
          sumvy=sumvy+vy(ii)
          sumvz=sumvz+vz(ii)
                    
       enddo	
         
cccccc corregimos las velocidades para que el momento total del sistema sea cero

       do 100 i=1,npart
         vx(i)=vx(i)-sumvx/dble(npart)
         vy(i)=vy(i)-sumvy/dble(npart)
         vz(i)=vz(i)-sumvz/dble(npart)
         sumvsq=sumvsq+vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i)
  100  continue
  
cccccc escalamos las velocidades para que la temperatura sea la deseada

       do 400 i=1,npart
         vx(i)=(3.0d0*dble(npart)*temp/sumvsq)**(1.0d0/2.0d0)*vx(i) 
         vy(i)=(3.0d0*dble(npart)*temp/sumvsq)**(1.0d0/2.0d0)*vy(i)
         vz(i)=(3.0d0*dble(npart)*temp/sumvsq)**(1.0d0/2.0d0)*vz(i)
     
         
 400   continue

       return      
       end
   
       
       
        Subroutine Ran_Init(Ijkl)
        Implicit None
        Integer Ijkl

Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C       This Is The Initialization Routine For Ran_Uniform()                  C
C       Note: The Seed Variable Should Be In The Range 0 <= 900 000 000       C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        Integer I,J,K,L,Ij,Kl
        Integer Ii,Jj,M
        Double Precision S,T

        Double Precision    U(97),C,Cd,Cm
        Integer I97,J97
        Logical Initialised
        Common  /Raset1/ U,C,Cd,Cm,I97,J97,Initialised

        Initialised = .False.

        If( Ijkl .Lt. 0  .Or.  Ijkl .Gt. 900 000 000) Then
            Write(6,*) 'The Random Number Seed Must Have A Value ',
     &                 'Between 0 And 900 000 000'
              Stop
        Endif

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
        End


        Function Ran_Uniform()
        Implicit None

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Generate A Uniformly Distributed Randomnumber Between 0 And 1    C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        Double Precision U(97),C,Cd,Cm,Uni,Ran_Uniform
        Integer I97,J97
        Logical Initialised

        Common /Raset1/ U,C,Cd,Cm,I97,J97,Initialised

        If (.Not. Initialised) Then
          Write(6,*) 'Ran_Uniform: Initialise Ran_Uniform With Ran_Init'
          Stop
        Endif

        Uni = U(I97) - U(J97)
        If( Uni .Lt. 0.0d0 ) Uni = Uni + 1.0d0
        U(I97) = Uni
        I97 = I97 - 1
        If(I97 .Eq. 0) I97 = 97
        J97 = J97 - 1
        If(J97 .Eq. 0) J97 = 97
        C = C - Cd
        If( C .Lt. 0.0d0 ) C = C + Cm
        Uni = Uni - C
        If( Uni .Lt. 0.0d0 ) Uni = Uni + 1.0d0
        Ran_Uniform = Uni

        Return
        End

 


      subroutine tcolinicial(npart,longx,longz,timbig,
     &rx,rz,vx,vz,
     &partner,tcol)
      
ccccc esta subrutina calcula todos los tiempos de colision de todas
ccccc las particulas(lo almacena en la matriz tcol()), asi como con que
ccccc particula colisiona cada una(lo almacena en la matriz partner())
      implicit none
      integer k,i,j,npart,partner(npart+2)
      double precision rx(npart),rz(npart),
     &vx(npart),vz(npart),
     &rxij,rzij,rijsq, 
     &vxij,vzij,vijsq,bij,disc,tij,tcol(npart+3),timbig

      double precision longx,longz,rparedx
ccccc inicializamos la matriz tcol en un tiempo muy grande

      do 100 k=1,npart+2
         tcol(k)=timbig
 100  continue

ccccc caculamos los vectores posicion  y velocidad relativa de los pares i,j

      do 200 i=1,npart-1
         do 300 j=i+1,npart
            rxij=rx(i)-rx(j) 
            
            rzij=rz(i)-rz(j)
            rxij=rxij-longx*anint(rxij/longx)
            

            rijsq=rxij**2.+rzij**2.
            vxij=vx(i)-vx(j)
            
            vzij=vz(i)-vz(j)
            vijsq=vxij**2.+vzij**2.
            bij=rxij*vxij+rzij*vzij

ccccc si es posible la colision entre las particulas (i,j) (bij<0 y disc>0) se ccccc calcula el tiempo de colision y se actualiza la matriz tcol()

            if (bij.lt.0) then

ccccc calculamos el discriminabte

                disc=bij**2-vijsq*(rijsq-1) 
                if (disc.gt.0) then
                   tij=(-bij-disc**(1.0d0/2.0d0))/vijsq
                      if(tij.lt.tcol(i)) then
                        tcol(i)=tij
                        partner(i)=j
                      endif
                      if (tij.lt.tcol(j))then
                         tcol(j)=tij
                         partner(j)=i
                      endif
                endif
             endif

 300     continue
c       print*,i,partner(i),tcol(i)

 200  continue      

cccc  Modificacion USF
cccc  Aqui vamos a ver los tiempos de colision de las particulas con las paredes horizontales
      do i=1,npart

         if (vz(i).gt.0) then

ccc  Calculo el tiempo de colision con la pared de arriba 
ccc  (suponiendo que esta en Lz=longz/2)

            tij=((longz-1.0d0)/2.0d0-rz(i))/vz(i)
           
c     print*,tij
            if (tij.lt.tcol(i))then
               tcol(i)=tij
               partner(i)=npart+1
            endif
            if (tcol(npart+1).gt.tij) then
               tcol(npart+1)=tij
               partner(npart+1)=i
               
            endif
         endif
         
         if (vz(i).lt.0) then

ccc  Calculo el tiempo de colision con la pared de abajo 
ccc  (suponiendo que esta en Lz=-longz/2)
            
            tij=(-(longz-1.0d0)/2.0d0-rz(i))/vz(i)
            
            if (tij.lt.tcol(i))then
               tcol(i)=tij
               partner(i)=npart+2
            endif
            if (tcol(npart+2).gt.tij) then
               tcol(npart+2)=tij
               partner(npart+2)=i

            endif
         endif
c                print*,'con pared',i,partner(i),tcol(i)
      enddo

cccccc   Calculo del tiempo de colision con la pared en elplano x=0 (para calcular presion)
c      do i=1,npart



         
c         if (rx(i)/vx(i).lt.-1.0d-15) then
c          if (rx(i).lt.rparedx.and.vx(i).gt.0) then
  
c            tij=(-rx(i)+rparedx)/vx(i)
c            print*,tij
            
c            if (tij.lt.tcol(i))then
c               tcol(i)=tij
c               partner(i)=npart+3
c            endif
c            if (tcol(npart+3).gt.tij) then
c               tcol(npart+3)=tij
c               partner(npart+3)=i
               
c            endif
c         endif


c         if (rx(i).gt.rparedx.and.vx(i).lt.0) then
  
c            tij=(-rx(i)+rparedx)/vx(i)
c            print*,tij,rx(i)+tij*vx(i),rx(i)+vx(i)/vx(i)*(-rx(i)+1)
            
c            if (tij.lt.tcol(i))then
c               tcol(i)=tij
c               partner(i)=npart+3
c            endif
c            if (tcol(npart+3).gt.tij) then
c               tcol(npart+3)=tij
c               partner(npart+3)=i
               
c            endif
c         endif   

c      enddo

cccc  Fin USF
c      do i=1,npart
c     write(77,*)i,partner(i),tcol(i)
c         print*,i,partner(i),tcol(i)
c      enddo
      return
      end

      subroutine tcolminimo(npart,tcol,partner,timbig,tcolmin,ii,jj)

ccccc Esta subrutina calcula el minimo de la matriz tcol() que llamaremos 
ccccc tcolmin, asi como la pareja que va a colisionar ii,jj
       
      implicit none
      integer i,ii,jj,npart,partner(npart+2)
      double precision tcol(npart+2),tcolmin,timbig
      tcolmin=timbig
      do 100 i=1,npart+2
         if (tcol(i).lt.tcolmin) then
            tcolmin=tcol(i)
            ii=i
            jj=partner(i)
         
         endif
 100  continue
cc      if (ii.gt.jj.and.ii.eq.npart+1) then
cccc  print*,ii,jj,tcolminc
c     c              call flush()
cc            endif
      return
      end























