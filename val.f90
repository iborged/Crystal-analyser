program xtal analyser
  implicit none
  integer :: nc, nmo, suma, ang
  integer :: i,j,m,num(1000),nas,n,suma2(1000), ang_raw
  real(kind=8) :: dx=0.0, dy=0.0, dz=0.0, dist=0.0, sx=0.0, sy=0.0, sz=0.0
  real(kind=8) :: sum_bv_moc(1000), sum_bv_momo(1000), sum_bv_cc(1000), sum_bv_cmo(1000)
  real(kind=8) :: bv=0.0, sum_bv_dev=0.0, op=0.0, dev=0.0, cal_eng=0.0
  real(kind=8) :: pi, c(3), vectores(1000,4),crd(48,3),a(3,3)
  !
  real(kind=8) :: alpha=0.0, norma1=0.0, norma2=0.0, normat=0.0, xx=0.0, yy=0.0, zz=0.0, pp=0.0
  !
  pi=4*atan(1.0)
  !
  nc=16
  nmo=32
  !
  !suma: es la sumatoria de los carbonos al cuadrado
  !
  !
  !  open(30,file="struct")          !este file en teoria tiene una estructura completa con un arreglo en particular
  !33 format(1x,f10.5,1x,f10.5,1x,f10.5)
  !
  do i=1,1000              !inicializando estos arreglos, es decir, metiendo ceros en cada espacio de la memoria.
     sum_bv_moc(i)=0.0
     sum_bv_momo(i)=0.0
     sum_bv_cc(i)=0.0
     sum_bv_cmo(i)=0.0
     num(i)=0
     suma2(i)=0            !esto es para el parametro de orden
     vectores(i,:)=0.0
  enddo
  !
  !
  open(unit=86,file='data')
  !open(unit=86,file='1.vasp')
  !
  read(86,*)
  read(86,*)
  do i=1,3
     read(86,*) a(i,1), a(i,2), a(i,3)
  enddo
  read(86,*)
  read(86,*)
  read(86,*)
  do i=1,48
     read(86,*) crd(i,:)
  enddo
  !
  !write(*,*) a(1,1),a(2,2),a(3,3)
  !stop
  !
  c(1)=1.947E-03    !c²  from Mo2C - tercer intento, con un r²=0.9404 
  c(2)=4.972E-01    !bv  from Mo2C - tercer intento, con un r²=0.9404 
  c(3)=3.501E-03    !ang from Mo2C - tercer intento, con un r²=0.9404 
  !
  !c(1)=0.0014*13.605668       !c²      (ev/c²)  |  mo2c - ib usando R, con un r²=0.957 que se obtuvieron de un sistema 2x2x2
  !c(2)=0.6912*13.605668       !bv      (ev/dev) |  mo2c - ib usando R, con un r²=0.957 que se obtuvieron de un sistema 2x2x2
  !c(3)=0.0069*13.605668       !angles  (ev/ang) |  mo2c - ib usando R, con un r²=0.957 que se obtuvieron de un sistema 2x2x2
  !
  sx=a(1,1)     !lattice x
  sy=a(2,2)     !lattice y
  sz=a(3,3)     !lattice z
  !
  !
  !                               ###########################################################################
  !                               ###########################  carbonos al  #################################
  !                               ###########################   cuadrado    #################################
  !                               ###########################################################################
  !
  do i=1,nmo
     n=0
     do j=nmo+1, nmo+nc
        dx=crd(i,1)-crd(j,1)
        if (dx.gt.sx/2.0) dx=dx-sx
        if (dx.lt.-sx/2.0) dx=dx+sx
        !
        dy=crd(i,2)-crd(j,2)
        if ( dy.gt.(sy/2.0) ) dy=dy-sy
        if ( dy.lt.(-sy/2.0)) dy=dy+sy
        !
        dz=crd(i,3)-crd(j,3)
        if ( dz.gt.(sz/2.0) ) dz=dz-sz
        if ( dz.lt.(-sz/2.0)) dz=dz+sz
        !
        dist=sqrt(dx**2+dy**2+dz**2)
        !
        if (dist.lt.2.4) then
           n=n+1
        endif
     enddo
     num(i)=n                 ! estoy guardando/almacenando cada valor de n en el arreglo num(i), eso es lo que se hace aca.
  enddo
  !
  suma=0
  do i=1,nmo
     suma=suma+num(i)**2     !aca esta haciendo la sumatoria de los carbonos al cuadrado
     suma2(i)=num(i)**2
     !  write(*,*)i,num(i)**2
  enddo
  !
  op=0.0
  op=sqrt(sum(suma2,dim=1)*(1.0/nmo))  !parametro de orden, root mean square (rms)
  !
  ! write(*,*)suma   !valor final de la sumatoria de los carbonos al cuadrado
  !
  !
  !                               ###########################################################################
  !                               ###########################              ##################################
  !                               ########################### bond valence ##################################
  !                               ###########################              ##################################
  !                               ###########################################################################
  !
  !
  !bond valence for mo-c -- i
  !r is the length of a bond between the two given atoms. in this code r is represented by "dist"
  !the bond valence has the property that its sum around each atom in a compound is equal to the valence (oxidation state) of that atom.
  !rmoc=1.877
  !
  !open(unit=1000, file="bv_moc.dat", action="write")
  bv=0.0
  do i=1,nmo                        !este bucle have bv para mo-c, usando de 1-32 estoy usando los mo
     !n=0
     bv=0.0
     do j=nmo+1,nmo+nc                    !aca estoy usando los c
        if (i.ne.j) then
           dx=crd(i,1)-crd(j,1)
           if (dx.gt.sx/2.0) dx=dx-sx
           if (dx.lt.-sx/2.0) dx=dx+sx
           !
           dy=crd(i,2)-crd(j,2)
           if (dy.gt.sy/2.0) dy=dy-sy
           if (dy.lt.-sy/2.0 ) dy=dy+sy
           !
           dz=crd(i,3)-crd(j,3)
           if (dz.gt.sz/2.0) dz=dz-sz
           if (dz.lt.-sz/2.0) dz=dz+sz
           !
           dist=sqrt(dx**2+dy**2+dz**2)
           !
           if (dist.lt.2.3) then
              bv=exp(((1.985-dist)*(1.0/0.37)))    !mo-c bond-valence parameter
              sum_bv_moc(i)=sum_bv_moc(i)+bv
              !n=n+1
              !write(111,*)i,(j-32)
           endif
        endif
     enddo
     !write(*,*)sum_bv_moc(i),i
     !write(111,*)i,j
  enddo
  !write(*,*)
  !
  !
  !bond valence for mo-mo -- ii
  !
  !open(unit=1001, file="bv_momo.dat", action="write")
  !write(1001,*)"    bv momo         mo atom"
  !
  bv=0.0
  do i=1,nmo    !este buble hace bv para mo-mo
     !n=0
     bv=0.0
     do j=1,nmo
        if (i.ne.j) then
           dx=crd(i,1)-crd(j,1)
           if (dx.gt.sx/2.0) dx=dx-sx
           if (dx.lt.-sx/2.0) dx=dx+sx
           !
           dy=crd(i,2)-crd(j,2)
           if (dy.gt.sy/2.0) dy=dy-sy
           if (dy.lt.-sy/2.0) dy=dy+sy
           !
           dz=crd(i,3)-crd(j,3)
           if (dz.gt.sz/2.0) dz=dz-sz
           if (dz.lt.-sz/2.0) dz=dz+sz
           !
           dist=sqrt(dx**2+dy**2+dz**2)
           !
           if (dist.lt.3.3) then
              bv=exp( ((2.615-dist)*(1.0/0.37)) )    !mo-mo bond-valence parameter
              sum_bv_momo(i)=sum_bv_momo(i)+bv
              !n=n+1
           endif
        endif
     enddo
     !write(*,*)sum_bv_momo(i),i
  enddo
  !write(*,*)
  !
  !
  !bond valence for c-c -- iii
  !rcc=1.540
  !
  !open(unit=1001, file="bv_cc.dat", action="write")
  !write(1001,*)"    bv cc         c atom"
  !
  bv=0.0
  do i=nmo+1,nmo+nc
     !n=0
     bv=0.0
     do j=nmo+1, nmo+nc                      !c atoms
        if (i.ne.j) then
           dx=crd(i,1)-crd(j,1)
           if (dx.gt.sx/2.0) dx=dx-sx
           if (dx.lt.-sx/2.0) dx=dx+sx
           !
           dy=crd(i,2)-crd(j,2)
           if (dy.gt.sy/2.0) dy=dy-sy
           if (dy.lt.-sy/2.0) dy=dy+sy
           !
           dz=crd(i,3)-crd(j,3)
           if (dz.gt.sz/2.0) dz=dz-sz
           if (dz.lt.-sz/2.0) dz=dz+sz
           !
           dist=sqrt(dx**2+dy**2+dz**2)
           !
           if (dist.lt.3.4) then
              bv=exp( ((1.540-dist)*(1.0/0.37)) )    !c-c bond-valence parameter
              sum_bv_cc(i)=sum_bv_cc(i)+bv
           endif
        endif
     enddo
     !write(*,*) sum_bv_cc(i),i
  enddo
  !write(*,*)
  !
  !bv for c-mo -- iv
  !
  !open(unit=1001, file="bv_cmo.dat", action="write")
  !write(1001,*)"    bv cmo         mo atom"
  !
  bv=0.0
  do i=nmo+1,nmo+nc             !atomos de carbono.
     !n=0
     bv=0.0
     do j=1,nmo             !atomos de mo
        if (i.ne.j) then
           dx=crd(i,1)-crd(j,1)
           if (dx.gt.sx/2.0) dx=dx-sx
           if (dx.lt.-sx/2.0) dx=dx+sx
           !
           dy=crd(i,2)-crd(j,2)
           if (dy.gt.sy/2.0) dy=dy-sy
           if (dy.lt.-sy/2.0) dy=dy+sy
           !
           dz=crd(i,3)-crd(j,3)
           if (dz.gt.sz/2.0) dz=dz-sz
           if (dz.lt.-sz/2.0) dz=dz+sz
           !
           dist=sqrt(dx**2+dy**2+dz**2)
           !
           if (dist.lt.2.3) then
              bv=exp( ((1.985-dist)*(1.0/0.37)) )    !c-mo bond-valence parameter
              sum_bv_cmo(i)=sum_bv_cmo(i)+bv
           endif
        endif
     enddo
     !write(*,*)sum_bv_cmo(i),i
  enddo
  !write(*,*)
  !
  sum_bv_dev=0.0
  do i=1,nmo
     sum_bv_dev=sum_bv_dev + ((sum_bv_moc(i) + sum_bv_momo(i)-6)**2)             !sum of all bv values for mo.
  enddo                                                                          !cual es el estado de oxidacion en el mo2c ??
  !
  !write(*,*)sum_bv_dev
  !
  !do i=1,32
  !   write(*,*) sum_bv_moc(i) + sum_bv_momo(i)                                  !aca imprime la valencia de cada mo individual
  !enddo
  !
  !write(*,*)
  dev=0.0
  dev=sqrt((sum_bv_dev/nmo))  !desviacion promedio de los atomos de mo.
  !write(*,*)dev
  !write(*,*)sum_bv_dev,dev
  !
  !write (*,*)suma,dev
  !
  !
  !                               ###########################################################################
  !                               ###########################              ##################################
  !                               ###########################  angles 180  ##################################
  !                               ###########################              ##################################
  !                               ###########################################################################
  !
  nas=0.0                       !el truco para el contador raro, tratar de entender.
  do i=1,nmo                    !almacenando las coordenadas de molybdeno
     n=0
     xx=0.0
     yy=0.0
     zz=0.0
     do j=nmo+1,nmo+nc               !atomos de carbono
        dx=crd(i,1)-crd(j,1)
        if (dx.gt.(sx/2.0)) dx=dx-sx
        if (dx.lt.(-sx/2.0)) dx=dx+sx
        !
        dy=crd(i,2)-crd(j,2)
        if (dy.gt.(sy/2.0)) dy=dy-sy
        if (dy.lt.(-sy/2.0)) dy=dy+sy
        !
        dz=crd(i,3)-crd(j,3)
        if (dz.gt.(sz/2.0)) dz=dz-sz
        if (dz.lt.(-sz/2.0)) dz=dz+sz
        !
        dist=sqrt(dx**2+dy**2+dz**2)
        !
        if (dist.lt.2.4) then
           n=n+1
           xx=dx !crd(i,1)-crd(j,1)         !al hacer esta resta de un punto menos otro punto, esto se convierte en un vector (imp. no olvidar).
           yy=dy !crd(i,2)-crd(j,2)         !por eso es que hago la resta de coordenada por coordenada, para obtener cada vector.
           zz=dz !crd(i,3)-crd(j,3)
           !42 format(f10.1,1x,f10.5,1x,f10.5,1x,f10.5)
           !write(*,42)i,xx,yy,zz
           nas=nas+1
           vectores(nas,:)=(/dble(i), xx, yy, zz/)   !aprender que este es el formato que se usa el / / para encerrar las cantidades deseadas.
           !write(*,*)i, n, vectores(nas,:)
           !write(*,42)vectores(nas,:)
        endif
     enddo
  enddo
  !
  ang=0     !si quito el comentario para 'ang=0' los angulos en el annealing son iguales a cero, asi que no quitar nunca.
  ang_raw=0
  do i=1,1000
     do n=1,1000
        pp=0.0
        alpha=0.0
        norma1=0.0
        norma2=0.0
        normat=0.0
        if (int(vectores(i,1)).eq.int(vectores(n,1))) then
           pp=dot_product((vectores(i,2:4)),(vectores(n,2:4)))
           norma1=sqrt(((vectores(i,2)**2) + (vectores(i,3)**2) + (vectores(i,4)**2)))
           norma2=sqrt(((vectores(n,2)**2) + (vectores(n,3)**2) + (vectores(n,4)**2)))
           normat=(norma1*norma2)
           alpha=(180/pi)*acos((pp/normat))
           !print*,alpha
           if (alpha.gt.171.and.alpha.lt.185) then
              !print*,alpha
              ang_raw=ang_raw+1
           endif
           ang=ang_raw/2
        endif
     enddo
  enddo
  !
  !
  !open(unit=1000, file="results.dat", action="write")
  !
  cal_eng=0.0
  !cal_eng=cal_eng + (c(1)*suma) + (c(2)*dev) + (c(3)*(ang) )
  cal_eng=cal_eng + (c(1)*suma) + (c(2)*dev) + (c(3)*(ang) - 4.644E03 )
  !
  !34 format(i4,1x,f10.5,1x,i6,1x,f10.6)
  write(*,'(i4,1x,f10.5,1x,i6,2x,f15.8)') suma, dev, ang, cal_eng
  !
  !write(*,*)suma, dev, ang/2, op
  !write(*,*)

endprogram xtal analyser
