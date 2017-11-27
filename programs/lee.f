c     gfortran -o lee lee.f
c      character*82 filename
      character*91 filename
      character*200 filein,fileout

      parameter (nmax=256**3)
      dimension x(nmax),y(nmax),z(nmax)      
      dimension vx(nmax),vy(nmax),vz(nmax)
      dimension pmass(nmax)
      dimension u(nmax),rho(nmax),el_ab(nmax),h_frac(nmax)
      dimension smooth(nmax),sfr(nmax),rmfs(nmax),fe_ab(nmax),o_ab(nmax)
      dimension epot(nmax)
c      dimension age(nmax)
      integer   id(nmax)
      real*8    rid(nmax) 
 
      integer*4 np(6)                  !Number of particles in each species
      integer*4 nall(6)                !Total number in each species
      real*8 massarr(6)                !Particle masses
      real*8 aexp,redshift             !Expansion factor, redshift
      integer*4 flagsfr,flagfeedback   !Star formation rate, feedback
      integer*4 flagcooling            !Cooling
      integer*4 NumFiles               !Number of files in snapshot
      real*8 BoxSize                   !Size of cosmological box
      real*8 omega0,omegalambda        !Matter, Vacuum energy density at z=0
      real*8 HubbleParam               !Value of Hubble's constant, h, in units of 100 kms^{-1} Mpc^{-1}
      real*4 a,omegaz,omegal
      character unused(256-6*4-6*4-6*8-2*8-4*4-4*8)   !Buffer

      character*3 :: istep_aux

      istep=497

      write(istep_aux, '(i3)') istep 
      filein = trim('/store/01/clues/B64_WM3_186592/LG/GAS_SFR/4096/SNAPS/snap_')//istep_aux

      !write(filein,'(''/store/01/clues/B64_WM3_186592/LG/GAS_SFR/4096/SNAPS/snap_'',i3.3)')istep
      write(*,*) 'Opening file=',filein
      open(1,file=filein,status='old',form='unformatted')
 
      read(1) np,massarr,aexp,redshift,flagsfr,flagfeedback
     &     ,nall,flagcooling,NumFiles,BoxSize,omega0,omegalambda
     &     ,HubbleParam
     &     ,unused
      write(*,*) 'Read header'
      stop

      nwmz=0
      do i=1,6
         if(massarr(i).eq.0)nwmz=nwmz+np(i)
         write(*,*) 'Mass in species',i,np(i),nall(i),massarr(i)
      end do
        write(*,*)'nwmz=',nwmz

        ngas=np(1)
        ndrk=np(2)
        nbn1=np(3)
        nbn2=np(4)
        nstr=np(5)
        nbn3=np(6)
        ntot=ngas+ndrk+nbn1+nbn2+nstr+nbn3
        nbar=ngas+nstr
        write(*,*)ntot,Nmax
 
      if(ntot.gt.Nmax) stop 'Need to increase Nmax'
 
      write(*,*) 'Number of particles: ',ntot
      write(*,*) 'aexp,z=',aexp,redshift
      write(*,*) 'BoxSize,HubbleParam,=',BoxSize,HubbleParam
      write(*,*) 'omega0,omegalambda,=',omega0,omegalambda
      write(*,*)
      read(1) (x(i),y(i),z(i),i=1,ntot)
      print *, 'Read positions'
      read(1) (vx(i),vy(i),vz(i),i=1,ntot)
      print *, 'Read velocities'
      read(1) (id(i), i=1,ntot)
      print *, 'Read identities'
      read(1) (pmass(i),i=1,nwmz)
      print *, 'Read masses'
      read(1) (u(i),i=1,ngas)
      print *, 'Read internal energy'
      read(1) (rho(i),i=1,ngas)
      print *, 'Read density'
      read(1) (el_ab(i),i=1,ngas)
      print *, 'Read free electron abundance'
      read(1) (h_frac(i),i=1,ngas)
      print *, 'Read Hydogen neutral fraction'
      read(1) (smooth(i),i=1,ngas)
      print *, 'Read smoothing lenght'
      read(1) (sfr(i),i=1,ngas)
      print *, 'Read sfr'
      read(1) (rmfs(i),i=1,nstr)
      print *, 'Read rmfs'
        if(jfile.ge.8.and.jfile.le.21)then
      read(1) (fe_ab(i),i=1,nbar)
      print *, 'Read fe_ab'
      read(1) (o_ab(i),i=1,nbar)
      print *, 'Read o_ab'
        endif
      read(1) (epot(i),i=1,ntot)
      print *, 'Read potential energy'
      close(1)

        if(jfile.ge.1.and.jfile.le.7)then
        do i=1,nbar
        fe_ab(i)=999999.
        o_ab(i) =888888.
        enddo
        endif
c-----------------------------------------------------------

      rkm2Mpc=3.24077885e-20
      Gyr2sec=3.1536e16
      h0=73.
      h_little=h0/100.
      h=h0*sqrt((1.+redshift)**3*omega0+omegalambda)

      do i=1,6
      massarr(i)=massarr(i)/h_little
      enddo

      do i=1,ntot
      rid(i)=1.0D0*id(i)
      pmass(i)=pmass(i)/h_little ! check with Franziska all these equations
      x(i)=x(i)*aexp/h_little
      y(i)=y(i)*aexp/h_little
      z(i)=z(i)*aexp/h_little
      vx(i)=vx(i)*sqrt(aexp)+h*x(i)/1000.
      vy(i)=vy(i)*sqrt(aexp)+h*y(i)/1000.
      vz(i)=vz(i)*sqrt(aexp)+h*z(i)/1000.
      enddo

c        a=aexp !age define with respect to the actual snapshot being read
c        a=1.0  !age define with respect to z=0, not to actual snapshot being read 
c      omegaz=omega0
c      omegal=omegalambda
c      timeGyr=alog(a**(3./2.)/sqrt(omegaz/omegal)+sqrt(1.+(a**(3./2.)/sqrt(omegaz/omegal))**2))*2./3./sqrt(omegal)/h0/rkm2Mpc/Gyr2sec
c      do i=1,nstr
c      age(i)=timeGyr-alog(rmfs(i)**(3./2.)/sqrt(omegaz/omegal)+sqrt(1.+(rmfs(i)**(3./2.)/sqrt(omegaz/omegal))**2))*2./3./sqrt(omegal)/h0/rkm2Mpc/Gyr2sec
c      enddo
c      write(*,*)timeGyr,omegaz,omegal,h0,rkm2Mpc,Gyr2sec
c      print *, 'Convert to physical units'
c----------------- WRITE ASCII -------------
      write(fileout,'(''../data/'',i2.2,''/snapshot_gas_'',i3.3,''.dat'')')ifile,istep
      open(2,file=fileout,status='unknown')
      do j=1,ngas
      i=j
      write(2,'(17e18.10)')pmass(i),x(i),y(i),z(i),vx(i),vy(i),vz(i),rid(i),u(i),rho(i),el_ab(i),h_frac(i),smooth(i),sfr(i),fe_ab(i),o_ab(i),epot(i)
      end do
      close(2)
      print *, 'Write gas'

      write(fileout,'(''../data/'',i2.2,''/snapshot_drk_'',i3.3,''.dat'')')ifile,istep
      open(2,file=fileout,status='unknown')
      do j=1,ndrk
      i=j+ngas
      write(2,'(9e18.10)')massarr(2),x(i),y(i),z(i),vx(i),vy(i),vz(i),rid(i),epot(i)
      end do
      close(2)
      print *, 'Write drk'

      write(fileout,'(''../data/'',i2.2,''/snapshot_bn1_'',i3.3,''.dat'')')ifile,istep
      open(2,file=fileout,status='unknown')
      do j=1,nbn1
      i=j+ngas+ndrk
      write(2,'(9e18.10)')massarr(3),x(i),y(i),z(i),vx(i),vy(i),vz(i),rid(i),epot(i)
      end do
      close(2)
      print *, 'Write bn1'

      write(fileout,'(''../data/'',i2.2,''/snapshot_bn2_'',i3.3,''.dat'')')ifile,istep
      open(2,file=fileout,status='unknown')
      do j=1,nbn2
      i=j+ngas+ndrk+nbn1
      write(2,'(9e18.10)')massarr(4),x(i),y(i),z(i),vx(i),vy(i),vz(i),rid(i),epot(i)
      end do
      close(2)
      print *, 'Write bn2'

      write(fileout,'(''../data/'',i2.2,''/snapshot_str_'',i3.3,''.dat'')')ifile,istep
      open(2,file=fileout,status='unknown')
      do j=1,nstr
      i=j+ngas+ndrk+nbn1+nbn2
      write(2,'(12e18.10)')pmass(ngas+j),x(i),y(i),z(i),vx(i),vy(i),vz(i),rid(i),rmfs(j),fe_ab(ngas+j),o_ab(ngas+j),epot(i)
      end do
      close(2)
      print *, 'Write str'

      write(fileout,'(''../data/'',i2.2,''/snapshot_bn3_'',i3.3,''.dat'')')ifile,istep
      open(2,file=fileout,status='unknown')
      do j=1,nbn3
      i=j+ngas+ndrk+nbn1+nbn2+nstr
      write(2,'(9e18.10)')pmass(nbar+j),x(i),y(i),z(i),vx(i),vy(i),vz(i),rid(i),epot(i)
      end do
      print *, 'Write bn3'
      close(2)

      write(fileout,'(''../data/'',i2.2,''/npar_'',i3.3,''.dat'')')ifile,istep
      open(2,file=fileout)
      write(2,*)ngas,ndrk,nstr,nbn1,nbn2,nbn3
      close(2)

      enddo
      end
