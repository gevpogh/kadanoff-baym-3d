!  KB3d.F
      program kb3d
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                                       !c
!c     Nov. 1998                                                         !c
!c                                                                       !c
!c     Authors: H. Sigurd Kohler, Nai-Hang Kwong                         !c
!c              University of Arizona                                    !c
!c                                                                       !c
!c              Hashim A. Yousif                                         !c
!c              University of Pittsburgh@bradford                        !c
!c                                                                       !c
!c     This code calculates the two-time non-equilibrium Green           !c
!c     functions for a homogeneous fermion system.  Starting             !c
!c     from a given input momentum distribution, the code solves         !c
!c     the Kadanoff-Baym transport equation, the self energy             !c
!c     being approximated at the second-order direct Born level.         !c
!c     Two illustrative systems are included : the electron gas          !c
!c     (modeling the conduction band electrons or valence band           !c
!c     holes in an excited semiconductor) and nuclear matter             !c
!c     (modeling heavy nuclei in collisions).                            !c
!c                                                                       !c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      include 'kb3dinc.f90'
!c
      pi=acos(-1.)
      eye=cmplx(0.,1.)
!c
      open (unit=10,file='kb3d.dat',form='formatted',status='old')
      open (unit=6,file='kb3d.out',status='new')
!c
      read(10,*) isy                   ! 1 : electrons; 2 : nucleons
!c
!c     electron parameters
!c
      read(10,*) aindb,eryd,amee,epsb,isck,akp
      read(10,*) ra,rhb,rgam,rgam1,rgam2
!c
!c     nucleon parameters
!c
      read(10,*) p0,pf
!c
      read(10,*) dpz,dft
      read(10,*) ntimes,time0         ! ntimes < or = ntl-1
      read(10,*) istf,mod_store,ifk   ! meanf
      read(10,*) nplot                ! nplot < or = 10
      if (nplot .gt. 0) read(10,*) (iplot(i),i=1,nplot) 
!c
      close(unit=10,status='keep')
!c
!c     set up grid labels and initialize common arrays and 
!c     fft coefficients
!c
      call init_grid(mod_store)
!c
      call init_var
!c
!c     the following call to initialize the fft routines
      
!c
      call cfft3di(nvec,nvec,nvec,ftcoef)
!c
!c     set up parameters and initial Green functions
!c
      if (isy .eq. 1) call init_elec
      if (isy .eq. 2) call init_nucl
!c
!c     calculate kinetic energy propagator
!c
      call epprop
!c
      call sigmas(1)

      if (ifk .eq. 1) call sig_fk(1)    ! meanf
!c
!c     test
!c 
!c      do 500 k=0,idim
!c      iq=irq(0,0,k)
!c      iqx=irq(k,0,0)
!c      write (*,*) k,smf(iq)/h2m2,smf(iqx)/h2m2
!c 500  continue
!c     
      time=time0
!c
!c     start time loop
!c
      do 1000 nt=1,ntimes
!c
      call tstepi(nt,1)
!c
      call out_stat(nt)
!c
      call sigmas(nt+1)
      if (ifk .eq. 1) call sig_fk(nt+1)    ! meanf
!c
      if (istf .eq. 1) then
      call tstepi(nt,2)
      call sigmas(nt+1)
      if (ifk .eq. 1) call sig_fk(nt+1)    ! meanf
      endif
!c
      if (isy .eq. 1 .and. isck .eq. 1) then
      call calakp(nt)
      call potcoul
      endif
!c
      time=time+dft
!c
1000  continue
!c
!c     write the parameters and the output of the program
 
      call pwrite(ntimes,time0,istf,mod_store)
!c
      stop
      end
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine init_grid(mod_store)
!c
!c     set up map between storage labels and momentum-space
!c     Cartesian grid points.  more details are given in section 3.3.a of the 
!c     long write-up
!c
!c     mod_store = 0 : full storage mode
!c                 1 : full storage in one octant
!c                 2 : spherical symmetry
!c                 3 : cylindrical symmetry
!c
      include 'kb3dinc.f90'

      dimension nfreq(0:idim**2),labelq(0:idim**2),&
               & nfreqc(-idim:idim,0:idim**2),&
               & labelc(-idim:idim,0:idim**2)
!c
      do 1 i=-idim,idim
      do 1 j=-idim,idim
      do 1 k=-idim,idim
      irq(i,j,k)=0
   1  continue
!c
!c     map (i,j,k) -> iq
!c
      if (mod_store .eq. 0) then
!c
      iqn=0
      do 11 i=-idim,idim
      do 11 j=-idim,idim
      do 11 k=-idim,idim
      ip2=i**2+j**2+k**2
      if(ip2.le.idim**2) then
      iqn=iqn+1
      irq(i,j,k)=iqn
      endif
  11  continue
!c
      else if (mod_store .eq. 1) then
!c
      iqn=0
      do 101 i=0,idim
      do 101 j=0,idim
      do 101 k=0,idim
      ip2=i**2+j**2+k**2
      if(ip2.le.idim**2) then
      iqn=iqn+1
      irq(i,j,k)=iqn
      endif
 101  continue
!c
      else if (mod_store .eq. 2) then
!c
      do 201 ip2=0,idim**2
      nfreq(ip2)=0
 201  continue
!c
      do 210 i=0,idim
      do 210 j=0,idim
      do 210 k=0,idim
      ip2=i**2+j**2+k**2
      if(ip2.le.idim**2) then
      nfreq(ip2)=nfreq(ip2)+1
      endif
 210  continue
!c
      iqn=0
      do 220 ip2=0,idim**2
      if (nfreq(ip2) .gt. 0) then
      iqn=iqn+1
      labelq(ip2)=iqn
      endif
 220  continue
!c
      do 230 i=0,idim
      do 230 j=0,idim
      do 230 k=0,idim
      ip2=i**2+j**2+k**2
      if(ip2.le.idim**2) then
      irq(i,j,k)=labelq(ip2)
      endif
 230  continue      
!c
      else if (mod_store .eq. 3) then
!c
      do 301 ipc2=0,idim**2
      do 301 k=-idim,idim
      nfreqc(k,ipc2)=0
 301  continue
!c
      do 310 i=0,idim
      do 310 j=0,idim
      ipc2=i**2+j**2
      do 311 k=-idim,idim
      ip2=ipc2+k**2
      if(ip2.le.idim**2) then
      nfreqc(k,ipc2)=nfreqc(k,ipc2)+1
      endif
 311  continue
 310  continue
!c
      iqn=0
      do 320 ipc2=0,idim**2
      do 320 k=-idim,idim
      if (nfreqc(k,ipc2) .gt. 0) then
      iqn=iqn+1
      labelc(k,ipc2)=iqn
      endif
 320  continue
!c
      do 330 i=0,idim
      do 330 j=0,idim
      ipc2=i**2+j**2
      do 330 k=-idim,idim
      ip2=ipc2+k**2
      if(ip2.le.idim**2) then
      irq(i,j,k)=labelc(k,ipc2)
      irq(-i,j,k)=irq(i,j,k)
      irq(i,-j,k)=irq(i,j,k)
      irq(-i,-j,k)=irq(i,j,k)
      endif
 330  continue      
!c
      endif      
!c
      if(iqn.ne.nq) then
      write(*,*) 'nq should be', iqn
      write(*,*) 'it is       ', nq
      stop
      endif
!c
      if (mod_store .eq. 1 .or. mod_store .eq. 2) then
        do 400 i=0,idim
        do 400 j=0,idim
        do 400 k=0,idim
        iq=irq(i,j,k)
        irq(i,j,-k)=iq
        irq(i,-j,k)=iq
        irq(i,-j,-k)=iq
        irq(-i,j,k)=iq
        irq(-i,j,-k)=iq
        irq(-i,-j,k)=iq
        irq(-i,-j,-k)=iq
 400    continue
      endif
!c
!c     iq -> (i,j,k)
!c
      do 1001 i=-idim,idim
      do 1001 j=-idim,idim
      do 1001 k=-idim,idim
      ip2=i**2+j**2+k**2
      if(ip2.le.idim**2) then
      iq=irq(i,j,k)
      ix(iq)=i
      jx(iq)=j
      kx(iq)=k
      endif
1001  continue
!c
      return
      end
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine init_var
!c
!c     initialize arrays and variables
!c
      include 'kb3dinc.f90'

      do 1 it=1,ntl
      do 1 jt=1,ntl
      do 1 iq=1,nq
      gre(iq,it,jt)=cmplx(0.,0.)
 1    continue
!c
      do 2 it=1,ntl
      do 2 iq=1,nq
      sgle(iq,it)=cmplx(0.,0.)
      sgla(iq,it)=cmplx(0.,0.)
 2    continue
!c
      do 3 iq=1,nq
      dpot(iq)=cmplx(0.,0.)
      smf(iq)=cmplx(0.,0.)          ! meanf
 3    continue
!c
      do 4 i=-idim,idim
      do 4 j=-idim,idim
      do 4 k=-idim,idim
      pot(k,j,i)=0.
 4    continue
!c
      do 5 i=1,nvec                ! meanf
      do 5 j=1,nvec
      do 5 k=1,nvec
      pot_fk(k,j,i)=cmplx(0.,0.)
 5    continue
!c
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine init_elec
!c
!c     set up parameters and initial Green functions for
!c     the electron gas system (semiconductor conduction band). the initial 
!c     distribution is given by eq.(12) of the long write-up.
!c     
      include 'kb3dinc.f90'

      h2m20=3.80998          ! eV-A^2
      h2m2=h2m20/amee
      anu=2.
      hbar=.658212           ! eV-fs
      dt=dft/hbar            ! eV^-1
      dpz=dpz/aindb          ! A^-1
      akp=akp/aindb
      dp3=(dpz/(2.*pi))**3
      coul=14.3997/epsb      ! eV-A
      vcd3d=coul*dpz*(pi+2.*.34773)/(pi**2)   ! meanf
      rgam22=rgam**2
!c
!c     initial Green functions
!c
      do 10 iq=1,nq
      i=ix(iq)
      j=jx(iq)
      k=kx(iq)
      pz=dpz*float(k)
      py=dpz*float(j)
      px=dpz*float(i)
      pl2=pz**2+py**2+px**2
      ree=h2m2*pl2
      rehh=h2m20*(rgam1-2.*rgam2)*pl2
      relh=h2m20*(rgam1+2.*rgam2)*pl2
      a1=(rhb-ree-rehh)**2
      rgh=exp(-a1/(2.*rgam22))
      a1=(rhb-ree-relh)**2
      rgl=exp(-a1/(2.*rgam22))
!c
      pl=sqrt(pl2)
      cth=pz/(pl+.000001)
      cth2=cth**2
      a1=(1.-cth2)/2.
      a2=(1.+3.*cth2)/6.
!c
      f=ra*(rgh*a1+rgl*a2)
      gre(iq,1,1)=eye*f
!c
 10   continue
!c
      if (isck .eq. 1) call calakp(1)
      call potcoul
      call potcoul_fk                   ! meanf
!c
      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine init_nucl
!c
!c     set up parameters and initial Green functions for
!c     the nucleon system (symmetric nuclear matter).
!c     the initial distribution is given by eq.(13) of the long write-up.
!c
      include 'kb3dinc.f90'

      h2m2=20.7355              ! MeV-fm^2
      anu=4.
      hbarc=197.327             ! MeV-fm
      dt=dft/hbarc              ! MeV^-1
      dp3=(dpz/(2.*pi))**3
      aindb=1.
!c
      do 10 iq=1,nq
      i=ix(iq)
      j=jx(iq)
      k=kx(iq)
      pz=dpz*float(k)
      py=dpz*float(j)
      px=dpz*float(i)
      apz=abs(pz)
      pp=sqrt((apz-p0)**2+py**2+px**2)
      f=0.
      if (pp .lt. pf) f=1.
      gre(iq,1,1)=eye*f
 10   continue
!c
      call potpaw
      call potnuc_fk               ! meanf
!c
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine epprop
!c
!c     calculate single particle energy propagator.
!c     more details are given in section 3.3.f of the long write-up
!c
      include 'kb3dinc.f90'

      common /ke_prop/ utk(nq),utkc(nq)
      complex utk,utkc
!c
!c     kinetic energy propagator
!c
      do 101 iq=1,nq
      i=ix(iq)
      j=jx(iq)
      k=kx(iq)
      ip2=i**2+j**2+k**2
      if (ip2 .gt. idim**2) then
        utk(iq)=cmplx(1.,0.)
        utkc(iq)=cmplx(0.,0.)
      else
        p2=dpz**2*float(ip2)
        ep=h2m2*p2
        utk(iq)=cexp(eye*dt*ep)
        if (abs(ep*dt) .le. 1.e-5) then
          utkc(iq)=eye*dt
        else
          utkc(iq)=(utk(iq)-1.)/ep
        end if
      end if
 101  continue
!c
      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tstepi(nt,ist)
!c
!c     Time-stepping the Green functions given the collision 
!c     integrals calculated in COLLIS.  the subroutine uses eqs. (15)-(17) of
!c     the long write-up.
!c
      include 'kb3dinc.f90'
!c
      common /ke_prop/ utk(nq),utkc(nq)
      common /gre_inc/ dgla(nq,ntl),dgle(nq,ntl)
      common /step12/ dgla0(nq,ntl),dgle0(nq,ntl),diag(nq)

      complex dgla,dgle,dgla0,dgle0,diag,utk,utkc
!c
      if (ist .eq. 1) then
!c
      do 6 it=1,ntl
      do 6 iq=1,nq
      dgle0(iq,it)=cmplx(0.,0.)
      dgla0(iq,it)=cmplx(0.,0.)
   6  continue
      do 7 iq=1,nq
      diag(iq)=cmplx(0.,0.)
   7  continue
!c
      call collis(nt)
!c
      do 11 it=1,nt
      do 11 iq=1,nq
      dgle0(iq,it)=dgle(iq,it)
      dgla0(iq,it)=dgla(iq,it)
  11  continue
      do 12 iq=1,nq
      diag(iq)=dgla(iq,nt)-dgle(iq,nt)
  12  continue
!c
      else if (ist .eq. 2) then
!c
      call collis(nt+1)
!c
      do 21 it=1,nt
      do 21 iq=1,nq
      dgle0(iq,it)=.5*(dgle0(iq,it)+dgle(iq,it))
      dgla0(iq,it)=.5*(dgla0(iq,it)+dgla(iq,it))
  21  continue
      do 22 iq=1,nq
      diag(iq)=.5*(diag(iq)+dgla(iq,nt+1)-dgle(iq,nt+1))
  22  continue
!c
      endif
!c
      do 31 it=1,nt
      do 32 iq=1,nq
      gre(iq,it,nt+1)=gre(iq,it,nt)*utk(iq)+utkc(iq)*dgle0(iq,it)
      gre(iq,nt+1,it)=gre(iq,nt,it)*conjg(utk(iq))+ &
                      & conjg(utkc(iq))*dgla0(iq,it)
  32  continue
  31  continue
      do 33 iq=1,nq
      gre(iq,nt+1,nt)=(-eye+gre(iq,nt,nt))*conjg(utk(iq))+ &
                     & conjg(utkc(iq))*dgla0(iq,nt)
      gre(iq,nt+1,nt+1)=gre(iq,nt,nt)-eye*dt*real(diag(iq))
!c
!c     dgla-dgle picks up an imaginary part in loop 303 in COLLIS
!c     because of the addition of the hermitian mean field which
!c     breaks the symmetry between I< and I>.  Taking the real
!c     part of diag corrects this.
!c
  33  continue
!c
      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine collis(mt)
!c
!c     This subroutine performs the time integration in the 
!c     calculation of the Kadanoff-Baym collision integrals
!c     given the sig<,> and G<,>. more details are given in section 3.3.g of the
!c     long write-up.
!c
!c     I> (mt,itp), itp=1,mt  : dgla
!c     I'< (it,mt), it=1,mt    : dgle
!c     

      include 'kb3dinc.f90'

      dimension tsum(nq)
      dimension gla(nq,ntl),gle(nq,ntl)
      common /gre_inc/ dgla(nq,ntl),dgle(nq,ntl)

      complex tsum,gla,gle,dgla,dgle
!c
      do 1 it=1,ntl
      do 1 iq=1,nq
      dgle(iq,it)=cmplx(0.,0.)
      dgla(iq,it)=cmplx(0.,0.)
   1  continue
!c
      if (mt .eq. 1) return
!c
!c     calculate collision integrals by trapezoidal rule
!c
!c     I'< first
!c
      do 100 it=1,mt
!c
      do 111 itb=1,mt
      if (itb .lt. it) then
        do 112 iq=1,nq
        gle(iq,itb)=-conjg(gre(iq,itb,it))
        gla(iq,itb)=gre(iq,it,itb)
 112    continue
      else if (itb .gt. it) then
        do 113 iq=1,nq
        gle(iq,itb)=gre(iq,it,itb)
        gla(iq,itb)=-conjg(gre(iq,itb,it))
 113    continue
      else if (itb .eq. it) then
        do 114 iq=1,nq
        gle(iq,itb)=gre(iq,it,itb)
        gla(iq,itb)=-eye+gle(iq,itb)
 114    continue
      end if
 111  continue
!c
      do 120 iq=1,nq
      tsum(iq)=cmplx(0.,0.)
 120  continue
!c
      if (it .gt. 1) then
      do 121 iq=1,nq
      tsum(iq)=tsum(iq)+.5*(gla(iq,1)*sgle(iq,1)+ &
                           & gla(iq,it)*sgle(iq,it))
 121  continue
      if (it .gt. 2) then
        do 122 itb=2,it-1
        do 123 iq=1,nq
        tsum(iq)=tsum(iq)+gla(iq,itb)*sgle(iq,itb)
 123    continue
 122    continue
      endif
      endif
!c
      if (mt .gt. 1) then
      do 131 iq=1,nq
      tsum(iq)=tsum(iq)+.5*(gle(iq,1)*conjg(sgla(iq,1))+ &
                     & gle(iq,mt)*conjg(sgla(iq,mt)))
 131  continue
      if (mt .gt. 2) then
        do 132 itb=2,mt-1
        do 133 iq=1,nq
        tsum(iq)=tsum(iq)+gle(iq,itb)*conjg(sgla(iq,itb))
 133    continue
 132    continue
      endif
      endif
!c
      if (it .lt. mt) then
      do 141 iq=1,nq
      tsum(iq)=tsum(iq)+.5*(gle(iq,it)*sgle(iq,it)+ &
                           & gle(iq,mt)*sgle(iq,mt))
 141  continue
      if (it .lt. mt-1) then
        do 142 itb=it+1,mt-1
        do 143 iq=1,nq
        tsum(iq)=tsum(iq)+gle(iq,itb)*sgle(iq,itb)
 143    continue
 142    continue
      endif
      endif
!c
      do 150 iq=1,nq
      dgle(iq,it)=tsum(iq)*dt
 150  continue
!c
 100  continue
!c
!c     end of I' calculation; now I>, t is fixed at mt
!c      
      do 200 itp=1,mt
!c
      do 211 itb=1,mt
      if (itb .lt. itp) then
        do 212 iq=1,nq
        gle(iq,itb)=gre(iq,itb,itp)
        gla(iq,itb)=-conjg(gre(iq,itp,itb))
 212    continue
      else if (itb .gt. itp) then
        do 213 iq=1,nq
        gle(iq,itb)=-conjg(gre(iq,itp,itb))
        gla(iq,itb)=gre(iq,itb,itp)
 213    continue
      else if (itb .eq. itp) then
        do 214 iq=1,nq
        gle(iq,itb)=gre(iq,itp,itb)
        gla(iq,itb)=-eye+gle(iq,itb)
 214    continue
      end if
 211  continue
!c
      do 220 iq=1,nq
      tsum(iq)=cmplx(0.,0.)
 220  continue
!c
      if (itp .gt. 1) then
      do 221 iq=1,nq
      tsum(iq)=tsum(iq)+.5*(gle(iq,1)*sgla(iq,1)+ &
                           & gle(iq,itp)*sgla(iq,itp))
 221  continue
      if (itp .gt. 2) then
        do 222 itb=2,itp-1
        do 223 iq=1,nq
        tsum(iq)=tsum(iq)+gle(iq,itb)*sgla(iq,itb)
 223    continue
 222    continue
      endif
      endif
!c
      if (mt .gt. 1) then
      do 231 iq=1,nq
      tsum(iq)=tsum(iq)+.5*(gla(iq,1)*conjg(sgle(iq,1))+ &
                    & gla(iq,mt)*conjg(sgle(iq,mt)))
 231  continue
      if (mt .gt. 2) then
        do 232 itb=2,mt-1
        do 233 iq=1,nq
        tsum(iq)=tsum(iq)+gla(iq,itb)*conjg(sgle(iq,itb))
 233    continue
 232    continue
      endif
      endif
!c
      if (itp .lt. mt) then
      do 241 iq=1,nq
      tsum(iq)=tsum(iq)+.5*(gla(iq,itp)*sgla(iq,itp)+ &
                           & gla(iq,mt)*sgla(iq,mt))
 241  continue
      if (itp .lt. mt-1) then
        do 242 itb=itp+1,mt-1
        do 243 iq=1,nq
        tsum(iq)=tsum(iq)+gla(iq,itb)*sgla(iq,itb)
 243    continue
 242    continue
      endif
      endif
!c
      do 250 iq=1,nq
      dgla(iq,itp)=tsum(iq)*dt
 250  continue
!c
 200  continue
!c
      do 260 iq=1,nq
      dpot(iq)=eye*conjg(dgle(iq,mt))
 260  continue
!c
!c     if ifk = 1, add mean field contributions
!c
      if (ifk .eq. 1) then     ! meanf
!c
      do 301 it=1,mt
      do 302 iq=1,nq
      dgla(iq,it)=dgla(iq,it)+smf(iq)*gre(iq,mt,it)
      dgle(iq,it)=dgle(iq,it)+gre(iq,it,mt)*smf(iq)
 302  continue
 301  continue
      do 303 iq=1,nq
      dgla(iq,mt)=dgla(iq,mt)-eye*smf(iq)
 303  continue
!c
      end if
!c
      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sig_fk(mt)
!c
!c     Calculates mean field self energies 
!c     using the Fast-Fourier Transform routine CFFT3D from 
!c     the SGI COMPLIB
!c
      include 'kb3dinc.f90'

      dimension wk(nvec,nvec,nvec)
!c
      complex wk

      mid=nvec/2+1
      if (idim .gt. mid) then
        write (*,*) 'Array allocation for FFT too small.'
        stop
      endif
!c
      do 1 iq=1,nq
      smf(iq)=cmplx(0.,0.)
 1    continue
!c
      do 116 if=1,nvec
      do 117 jf=1,nvec
      do 118 kf=1,nvec
      wk(kf,jf,if)=cmplx(0.,0.)
 118  continue
 117  continue
 116  continue
!c
      do 121 i=-idim,idim
        if=i+mid
        if (if .eq. nvec+1) if=1
      do 121 j=-idim,idim
        jf=j+mid
        if (jf .eq. nvec+1) jf=1
      do 121 k=-idim,idim
        kf=k+mid
        if (kf .eq. nvec+1) kf=1
!c
        ip2=i**2+j**2+k**2
        if (ip2 .le. idim**2) then
        iq=irq(i,j,k)
        wk(kf,jf,if)=eye*gre(iq,mt,mt)
        endif
!c
 121  continue
!c
      call cfft3d(-1,nvec,nvec,nvec,wk,nvec,nvec,ftcoef)
!c
      do 301 i=1,nvec
      do 302 j=1,nvec
      do 303 k=1,nvec
      wk(k,j,i)=wk(k,j,i)*pot_fk(k,j,i)
 303  continue
 302  continue
 301  continue
!c
      call cfft3d(1,nvec,nvec,nvec,wk,nvec,nvec,ftcoef)
!c
      do 401 iq=1,nq
      if=ix(iq)+mid
      jf=jx(iq)+mid
      kf=kx(iq)+mid
      i=mod(if+mid-1,nvec)
      j=mod(jf+mid-1,nvec)
      k=mod(kf+mid-1,nvec)
      if (i .eq. 0) i=nvec
      if (j .eq. 0) j=nvec
      if (k .eq. 0) k=nvec
      smf(iq)=wk(k,j,i)*dp3/float(nvec*nvec*nvec)
 401  continue
!c
      return
      end
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sigmas(mt)
!c
!c     Calculates self energies in the second-order direct Born 
!c     approximation using the 3-D Fast-Fourier-Transform
!c     routine CFFT3D from the sgi library or CCFFT3D of the cray library.
!c
      include 'kb3dinc.f90'

      dimension cgla(nvec,nvec,nvec),cgle(nvec,nvec,nvec), &
               & h(nvec,nvec,nvec),h1(nvec,nvec,nvec)
      dimension potl(nvec+1,nvec+1,nvec+1)
!c
      complex cgla,cgle,h,h1
!c
      fact=anu*dp3**2
      scale=float(nvec*nvec*nvec)
      
      mid=nvec/2+1
      if (idim .gt. mid) then
      write (*,*) 'Array allocation for FFT too small.'
      stop
      endif
!c
      do 11 if=1,nvec+1
      do 12 jf=1,nvec+1
      do 13 kf=1,nvec+1
      potl(kf,jf,if)=0.
  13  continue
  12  continue
  11  continue
!c
      do 21 i=-idim,idim
      do 22 j=-idim,idim
      do 23 k=-idim,idim
      potl(k+mid,j+mid,i+mid)=pot(k,j,i)
  23  continue
  22  continue
  21  continue
!c
      do 100 itb=1,mt
!c
      do 111 if=1,nvec
      do 112 jf=1,nvec
      do 113 kf=1,nvec
      cgla(kf,jf,if)=cmplx(0.,0.)
      cgle(kf,jf,if)=cmplx(0.,0.)
 113  continue
 112  continue
 111  continue
!c
      do 121 i=-idim,idim
        if=i+mid
        if (if .eq. nvec+1) if=1
      do 121 j=-idim,idim
        jf=j+mid
        if (jf .eq. nvec+1) jf=1
      do 121 k=-idim,idim
        kf=k+mid
        if (kf .eq. nvec+1) kf=1
!c
        ip2=i**2+j**2+k**2
        if (ip2 .le. idim**2) then
        iq=irq(i,j,k)
        cgla(kf,jf,if)=gre(iq,mt,itb)
        cgle(kf,jf,if)=gre(iq,itb,mt)
        endif
!c
 121  continue
!c
      if (itb .eq. mt) then
        do 122 i=-idim,idim
          if=i+mid
          if (if .eq. nvec+1) if=1
        do 122 j=-idim,idim
          jf=j+mid
          if (jf .eq. nvec+1) jf=1
        do 122 k=-idim,idim
          kf=k+mid
          if (kf .eq. nvec+1) kf=1
!c
          ip2=i**2+j**2+k**2
          if (ip2 .le. idim**2) then
          iq=irq(i,j,k)
          cgla(kf,jf,if)=-eye+gre(iq,mt,itb)
          endif
!c
 122    continue
      endif
!c
      call cfft3d(-1,nvec,nvec,nvec,cgle,nvec,nvec,ftcoef)
      call cfft3d(-1,nvec,nvec,nvec,cgla,nvec,nvec,ftcoef)
!c
      do 131 if=1,nvec      
      i=nvec+2-if
      if (if .eq. 1) i=1
      do 132 jf=1,nvec  
      j=nvec+2-jf
      if (jf .eq. 1) j=1    
      do 133 kf=1,nvec
      k=nvec+2-kf
      if (kf .eq. 1) k=1    
      h1(kf,jf,if)=cgle(kf,jf,if)*cgla(k,j,i)
 133  continue
 132  continue
 131  continue
!c 
      call cfft3d(1,nvec,nvec,nvec,h1,nvec,nvec,ftcoef)
!c
      do 150 i=1,nvec
      do 150 j=1,nvec
      do 150 k=1,nvec
      h1(k,j,i)=h1(k,j,i)/scale
 150  continue
!c
      do 201 if=1,nvec
      i=mod(if+mid-1,nvec)
      if (i .eq. 0) i=nvec
      do 202 jf=1,nvec      
      j=mod(jf+mid-1,nvec)
      if (j .eq. 0) j=nvec
      do 203 kf=1,nvec      
      k=mod(kf+mid-1,nvec)
      if (k .eq. 0) k=nvec
      h(kf,jf,if)=h1(k,j,i)*potl(kf,jf,if)
 203  continue
 202  continue
 201  continue
!c 
      call cfft3d(-1,nvec,nvec,nvec,h,nvec,nvec,ftcoef)
!c
      do 301 if=1,nvec      
      i=nvec+2-if
      if (if .eq. 1) i=1
      do 302 jf=1,nvec      
      j=nvec+2-jf
      if (jf .eq. 1) j=1    
      do 303 kf=1,nvec      
      k=nvec+2-kf
      if (kf .eq. 1) k=1    
      cgle(kf,jf,if)=cgle(kf,jf,if)*h(kf,jf,if)
      cgla(kf,jf,if)=cgla(kf,jf,if)*h(k,j,i)
 303  continue
 302  continue
 301  continue
!c 
      call cfft3d(1,nvec,nvec,nvec,cgle,nvec,nvec,ftcoef)
      call cfft3d(1,nvec,nvec,nvec,cgla,nvec,nvec,ftcoef)
!c
      do 401 iq=1,nq
      if=ix(iq)+mid
      jf=jx(iq)+mid
      kf=kx(iq)+mid
      i=mod(if+mid-1,nvec)
      j=mod(jf+mid-1,nvec)
      k=mod(kf+mid-1,nvec)
      if (i .eq. 0) i=nvec
      if (j .eq. 0) j=nvec
      if (k .eq. 0) k=nvec
      sgle(iq,itb)=fact*cgle(k,j,i)/scale
      sgla(iq,itb)=fact*cgla(k,j,i)/scale
 401  continue
!c
 100  continue      
!c     
      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine potpaw
!c
!c     it calculates the square of the Gaussian potential given by eq. (18) of 
!c     the long write-up and also the Fourier Transform of the Gaussian
!c     potential for use in the calculation of the mean field.
!c
      include 'kb3dinc.f90'
!c
      data eta,vp /0.57,453./
!c
      vpaw=pi**3*eta**6*vp**2
      pawc=exp(-.5*(eta*dpz)**2)
!c
      do 1 i=-idim,idim
      i2=i**2
      do 1 j=-idim,idim
      j2=j**2
      do 1 k=-idim,idim
      k2=k**2
      ip2=k2+j2+i2
      pot(k,j,i)=vpaw*pawc**ip2
   1  continue
!c
      return
      end
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine potnuc_fk
!c
!c     Calculates Fourier Transform of Das Gupta et al's potential
!c     for use in the Fock part of the mean field
!c
      include 'kb3dinc.f90'
!c
      mid=nvec/2+1
      if (idim .gt. mid) then
        write (*,*) 'Array allocation for FFT too small.'
        stop
      endif
!c
      pfz=1.3333
      grho=.16
      cdas=64.95*anu
      flamb=(1.58*pfz)**2
      do 11 if=1,nvec
      do 12 jf=1,nvec
      do 13 kf=1,nvec
      i=if-mid
      j=jf-mid
      k=kf-mid
      ip2=i**2+j**2+k**2
      das=(2.*cdas/grho)/(1.+dpz**2*float(ip2)/flamb)
      pot_fk(kf,jf,if)=cmplx(das,0.)
 13   continue
 12   continue
 11   continue
!c
!c     test
!c
!c      do 500 k=1,nvec
!c      write (*,*) k,pot_fk(mid,mid,k)
!c 500  continue
!c
      call cfft3d(-1,nvec,nvec,nvec,pot_fk,nvec,nvec,ftcoef)
!c
!c     test
!c
!c      do 501 k=1,nvec
!c      write (*,*) k,pot_fk(mid,mid,k)
!c 501  continue
!c
      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine potcoul
!c
!c     this subroutine calculates the screened Coulomb potential given by
!c     eq. (19) of the long write-up.
!c
      include 'kb3dinc.f90'
!c
      akp2=akp**2
!c
      do 1 i=-idim,idim
      do 1 j=-idim,idim
      do 1 k=-idim,idim
      ip2=k**2+j**2+i**2
      fp2=dpz**2*float(ip2)
      pot(k,j,i)=(4.*pi*coul/(fp2+akp2))**2
   1  continue
!c
      return
      end      
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine potcoul_fk
!c
!c     Calculates Fourier Transform of unscreened Coulomb potential
!c     for use in the Fock part of the mean field
!c
      include 'kb3dinc.f90'
!c
      mid=nvec/2+1
      if (idim .gt. mid) then
        write (*,*) 'Array allocation for FFT too small.'
        stop
      endif
!c
      do 11 if=1,nvec
      do 12 jf=1,nvec
      do 13 kf=1,nvec
      i=if-mid
      j=jf-mid
      k=kf-mid
      ip2=i**2+j**2+k**2
      if (ip2 .eq. 0) vcoul=vcd3d/dp3
      if (ip2 .ne. 0) vcoul=4.*pi*coul/(dpz**2*float(ip2))
      pot_fk(kf,jf,if)=cmplx(vcoul,0.)
 13   continue
 12   continue
 11   continue
!c
      call cfft3d(-1,nvec,nvec,nvec,pot_fk,nvec,nvec,ftcoef)
!c
      return
      end
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calakp(nt)
!c
!c     calculates self-consistent static screening constant (kappa)
!c     using eq. (20) of the long write-up.
!c
      include 'kb3dinc.f90'
!c
      akp=.0
!c
      do 1 i=-idim,idim
      do 1 j=-idim,idim
      do 1 k=-idim,idim
      ip2=i**2+j**2+k**2
      if (ip2 .le. idim**2) then
      iq=irq(i,j,k)
      a1=-real(eye*gre(iq,nt,nt))
      p2=float(i**2+j**2+k**2)
      if(p2.eq..0) akp=akp+a1*2.*pi
      if(p2.gt..0) akp=akp+a1/p2
      endif
 1    continue
      akp=akp*dpz/(4.*pi)
      akp=akp*coul*anu/(h2m2*pi)
      akp=sqrt(akp)
!c
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine out_stat(nt)
!c
!c     calculates statistics of the distribution : 
!c     density, quadrupole moment and energies using eqs.(21)-25
!c     of the long write-up.
!c
      include 'kb3dinc.f90'
      complex pe_corr
!c
      fke=0.
      pe_corr=cmplx(0.,0.)
      pe_mf=0.
      densy=0.
      qua=0.
!c
      do 100 i=-idim,idim
      do 100 j=-idim,idim
      do 100 k=-idim,idim
      ip2=i**2+j**2+k**2
      if ( ip2 .le. idim**2 ) then
      iq=irq(i,j,k)
      a1=-real(eye*gre(iq,nt,nt))
      pz2=(float(k)*dpz)**2
      py2=(float(j)*dpz)**2
      px2=(float(i)*dpz)**2
      pa2=py2+px2
      p2=pa2+pz2
      pe_corr=pe_corr+dpot(iq)
      pe_mf=pe_mf-real(eye*gre(iq,nt,nt)*smf(iq))
      fke=fke+p2*a1
      densy=densy+a1
      qua=qua+(2.*pz2-pa2)*a1
      endif
 100  continue
!c
      densy=anu*dp3*densy    
      qua=anu*dp3*sqrt(5./(16.*pi))*qua/densy
      fke=anu*dp3*h2m2*fke
      pe_corr=anu*dp3*.5*pe_corr
      pe_mf=anu*dp3*.5*pe_mf
!c
      if ( isy .eq. 1 ) then 
        omega_pl=sqrt(4.*pi*coul*densy*h2m2*2.)/hbar
        per_pl=2.*pi/omega_pl        
        densy=densy*aindb**3
        fke=fke*aindb**3/eryd
        pe_corr=pe_corr*aindb**3/eryd
        pe_mf=pe_mf*aindb**3/eryd
        qua=qua*aindb**2
      else if ( isy .eq. 2 ) then
        fke=fke/densy
        pe_corr=pe_corr/densy
        pe_mf=pe_mf/densy
      endif
!c
      oke(nt)=fke
      ope_corr(nt)=real(pe_corr)
      ope_mf(nt)=pe_mf
      ote(nt)=fke+ope_corr(nt)+ope_mf(nt)
      oden(nt)=densy
      oqua(nt)=qua
!c
      if ( nplot .gt. 0 ) then

      do 200 ip=1,nplot
      if ( nt .eq. iplot(ip) ) then
      nchn=50+ip
      do 210 i=0,idim      
      iqx=irq(i,0,0)
      iqz=irq(0,0,i)
      fx=-eye*gre(iqx,nt,nt)
      fz=-eye*gre(iqz,nt,nt)
      write (nchn,9990) float(i)*dpz*aindb,fx,fz
 210  continue
      endif
 200  continue

 9990 format (3e13.5)
!c
      endif
!c
      return
      end
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine pwrite(ntimes,time0,istf,mod_store)
!c     
!c      this subroutine writes the input parameters and the output
!c
      include 'kb3dinc.f90'

      if( isy .eq. 1 ) then
      dpz=dpz*aindb
      akp=akp*aindb
      endif
!c 
      write(6,*)
      write(6,*)'              The input parameters of the program:'
      write(6,*)'              ___________________________________'
      write(6,*)
      write(6,10) dpz,dft
 10   format(5x,'dpz   =',e13.6,3x,'dft   =',e13.6)
      write(6,20) ntimes,time0
 20   format(5x,'ntimes=',i3,13x,'time0 =',e13.6)
      write(6,27) 2*idim+1
 27   format(5x,'number momentum mesh points per dimension=',i3)
      
      if( isy .eq. 1) then
      write(6,*)
      write(6,40)aindb
 40   format(5x,'exciton first radius (aindb)=',e13.6)
      write(6,41)eryd
 41   format(5x,'Rydberg in eV (eryd) ='e13.6)
      write(6,42)amee
 42   format(5x,'electron effective mass (amee) =',e13.6)
      write(6,50) epsb
 50   format(5x,'background dielectric constant (epsb) =',e13.6)
      endif
      if ( isy .eq. 1) then
      write(6,*)
      write(6,61)
 61   format(5x,'initial distribution parameters:')      
      write(6,60) ra,rhb,rgam
 60   format(5x,'ra   =',e13.6,3x,'rhb  =',e13.6,3x,'rgam =',e13.6)
      write(6,70) rgam1,rgam2
 70   format(5x,'rgam1=',e13.6,3x,'rgam2=',e13.6)
      endif
!c
      if( isy .eq. 2) then
      write(6,*)
      write(6,81)
 81   format(5x,'initial distribution parameters:')
!c
      write(6,80)p0,pf
 80   format(5x,'p0=',e13.6,3x,'pf=',e13.6)
      endif
      write(6,*)
      write(6,19)
 19   format(5x,'options selected:')
      if( mod_store .eq. 0) then
      write(6,21)mod_store
  21  format('full storage mode,mod_store='i1)
      endif
      if( mod_store .eq. 1) then
      write(6,22)mod_store
  22  format(5x,'full storage in one octant,mod_store='i1)
      endif
      if( mod_store .eq. 2) then
      write(6,23)mod_store
  23  format(5x,'spherical symmetry,mod_store=',i1)
      endif
      if( mod_store .eq. 3) then
      write(6,24)mod_store
  24  format(5x,'cylindrical symmetry,mod_store is set to ',i1)
      end if
      if( istf .eq. 1) then
      write(6,25)istf
  25  format(5x,'two passes per time step,istf= ',i1)
      endif
      if( isy .eq. 1) then
      write(6,26)isy
  26  format(5x,'electrons calculations,isy= ',i1)
      endif     
      if( isy .eq. 2) then
      write(6,28)isy
  28  format(5x,'nuclear calculations, isy= ',i1)
      endif
      if( ifk .eq. 1) then
      write(6,29)ifk
  29  format(5x,'mean field included,ifk= ',i1)
      else
      write(6,30)ifk
  30  format(5x,'mean field not included,ifk= ',i1)
      endif
      if( isy .eq. 1 .and. isck .eq. 1) then
      write(6,71)akp
 71   format(5x,'self-consistent static screening const.=',e13.6)
      write(6,78)isck
 78   format(5x,'the above constant is calculated, isck= ',i1)
      endif
      if( isy .eq. 1 .and. isck .eq. 0) then
      write(6,75)akp
 75   format(5x,'self-consistent static screening const.=',e13.6)
      write(6,72)isck
 72   format(5x,'the above constant is not calculated,isck= ',i1)
      endif
      if( isy .eq. 1) write(6,73) omega_pl
 73   format(5x,'plasma frequency in fs^-1=',e13.6)
!c
      write(6,*)
      write(6,47)
 47   format(5x,'dimension statemement parameters:')
      write(6,13)nvec,idim
 13   format(5x,'nvec=',i5,5x,'idim=',i5)
      write(6,14)nq,ntl
 14   format(5x,'nq  =',i5,5x,'ntl =',i5)
      write(6,*)
      write(6,90)
 90   format(15x,'The output of the program')
      write(6,100)
 100  format(15x,'--------------------------')
      write(6,110)
 110  format(3x,'Time',4x,'K. Energy',3x,'P. Energy',3x,'P. Energy', &
        & 3x,'T. Energy',3x,'Density',2x,'Quad. Moment'/25x,'(MF)',7x, &
       & '(Corr.)')
      do 120 n=1,ntimes
      tt=time0+dft*float(n-1)
      write(6,121)tt,oke(n),ope_mf(n),ope_corr(n),ote(n), &
       & oden(n),oqua(n)
 121  format(f8.3,6(2x,e10.4))
 120  continue
      return
      end
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cfft3di(m1,m2,m3,ftc)
!c
!c     this subroutine initializes subroutine ccfft3d of the cray computer.  if 
!c     this subroutine is not available (or not using the a cray computer),
!c     it must be deleted.
!c
      include 'kb3dinc.f90'
      complex ftc(nvec+15+nvec+15+nvec+15)
      dimension work(2*nvec*nvec*nvec),table(2*nvec+nvec+nvec)
      dimension isys(4)
      isys(1)=3
      isys(2)=0
      isys(3)=0
      isys(4)=0
      call ccfft3d(0,nvec,nvec,nvec,1.,arr,1,1,arr, &
                   & nvec,nvec,table,work,isys)
      return
       end
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cfft3d(index,m1,m2,m3,array,m4,m5,ftc)
!c
!c     this subroutine calls subroutine ccfft3d of the cray computer to perform
!c     the FFT calculations.  if this subroutine is not available 
!c     (or not using the a cray computer), it must be deleted.

      include 'kb3dinc.f90'
      complex ftc(nvec+15+nvec+15+nvec+15)
      complex array(nvec,nvec,nvec)
      dimension work(2*nvec*nvec*nvec),table(2*nvec+nvec+nvec)
      dimension isys(4)
      isys(1)=3
      isys(2)=0
      isys(3)=0
      isys(4)=0
      call ccfft3d(-index,nvec,nvec,nvec,1.0, &
                   & array,nvec,nvec,array,nvec,nvec,table,work,isys) 
      return
       end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
