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
