!C    KB3DINC.F90
!c     kb3dinc.f90
!c
!c     include file to be used in kb3d.f90
!c
      parameter (ntl=61,nq=4620,idim=21,nvec=60,mplot=10)
      complex gre,ftcoef,sgla,sgle,dpot,eye,smf,pot_fk
      common /system_id/ isy
      common /param/ anu,h2m2,eye,pi
      common /par_elec/ hbar,aindb,eryd,amee,epsb,coul,akp,vcd3d,isck, &
                        & ra,rhb,rgam,rgam1,rgam2
      common /par_nucl/ hbarc,p0,pf
      common /grid_size/  dft,dt,dpz,dp3
      common /grid_label/ ix(nq),jx(nq),kx(nq),&
                          & irq(-idim:idim,-idim:idim,-idim:idim)
      common /green/ gre(nq,ntl,ntl)
      common /fft_coef/ ftcoef(nvec+15+nvec+15+nvec+15)
      common /sigma/ sgla(nq,ntl),sgle(nq,ntl)
      common /sig_mf/ ifk,smf(nq)
      common /pot_sqd/  pot(-idim:idim,-idim:idim,-idim:idim)
      common /pot_mf/ pot_fk(nvec,nvec,nvec)
      common /pot_a/ dpot(nq)
      common /output/ oke(ntl),ope_corr(ntl),ope_mf(ntl),&
                      & ote(ntl),oden(ntl),oqua(ntl),omega_pl
      common /plot/ nplot,iplot(mplot)

