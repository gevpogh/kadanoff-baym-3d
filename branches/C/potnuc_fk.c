     void potnuc_fk()
     {
/*     Calculates Fourier Transform of Das Gupta et al's potential
!c     for use in the Fock part of the mean field
*/

      mid=nvec/2+1;
      if(idim > mid)
      {
        printf("Array allocation for FFT too small.");
        stop
      }

      pfz=1.3333;
      grho=0.16;
      cdas=64.95*anu;
      flamb=pow((1.58*pfz),2);
      for(if=1;if<=nvec;if++)
      {
      for(jf=1;if<=nvec;if++)
      {
      for(kf=1;kf<=nvec;kf++)
      {
      i=if-mid;
      j=jf-mid;
      k=kf-mid;
      ip2=pow(i,2)+pow(j,2)+pow(k,2);
      das=(2.0*cdas/grho)/(1.0+pow(dpz,2)*float(ip2)/flamb);
      pot_fk(kf,jf,if)=cmplx(das,0.0);
      }}}

//     test

/*      do 500 k=1,nvec
        write (*,*) k,pot_fk(mid,mid,k)
   500  continue
*/
      cfft3d(-1,nvec,nvec,nvec,pot_fk,nvec,nvec,ftcoef);

/*     test
!c
!c      do 501 k=1,nvec
!c      write (*,*) k,pot_fk(mid,mid,k)
!c 501  continue
*/
      return;
      }
