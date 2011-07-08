     void potcoul_fk()
     {
/*    Calculates Fourier Transform of unscreened Coulomb potential
!c     for use in the Fock part of the mean field
*/
       
      int jf,kf,if;
      mid=nvec/2+1;
      if(idim > mid) 
      {
        printf("Array allocation for FFT too small.");
        stop
      }

      for(if=1;if<=nvec;if++)
      {
      for(jf=1;jf<=nvec;jf++)
      
      {
      for(kf=1;kf<=nvec;kf++)
      {
      i=if-mid;
      j=jf-mid;
      k=kf-mid;
      ip2=pow(i,2)+pow(j,2)+pow(k,2);
      if(ip2 == 0) vcoul=vcd3d/dp3;
      if(ip2 != 0) vcoul=4.0*pi*coul/(pow(dpz,2)*float(ip2));
      pot_fk(kf,jf,if)=cmplx(vcoul,0.0);
      }}}

      call cfft3d(-1,nvec,nvec,nvec,pot_fk,nvec,nvec,ftcoef);

      return;
     }
