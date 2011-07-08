     void sig_fk(mt)
    {

/*     Calculates mean field self energies 
       using the Fast-Fourier Transform routine CFFT3D from 
       the SGI COMPLIB
*/
      
      dimension wk(nvec,nvec,nvec)
      complex wk

      mid=nvec/2+1;
      if (idim > mid) 
      {
        printf("Array allocation for FFT too small.");
        exit();
      }

      for(iq=1;iq<=nq;iq++)
      {
      smf(iq)=cmplx(0.,0.);
      }

      for(if=1;if<=nvec;if++)
      {
      for(jf=1;jf<=nvec;jf++)
      {
      for(kf=1;kf<=nvec;kf++)
      {
      wk(kf,jf,if)=cmplx(0.,0.);
      }}}

      for(i=-idim;i<=idim;i++)
      {
        if=i+mid;
        if(if == nvec+1) if=1;
      for(j=-idim;j<=idim;j++)
      {
        jf=j+mid;
        if(jf == nvec+1) jf=1;
      for(k=-idim;k<=idim;k++)
      {
        kf=k+mid;
        if(kf == nvec+1) kf=1;
        ip2=pow(i,2)+pow(j,2)+pow(k,2);
        if(ip2 <= pow(idim,2))
        {
        iq=irq(i,j,k);
        wk(kf,jf,if)=eye*gre(iq,mt,mt);
        }
      }

      cfft3d(-1,nvec,nvec,nvec,wk,nvec,nvec,ftcoef);

      for(i=1;i<=nvec;i++)
      {
      for(i=1;i<=nvec;i++)
      {
      for(i=1;i<=nvec;i++)
      {
      wk(k,j,i)=wk(k,j,i)*pot_fk(k,j,i);
      }}}

      cfft3d(1,nvec,nvec,nvec,wk,nvec,nvec,ftcoef);

      for(iq=1;iq<=nq;iq++)
      {
      if=ix(iq)+mid;
      jf=jx(iq)+mid;
      kf=kx(iq)+mid;
      i=mod(if+mid-1,nvec);
      j=mod(jf+mid-1,nvec);
      k=mod(kf+mid-1,nvec);
      if(i == 0) i=nvec;
      if(j == 0) j=nvec;
      if(k == 0) k=nvec;
      smf(iq)=wk(k,j,i)*dp3/float(nvec*nvec*nvec);
      }
      return();
      
   }  
