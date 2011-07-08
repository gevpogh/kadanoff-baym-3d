 void sigmas(mt)
 {
/*     Calculates self energies in the second-order direct Born 
       approximation using the 3-D Fast-Fourier-Transform
       routine CFFT3D from the sgi library or CCFFT3D of the cray library.
*/
      
      float cgla(nvec,nvec,nvec),cgle(nvec,nvec,nvec), h(nvec,nvec,nvec),h1(nvec,nvec,nvec);
      float potl(nvec+1,nvec+1,nvec+1);
      int if,jf,kf;
      complex cgla,cgle,h,h1

      fact=anu*pow(dp3,2);
      scale=float(nvec*nvec*nvec);
      
      mid=nvec/2+1;
      if(idim > mid)
      {
      printf("Array allocation for FFT too small.");
      exit();
      }

      for(if=1;if<=nvec+1;if++)
      {
      for(jf=1;jf<=nvec+1;jf++)
      {
      for(kf=1;kf<=nvec+1;kf++)
      {
      potl(kf,jf,if)=0.0;
      }}}

      for(i=-idim;i<=idim;i++)
      {
      for(j=-idim;j<=idim;j++)
      {
      for(k=-idim;k<=idim;k++)
      {
      potl(k+mid,j+mid,i+mid)=pot(k,j,i);
      }}}
  
      for(itb=1;itb<=mt;itb++)
      {
      for(if=1;if<=nvec;if++)
      {
      for(jf=1;jf<=nvec;jf++)
      {
      for(kf=1;kf<=nvec;kf++)
      {
      cgla(kf,jf,if)=cmplx(0.0,0.0);
      cgle(kf,jf,if)=cmplx(0.0,0.0);
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
        if (ip2 <= pow(idim,2)) 
        {
        iq=irq(i,j,k);
        cgla(kf,jf,if)=gre(iq,mt,itb);
        cgle(kf,jf,if)=gre(iq,itb,mt);
        }
        }}}

      if(itb == mt) 
       {
        for(i=-idim;i<=idim;i++)
        {
          if=i+mid;
          if (if == nvec+1) if=1;
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
          cgla(kf,jf,if)=-eye+gre(iq,mt,itb);
          }
       }}}
        }

       cfft3d(-1,nvec,nvec,nvec,cgle,nvec,nvec,ftcoef);
       cfft3d(-1,nvec,nvec,nvec,cgla,nvec,nvec,ftcoef);

      for(if=1;if<=nvec;if++)
      {      
      i=nvec+2-if;
      if(if == 1) i=1;
      for(jf=1;jf<=nvec;jf++)
      {  
      j=nvec+2-jf;
      if(jf == 1) j=1;   
      for(kf=1;kf<=nvec;kf++)
      {
      k=nvec+2-kf;
      if (kf == 1) k=1;    
      h1(kf,jf,if)=cgle(kf,jf,if)*cgla(k,j,i);
      }}}
      
      cfft3d(1,nvec,nvec,nvec,h1,nvec,nvec,ftcoef);

      for(i=1;i<=nvec;i++)
      {
      for(j=1;j<=nvec;j++)
      {
      for(k=1;k<=nvec;k++)
      {
      h1(k,j,i)=h1(k,j,i)/scale;
      }}}

      for(if=1;if<=nvec;if++)
      {
      i=mod(if+mid-1,nvec);
      if(i == 0) i=nvec;
      for(jf=1;jf<=nvec;jf++)
      {      
      j=mod(jf+mid-1,nvec);
      if(j == 0) j=nvec;
      for(kf=1;kf<=nvec;kf++)
      {      
      k=mod(kf+mid-1,nvec);
      if(k == 0) k=nvec;
      h(kf,jf,if)=h1(k,j,i)*potl(kf,jf,if);
      }}}
 
      cfft3d(-1,nvec,nvec,nvec,h,nvec,nvec,ftcoef);

      for(if=1;if<=nvec;if++)
      {      
      i=nvec+2-if;
      if(if == 1) i=1;
      for(jf=1;jf<=nvec;jf++)
      {     
      j=nvec+2-jf;
      if(jf == 1) j=1;    
      for(kf=1;kf<=nvec;kf++)
      {      
      k=nvec+2-kf;
      if(kf == 1) k=1; 
      cgle(kf,jf,if)=cgle(kf,jf,if)*h(kf,jf,if);
      cgla(kf,jf,if)=cgla(kf,jf,if)*h(k,j,i);
      }}}
 
      cfft3d(1,nvec,nvec,nvec,cgle,nvec,nvec,ftcoef);
      cfft3d(1,nvec,nvec,nvec,cgla,nvec,nvec,ftcoef);

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
      sgle(iq,itb)=fact*cgle(k,j,i)/scale;
      sgla(iq,itb)=fact*cgla(k,j,i)/scale;
      }

  }     
      return;
  }
