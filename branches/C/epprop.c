      void epprop()
      {
      /* calculate single particle energy propagator.
         more details are given in section 3.3.f of the long write-up
     */
      

      common /ke_prop/ utk(nq),utkc(nq)
      complex utk,utkc
    /*
      kinetic energy propagator
    */
    
      for(iq=1;iq<=nq;iq++)
      {
      i=ix(iq);
      j=jx(iq);
      k=kx(iq);
      ip2=pow(i,2)+pow(j,2)+pow(k,2);
      if (ip2 > pow(idim,2)) 
      {
        utk(iq)=cmplx(1.0,0.0);
        utkc(iq)=cmplx(0.0,0.0);
      }  
      else
       {
        p2=pow(dpz,2)*float(ip2);
        ep=h2m2*p2; 
        utk(iq)=cexp(eye*dt*ep);
           if (abs(ep*dt) < 1.e-5) 
             utkc(iq)=eye*dt;
           else
             utkc(iq)=(utk(iq)-1.)/ep;
        }
      }

      return;
    }
