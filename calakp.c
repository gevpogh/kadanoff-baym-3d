     void calakp(nt)
     {
/*     calculates self-consistent static screening constant (kappa)
!c     using eq. (20) of the long write-up.
*/
      float akp=0.0;
      int i,j,k;
      for(i=-idim;i<=idim;i++)
      {
      for(j=-idim;j<=idim;j++)
      {
      for(k=-idim;k<=idim;k++)
      {
      ip2=pow(i,2)+pow(j,2)+pow(k,2);
      if(ip2 <= pow(idim,2)) 
      {
      iq=irq(i,j,k);
      a1=-float(eye*gre(iq,nt,nt));
      p2=float(pow(i,2)+pow(j,2)+pow(k,2));
      if(p2 == 0.0) akp=akp+a1*2.0*pi;
      if(p2 == 0.0) akp=akp+a1/p2;
      }
      }}}
      akp=akp*dpz/(4.0*pi);
      akp=akp*coul*anu/(h2m2*pi);
      akp=sqrt(akp);
      return;
     }
