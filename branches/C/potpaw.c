     void potpaw()
     {
/*     it calculates the square of the Gaussian potential given by eq. (18) of 
!c     the long write-up and also the Fourier Transform of the Gaussian
!c     potential for use in the calculation of the mean field.
*/
      int i,j,k;
       float eta,vp;
       eta=0.57;
       vp=453.0;

      vpaw=pow(pi,3)*pow(eta,6)*pow(vp,2);
      pawc=exp(-.5*pow((eta*dpz),2));

      for(i=-idim;i<=idim;i++)
      {
      i2=pow(i,2);
      for(j=-idim;j<=idim;j++)
      {
      j2=pow(j,2);
      for(k=-idim;k<=idim;k++)
      {
      k2=pow(k,2);
      ip2=k2+j2+i2;
      pot(k,j,i)=vpaw*pow(pawc,ip2);
      }}}
      return;
     }      

