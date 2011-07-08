    void potcoul()
   {
/*    this subroutine calculates the screened Coulomb potential given by
!c     eq. (19) of the long write-up. 
*/
      
      int i,j,k;
      akp2=pow(akp,2);
      for(i=-idim;i<=idim;i++)
      {
      for(j=-idim;j<=idim;j++)
      {
      for(k=-idim;k<=idim;k++)
      {
      ip2=pow(k,2)+pow(j,2)+pow(i,2);
      fp2=pow(dpz,2)*float(ip2);
      pot(k,j,i)=pow((4.0*pi*coul/(fp2+akp2)),2);
      }}}

      return;
   }     
