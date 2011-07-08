      void out_stat(nt)
      {
/*     calculates statistics of the distribution : 
!c     density, quadrupole moment and energies using eqs.(21)-25
!c     of the long write-up.
*/
      
      complex pe_corr;
      float fke=0.0;
      pe_corr=cmplx(0.,0.);
      float pe_mf=0.0;
      float densy=0.0;
      float qua=0.0;
      int i,j,k;

      for(i=-idim;i<=idim;i++)
      {
      for(j=-idim;j<=idimj++)
      {
      for(k=-idim;k<=idim;k++)
      {
      ip2=pow(i,2)+pow(j,2)+pow(k,2);
      if( ip2 <= pow(idim,2) ) 
      {
      iq=irq(i,j,k);
      a1=-float(eye*gre(iq,nt,nt));
      pz2=pow((float(k)*dpz),2);
      py2=pow((float(j)*dpz),2);
      px2=pow((float(i)*dpz),2);
      pa2=py2+px2;
      p2=pa2+pz2;
      pe_corr=pe_corr+dpot(iq);
      pe_mf=pe_mf-float(eye*gre(iq,nt,nt)*smf(iq));
      fke=fke+p2*a1;
      densy=densy+a1;
      qua=qua+(2.0*pz2-pa2)*a1;
      }}}}
 
      densy=anu*dp3*densy; 
      qua=anu*dp3*sqrt(5.0/(16.0*pi))*qua/densy;
      fke=anu*dp3*h2m2*fke;
      pe_corr=anu*dp3*.5*pe_corr;
      pe_mf=anu*dp3*.5*pe_mf;

      if(isy == 1) 
      { 
        omega_pl=sqrt(4.*pi*coul*densy*h2m2*2)/hbar;
        per_pl=2.*pi/omega_pl;
        densy=densy*pow(aindb,3);
        fke=fke*pow(aindb,3)/eryd;
        pe_corr=pe_corr*pow(aindb,3)/eryd;
        pe_mf=pe_mf*pow(aindb,3)/eryd;
        qua=qua*pow(aindb,2);
      } 
      else if(isy == 2)
      {
        fke=fke/densy;
        pe_corr=pe_corr/densy;
        pe_mf=pe_mf/densy;
      }

      oke(nt)=fke;
      ope_corr(nt)=float(pe_corr)
      ope_mf(nt)=pe_mf;
      ote(nt)=fke+ope_corr(nt)+ope_mf(nt);
      oden(nt)=densy;
      oqua(nt)=qua;

      if(nplot > 0 ) 
      {
      for(ip=1;ip<=nplot;ip++)
      {
      if(nt == iplot(ip)) 
      {
      nchn=50+ip;
      for(i=0;i<=idim;i++)
      {      
      iqx=irq(i,0,0);
      iqz=irq(0,0,i);
      fx=-eye*gre(iqx,nt,nt);
      fz=-eye*gre(iqz,nt,nt);
      write (nchn,9990) float(i)*dpz*aindb,fx,fz
      }}}

 9990 format (3e13.5)
      }
      
      return;
     }
