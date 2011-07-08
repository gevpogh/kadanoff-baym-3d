  void tstepi(nt,ist)
 {
  /*   Time-stepping the Green functions given the collision 
      integrals calculated in COLLIS.  the subroutine uses eqs. (15)-(17) of
      the long write-up.
  */
    

      common /ke_prop/ utk(nq),utkc(nq)
      common /gre_inc/ dgla(nq,ntl),dgle(nq,ntl)
      common /step12/ dgla0(nq,ntl),dgle0(nq,ntl),diag(nq)

      complex dgla,dgle,dgla0,dgle0,diag,utk,utkc

      if(ist == 1)
      {
      for(it=1;it<=ntl;it++)
      {
      for(iq=1;iq<=nq;iq++)
      {
      dgle0(iq,it)=cmplx(0.,0.);
      dgla0(iq,it)=cmplx(0.,0.);
      }}
      for(iq=1;iq<=nq;iq++)
      {
      diag(iq)=cmplx(0.,0.);
      }

      collis(nt);

      for(it=1;it<=nt;it++)
      {
      for(iq=1;iq<=nq;iq++)
      {
      dgle0(iq,it)=dgle(iq,it);
      dgla0(iq,it)=dgla(iq,it);
      }}
      
      for(iq=1;iq<=nq;iq++)
      {
      diag(iq)=dgla(iq,nt)-dgle(iq,nt);
      }
      }
      else if(ist == 2) 
      {
      collis(nt+1);

      for(it=1;it<=nt;it++)
      {
      for(iq=1;iq<=nq;iq++)
      {
      dgle0(iq,it)=0.5*(dgle0(iq,it)+dgle(iq,it));
      dgla0(iq,it)=0.5*(dgla0(iq,it)+dgla(iq,it));
      }}
      for(iq=1;iq<=nq;iq++)
      {
      diag(iq)=0.5*(diag(iq)+dgla(iq,nt+1)-dgle(iq,nt+1));
      }

      }

      for(it=1;it<=nt;it++)
      {
      for(iq=1;iq<=nq;iq++)
      {
      gre(iq,it,nt+1)=gre(iq,it,nt)*utk(iq)+utkc(iq)*dgle0(iq,it);
      gre(iq,nt+1,it)=gre(iq,nt,it)*conjg(utk(iq))+ conjg(utkc(iq))*dgla0(iq,it);
      }
      }
  
      for(iq=1;iq<=nq;iq++)
      {
      gre(iq,nt+1,nt)=(-eye+gre(iq,nt,nt))*conjg(utk(iq))+ conjg(utkc(iq))*dgla0(iq,nt);
      gre(iq,nt+1,nt+1)=gre(iq,nt,nt)-eye*dt*real(diag(iq));
  /*
    dgla-dgle picks up an imaginary part in loop 303 in COLLIS
    because of the addition of the hermitian mean field which
    breaks the symmetry between I< and I>.  Taking the real
    part of diag corrects this.
  */
      }

      return;
   }
