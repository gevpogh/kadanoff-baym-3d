  void collis(mt)
  {
/*       This subroutine performs the time integration in the 
       calculation of the Kadanoff-Baym collision integrals
       given the sig<,> and G<,>. more details are given in section 3.3.g of the
      long write-up.
*/    
//     I> (mt,itp), itp=1,mt  : dgla
//     I'< (it,mt), it=1,mt    : dgle
     
      dimension tsum(nq)
      dimension gla(nq,ntl),gle(nq,ntl)
      common /gre_inc/ dgla(nq,ntl),dgle(nq,ntl)

      complex tsum,gla,gle,dgla,dgle

      for(it=1;it<=ntl;it++)
      {
      for(iq=1;iq<=nq;iq++)
      {
      dgle(iq,it)=cmplx(0.,0.);
      dgla(iq,it)=cmplx(0.,0.);
      }}

      if(mt == 1)
       return;
/*
     calculate collision integrals by trapezoidal rule
     I'< first
*/

      for(it=1;it<=mt;it++)
      {
      for(itb=1;itb<=mt;itb++)
      {
      if(itb < it) 
      {
        for(iq=1;iq<=nq;iq++)
        {
        gle(iq,itb)=-conjg(gre(iq,itb,it));
        gla(iq,itb)=gre(iq,it,itb);
        }
      }
      else if(itb > it)
      {
        for(iq=1;iq<=nq;iq++)
        {
        gle(iq,itb)=gre(iq,it,itb);
        gla(iq,itb)=-conjg(gre(iq,itb,it));
        }
      } 
      else if(itb == it) 
      {
        for(iq=1;iq<=nq;iq++)
        {
        gle(iq,itb)=gre(iq,it,itb)
        gla(iq,itb)=-eye+gle(iq,itb)
        }
      }  
      }

      for(iq=1;iq<=nq;iq++)
      {
      tsum(iq)=cmplx(0.,0.);
      }

      if(it > 1) 
      {
      for(iq=1;iq<=nq;iq++)
      {
      tsum(iq)=tsum(iq)+.5*(gla(iq,1)*sgle(iq,1)+ gla(iq,it)*sgle(iq,it))
      }
      if (it > 2) 
      {
        for(itb=2;itb<=it-1;itb++)
        {
        for(iq=1;iq<=nq;iq++)
        {
        tsum(iq)=tsum(iq)+gla(iq,itb)*sgle(iq,itb);
        }}
      }
      }

      if(mt > 1) 
      {
      for(iq=1;iq<=nq;iq++)
      {
      tsum(iq)=tsum(iq)+.5*(gle(iq,1)*conjg(sgla(iq,1))+ gle(iq,mt)*conjg(sgla(iq,mt)))
      }
      if (mt > 2) 
      {
        for(itb=2;itb<=mt-1;itb++)
        {
        for(iq=1;iq<=nq;iq++)
        {
        tsum(iq)=tsum(iq)+gle(iq,itb)*conjg(sgla(iq,itb));
        }}
      }
      }

      if(it < mt) 
      {
      for(iq=1;iq<=nq;iq++)
      {
      tsum(iq)=tsum(iq)+.5*(gle(iq,it)*sgle(iq,it)+ gle(iq,mt)*sgle(iq,mt));
      }
 
      if (it < mt-1) 
      {
        for(itb=it+1;itb<=mt-1;itb++)
        {
        for(iq=1;iq<=nq;iq++)
        {
        tsum(iq)=tsum(iq)+gle(iq,itb)*sgle(iq,itb);
        }}
      }
      }

      for(iq=1;iq<=nq;iq++)
      {
      dgle(iq,it)=tsum(iq)*dt;
      }
   }
/*
    end of I' calculation; now I>, t is fixed at mt
*/      

      for(itp=1;ipt<=mt;ipt++)
      {
      for(itb=1;ibt<=mt;ibt++)
      {
      if(itb < itp) 
      {
        for(iq=1;iq<=nq;iq++)
        {
        gle(iq,itb)=gre(iq,itb,itp);
        gla(iq,itb)=-conjg(gre(iq,itp,itb));
        }
      } 
      else if(itb > itp) 
      {
        for(iq=1;iq<=nq;iq++)
        {
        gle(iq,itb)=-conjg(gre(iq,itp,itb));
        gla(iq,itb)=gre(iq,itb,itp);
        }
      }
      else if(itb == itp) 
      {
        for(iq=1;iq<=nq;iq++)
        {
        gle(iq,itb)=gre(iq,itp,itb);
        gla(iq,itb)=-eye+gle(iq,itb);
        }
      }
      }

      for(iq=1;iq<=nq;iq++)
      {
      tsum(iq)=cmplx(0.,0.);
      }

      if(itp > 1) 
      {
      for(iq=1;iq<=nq;iq++)
      {
      tsum(iq)=tsum(iq)+.5*(gle(iq,1)*sgla(iq,1)+ gle(iq,itp)*sgla(iq,itp));
      }
 
      if(itp > 2) 
      {
      for(itb=2;itb<=itp-1;itb++)
      {
        for(iq=1;iq<=nq;iq++)
        {
        tsum(iq)=tsum(iq)+gle(iq,itb)*sgla(iq,itb);
        }
      }  
      }
      }
      }

      if(mt > 1) 
      {
      for(iq=1;iq<=nq;iq++)
      {
      tsum(iq)=tsum(iq)+.5*(gla(iq,1)*conjg(sgle(iq,1))+ gla(iq,mt)*conjg(sgle(iq,mt)))
      }
      if(mt > 2) 
      {
        for(itb=2;itb<=mt-1;itb++)
        {
        for(iq=1;iq<=nq;iq++)
        {
        tsum(iq)=tsum(iq)+gla(iq,itb)*conjg(sgle(iq,itb));
        }}
      }
      }

      if(itp < mt) 
      {
      for(iq=1;iq<=nq;iq++)
      {
      tsum(iq)=tsum(iq)+.5*(gla(iq,itp)*sgla(iq,itp)+ gla(iq,mt)*sgla(iq,mt));
      }
      if (itp < mt-1) 
      {
        for(itb=itp+1;itb<=mt-1;itb++)
        {
        for(iq=1;iq<=nq;iq++)
        {
          tsum(iq)=tsum(iq)+gla(iq,itb)*sgla(iq,itb);
        }}  
      }
      }

      for(iq=1;iq<=nq;iq++)
      {
      dgla(iq,itp)=tsum(iq)*dt;
      }
   }

      for(iq=1;iq<=nq;iq++)
      {
        dpot(iq)=eye*conjg(dgle(iq,mt));
      }
/*
!c     if ifk = 1, add mean field contributions
*/
      if (ifk == 1)        //meanf
      {
      for(it=1;it<=mt;it++)
      {
      for(iq=1;iq<=nq;iq++)
      {
      dgla(iq,it)=dgla(iq,it)+smf(iq)*gre(iq,mt,it);
      dgle(iq,it)=dgle(iq,it)+gre(iq,it,mt)*smf(iq);
      }}
      
      for(iq=1;iq<=nq;iq++)
      {
      dgla(iq,mt)=dgla(iq,mt)-eye*smf(iq);
      }
      }

      return;
  }

