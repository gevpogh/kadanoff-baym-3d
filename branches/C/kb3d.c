#include<conio.h>
#include<stdio.h>
#include<math.h>
#include 'kb3dinc.c'
#include<complex.h>
#define pi acos(-1)


/*
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                                       !c
!c     Nov. 1998                                                         !c
!c                                                                       !c
!c     Authors: H. Sigurd Kohler, Nai-Hang Kwong                         !c
!c              University of Arizona                                    !c
!c                                                                       !c
!c              Hashim A. Yousif                                         !c
!c              University of Pittsburgh@bradford                        !c
!c                                                                       !c
!c     This code calculates the two-time non-equilibrium Green           !c
!c     functions for a homogeneous fermion system.  Starting             !c
!c     from a given input momentum distribution, the code solves         !c
!c     the Kadanoff-Baym transport equation, the self energy             !c
!c     being approximated at the second-order direct Born level.         !c
!c     Two illustrative systems are included : the electron gas          !c
!c     (modeling the conduction band electrons or valence band           !c
!c     holes in an excited semiconductor) and nuclear matter             !c
!c     (modeling heavy nuclei in collisions).                            !c
!c                                                                       !c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*/

void init_grid(int);
void init_var();
void init_elec();
void init_nucl();
void epprop();
void tstepi(int,int);
void collis(int);
void sig_fk(int);
void sigmas(int);
void potpaw();
void potnuc_fk();
void potcoul();
void potcoul_fk();
void calakp(int);
void out_stat(int);
void pwrite(int,int,int,int);
void cfft3di(int,int,int,complex);
void cfft3d(int,int,int,int,complex,int,int,complex);

int main()
{
    complex eye(0.0,1.0);
    FILE *fp,*fc;
    fp=fopen("kb3d.dat","r");
    fc=fopen("kb3d.out","a");
    
    //conversion of data types required (lines 33 till 49 not added)
        
    if(nplot > 0)
    {}
    
    fclose(fp);
    
    //set up grid labels and initialize common arrays and fft coefficients
    
    init_grid(mod_store);
    init_var();
    
    //set up grid labels and initialize common arrays and fft coefficients
   
    cfft3di(nvec,nvec,nvec,ftcoef);
   
   //set up parameters and initial Green functions
 
      if(isy == 1) 
       init_elec();
      if(isy == 2) 
       init_nucl();
    
   //calculate kinetic energy propagator
    
    epprop();
    sigmas(1);
     
     if(ifk == 1)
      sig_fk(1);
      
/* 

!c     test
!c 
!c      do 500 k=0,idim
!c      iq=irq(0,0,k)
!c      iqx=irq(k,0,0)
!c      write (*,*) k,smf(iq)/h2m2,smf(iqx)/h2m2
!c 500  continue

*/      

    time=time0;
  
  // start time loop
    
 for(nt=1;nt<=ntimes;nt++)
  {
            tstepi(nt,1);
            out_stat(nt);
            sigmas(nt+1);
      if(ifk == 1) 
        sig_fk(nt+1);    

      if(istf == 1) 
      {
      tstepi(nt,2);
      sigmas(nt+1);
      if(ifk == 1) sig_fk(nt+1);    
      }

      if((isy == 1) && (isck == 1)) 
      {
       calakp(nt);
       potcoul();
      }

      time=time+dft
  }
 
 // write the parameters and the output of the program
 
   pwrite(ntimes,time0,istf,mod_store);
 
 }

  
 void init_grid(int mod_store)
  {
    /*
    set up map between storage labels and momentum-space
!c     Cartesian grid points.  more details are given in section 3.3.a of the 
!c     long write-up
!c
!c     mod_store = 0 : full storage mode
!c                 1 : full storage in one octant
!c                 2 : spherical symmetry
!c                 3 : cylindrical symmetry  
   */   
   
      int i,j,k;
      int nfreq(idim*idim+1),labelq(idim*idim+1),nfreqc(2*idim+1,idim*idim+1),labelc(2*idim+1,idim*idim+1);
      int iqn,ip2;
      for(i=-idim;i<=idim;i++)
      {
       for(j=-idim;j<=idim;j++)
        {
        for(k=-idim;k<=idim;k++)
         {
           irq(i,j,k)=0;
         } 
        } 
      } 
      
      //map (i,j,k) -> iq
      
      if(mod_store == 0)
      {
      iqn=0;
      for(i=-idim;i<=idim;i++)
      {
       for(j=-idim;j<=idim;j++)
        {
        for(k=-idim;k<=idim;k++)
        {
          ip2=pow(i,2)+pow(j,2)+pow(k,2);
          if(ip2 <= pow(idim,2)) 
          {
                 iqn=iqn+1;
                 irq(i,j,k)=iqn;
          }
        }}}}
  
      else if(mod_store == 1) then
      {
      iqn=0;
      for(i=0;i<=idim;i++)
      {
        for(j=0;j<=idim;j++)
         {
          for(k=0,k<=idim;k++)
           {
           ip2=pow(i,2)+pow(j,2)+pow(k,2);
           if(ip2 <= pow(idim,2))
           {
                iqn=iqn+1;
                irq(i,j,k)=iqn;
           }
      }}}

      else if(mod_store == 2) then
      {
      for(ip2=0;ip2<=pow(idim,2);ip2++)
      {
      nfreq(ip2)=0;
      }

      for(i=0;i<=idim;i++)
       {
       for(j=0;j<=idim;j++)
        {
        for(k=0;k<=idim;++)
          {
           ip2=pow(i,2)+pow(j,2)+pow(k,2);
           if(ip2 <= pow(idim,2)) 
           {
            nfreq(ip2)=nfreq(ip2)+1;
           }
        }}}

      iqn=0;
      for(ip2=0;ip2<=pow(idim,2);ip2++)
      {
      if (nfreq(ip2) > 0)
      {
      iqn=iqn+1;
      labelq(ip2)=iqn;
      }
      }

      for(i=0;i<=idim;i++)
      {
       for(j=0;j<=idim;j++)
        {
          for(k=0;k<=idim;k++)
          {
            ip2=pow(i,2)+pow(j,2)+pow(k,2);
            if(ip2 <= pow(idim,2)) 
            {
             irq(i,j,k)=labelq(ip2);
            } 
        }}}      
      }
      
      else if(mod_store == 3) 
      {
       for(ipc2=0;ipc2<=pow(idim,2);ipc2++)
       {
         for(k=-idim;k<=idim;k++)
           {
             nfreqc(k,ipc2)=0;
       }}
       
      for(i=0;i<=idim;i++)
      {
      for(j=0;j<=idim;j++)
      {
        ipc2=pow(i,2)+pow(j,2);
        for(k=-idim;k<=idim;k++)
        {
          ip2=ipc2+pow(k,2);
          if(ip2 <= pow(idim,2)) 
          {
           nfreqc(k,ipc2)=nfreqc(k,ipc2)+1;
          }
      }}}

      iqn=0;
      for(ipc2=0;ipc2<=pow(idim,2);ipc2++)
      {
       for(k=-idim;k<=idim;k++)
        {
         if (nfreqc(k,ipc2) > 0) 
          {
           iqn=iqn+1;
           labelc(k,ipc2)=iqn;
           }
      }}

      for(i=0;i<=idim;i++)
      {
        for(j=0;j<=idim;j++)
         {
           ipc2=pow(i,2)+pow(j,2);
         for(k=-idim;k<=idim;k++)
         {
           ip2=ipc2+pow(k,2);
           if(ip2 <= pow(idim,2))
           {
             irq(i,j,k)=labelc(k,ipc2);
             irq(-i,j,k)=irq(i,j,k);
             irq(i,-j,k)=irq(i,j,k);
             irq(-i,-j,k)=irq(i,j,k);
             }
      }}}      
 
  }      

      if(iqn != nq) 
      {
      printf("nq should be %d",iqn);
      printf("it is %d",nq);
      exit();
      }

      if((mod_store == 1) || (mod_store == 2))
      {
        for(i=0;i<=idim;i++)
        {
        for(j=0;j<=idim;j++)
        {
        for(k=0;k<=idim;k++)
        {
          iq=irq(i,j,k);
          irq(i,j,-k)=iq;
          irq(i,-j,k)=iq;
          irq(i,-j,-k)=iq;
          irq(-i,j,k)=iq;
          irq(-i,j,-k)=iq;
          irq(-i,-j,k)=iq;
          irq(-i,-j,-k)=iq;
       }}}
      }

      // iq -> (i,j,k)
      
      for(i=-idim;i<=idim;i++)
       {
        for(j=-idim;j<=idim;j++)
         {
          for(k=-idim;k<=idim;k++)
          {
            ip2=pow(i,2)+pow(j,2)+pow(k,2);
            if(ip2 <= pow(idim,2)) 
            {
             iq=irq(i,j,k)
             ix(iq)=i
             jx(iq)=j
             kx(iq)=k
            }
        }}}

      return;
 }
    
   void init_var()
  {
//     initialize arrays and variables
      
      int it,jt,iq;
      int i,j,k;
      for(it=1;it<=ntl;it++)
      {
       for(jt=1;jt<=ntl;jt++)
        {
         for(iq=1;iq<=nq;iq++)
          {
            gre(iq,it,jt)=cmplx(0.,0.);  //change this line
      }}}

      for(it=1;it<=ntl;it++)
      {
        for(iq=1;iq<=nq;iq++)
         {
           sgle(iq,it)=cmplx(0.,0.);
           sgla(iq,it)=cmplx(0.,0.);
      }}

      for(iq=1;iq<=nq;iq++)
      {
      dpot(iq)=cmplx(0.,0.);
      smf(iq)=cmplx(0.,0.);        //meanf
      }

      for(i=-idim;i<=idim;i++)
      {
      for(j=-idim;j<=idim
      {
       for(k=-idim;k<=idim;k++)
      {
        pot(k,j,i)=0.;
      }}}

      for(i=1;i<=nvec;i++)       // meanf
      {
       for(j=1;j<=nvec;j++)
      {
       for(k=1;k<=nvec;k++)
      {
       pot_fk(k,j,i)=cmplx(0.,0.);
      }}}

      return;
   }

    void init_elec()
    { 
  /*
     set up parameters and initial Green functions for
     the electron gas system (semiconductor conduction band). the initial 
     distribution is given by eq.(12) of the long write-up.
  */   
      

      h2m20=3.80998;         //eV-A^2
      h2m2=h2m20/amee;
      anu=2.0;
      hbar=0.658212;           // eV-fs
      dt=dft/hbar;           //eV^-1
      dpz=dpz/aindb;          //A^-1
      akp=akp/aindb;
      dp3=pow((dpz/(2.*pi)),3);
      coul=14.3997/epsb;      //eV-A
      vcd3d=coul*dpz*(pi+2.*.34773)/(pow(pi,2));    //meanf
      rgam22=pow(rgam,2);
  /*
     initial Green functions
  */
  
      for(iq=1;iq<=nq;iq++)
      {
      i=ix(iq);
      j=jx(iq);
      k=kx(iq);
      pz=dpz*float(k);
      py=dpz*float(j);
      px=dpz*float(i);
      pl2=pow(pz,2)+pow(py,2)+pow(px,2);
      ree=h2m2*pl2;
      rehh=h2m20*(rgam1-2.*rgam2)*pl2;
      relh=h2m20*(rgam1+2.*rgam2)*pl2;
      a1=pow((rhb-ree-rehh),2);
      rgh=exp(-a1/(2.0*rgam22));
      a1=pow((rhb-ree-relh),2);
      rgl=exp(-a1/(2.0*rgam22));

      pl=sqrt(pl2);
      cth=pz/(pl+0.000001);
      cth2=pow(cth,2);
      a1=(1.0-cth2)/2.;
      a2=(1.0+3.0*cth2)/6.;

      f=ra*(rgh*a1+rgl*a2);
      gre(iq,1,1)=eye*f;
      }

      if(isck == 1) 
      {
              calakp(1)
      }
      potcoul();
      potcoul_fk();                 // meanf

      return;
 }

   void init_nucl()
   {
   /*   set up parameters and initial Green functions for
        the nucleon system (symmetric nuclear matter).
        the initial distribution is given by eq.(13) of the long write-up.
   */
     

      h2m2=20.7355;              // MeV-fm^2
      anu=4.0;
      hbarc=197.327;             // MeV-fm
      dt=dft/hbarc;              // MeV^-1
      dp3=pow((dpz/(2.0*pi)),3);
      aindb=1.;

      for(iq=1;iq<=nq;iq++)
      {
      i=ix(iq);
      j=jx(iq);
      k=kx(iq);
      pz=dpz*float(k);
      py=dpz*float(j);
      px=dpz*float(i);
      apz=abs(pz);
      pp=sqrt(pow((apz-p0),2)+pow(py,2)+pow(px,2));
      f=0.0;
      if (pp < pf) 
      f=1.0;
      gre(iq,1,1)=eye*f;
      }

      potpaw();
      potnuc_fk();               //meanf

      return;
   }  

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

          void potnuc_fk()
     {
/*     Calculates Fourier Transform of Das Gupta et al's potential
!c     for use in the Fock part of the mean field
*/

      mid=nvec/2+1;
      if(idim > mid)
      {
        printf("Array allocation for FFT too small.");
        stop
      }

      pfz=1.3333;
      grho=0.16;
      cdas=64.95*anu;
      flamb=pow((1.58*pfz),2);
      for(if=1;if<=nvec;if++)
      {
      for(jf=1;if<=nvec;if++)
      {
      for(kf=1;kf<=nvec;kf++)
      {
      i=if-mid;
      j=jf-mid;
      k=kf-mid;
      ip2=pow(i,2)+pow(j,2)+pow(k,2);
      das=(2.0*cdas/grho)/(1.0+pow(dpz,2)*float(ip2)/flamb);
      pot_fk(kf,jf,if)=cmplx(das,0.0);
      }}}

//     test

/*      do 500 k=1,nvec
        write (*,*) k,pot_fk(mid,mid,k)
   500  continue
*/
      cfft3d(-1,nvec,nvec,nvec,pot_fk,nvec,nvec,ftcoef);

/*     test
!c
!c      do 501 k=1,nvec
!c      write (*,*) k,pot_fk(mid,mid,k)
!c 501  continue
*/
      return;
      }

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

      void cfft3di(m1,m2,m3,ftc)
      {
/*     this subroutine initializes subroutine ccfft3d of the cray computer.  if 
!c     this subroutine is not available (or not using the a cray computer),
!c     it must be deleted.
!*/
      
      complex ftc(nvec+15+nvec+15+nvec+15);
      float work(2*nvec*nvec*nvec),table(2*nvec+nvec+nvec);
      int isys(4);
      isys(1)=3;
      isys(2)=0;
      isys(3)=0;
      isys(4)=0;
      ccfft3d(0,nvec,nvec,nvec,1.,arr,1,1,arr,nvec,nvec,table,work,isys)
      return;
      }

            void cfft3d(index,m1,m2,m3,array,m4,m5,ftc)
      {
/*     this subroutine calls subroutine ccfft3d of the cray computer to perform
!c     the FFT calculations.  if this subroutine is not available 
!c     (or not using the a cray computer), it must be deleted.
*/

      complex ftc(nvec+15+nvec+15+nvec+15);
      complex array(nvec,nvec,nvec);
      
      float work(2*nvec*nvec*nvec),table(2*nvec+nvec+nvec);
      int isys(4);
      isys(1)=3;
      isys(2)=0;
      isys(3)=0;
      isys(4)=0;
      ccfft3d(-index,nvec,nvec,nvec,1.0,array,nvec,nvec,array,nvec,nvec,table,work,isys) 
      return;
      }





