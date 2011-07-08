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
