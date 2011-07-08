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
