    init_var()
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
