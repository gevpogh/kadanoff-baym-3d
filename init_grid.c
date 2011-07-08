  
  init_grid(mod_store)
  {
     int i,j,k;
      int nfreq(0:idim**2),labelq(0:idim**2),nfreqc(-idim:idim,0:idim**2),labelc(-idim:idim,0:idim**2);
     
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
    
