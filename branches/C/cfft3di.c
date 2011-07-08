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

