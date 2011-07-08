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

