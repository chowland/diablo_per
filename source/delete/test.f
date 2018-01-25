      program test


      complex A
      real B

      real u(2,2)
 
      A=(1.0,2.0)
      b=3.0
      A=A*b

      write(*,*) 'A,b,A/b: ',A,b,A*b
      u(:,:)=1;
      write(*,*) 'sum(u): ',sum(u)

      stop
      end

