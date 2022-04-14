      SUBROUTINE GQ11(XMIN,XMAX,N,X,Y,W)
C                                                   -=-=- gq11
      Implicit Double Precision(A-H, O-Z)

      DIMENSION X(1),Y(1),V(6),R(6) 
      DATA V/0.03376524,0.16939531,0.38069041,0.61930959,0.83060469,    ! GAUS 01 
     10.96623476/,R/0.08566225,0.18038079,0.23395697,0.23395697,        ! GAUS 01 
     20.18038079,0.08566225/
      IF(N) 250,250,150 
  150 K=0 
      NSAVE=N 
      D=(XMAX-XMIN)/N 
      XL=XMIN-D 
      DO 200 I=1,N
      XL=XL+D 
      DO 200 J=1,6
      K=K+1 
  200 X(K)=XL+D*V(J)
      RETURN
  250 W=0.
      K=0 
      DO 300 I=1,NSAVE
      DO 300 J=1,6
      K=K+1 
  300 W=W+Y(K)*R(J) 
      W=W*D 
      RETURN
      END 
