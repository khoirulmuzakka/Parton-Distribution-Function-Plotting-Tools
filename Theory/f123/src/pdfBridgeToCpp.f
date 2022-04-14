C============================================================================= 
      FUNCTION PDF(iSet, IHADRONin, IPARTONin, x, q, iRet)
C============================================================================= 
C      Front end program for  cteq4 pdf's: Fred Olness  6/20/96
C      Modify for CTEQ6.6: Fred Olness  6/2/11
C
C============================================================================= 
      Implicit Double Precision (A-H, O-Z)
      character namein*64
      character name*64
c      character(len=:), allocatable :: namein
      double precision f(-6:6)

      DATA IPARM,IPRINT /1,1/
      SAVE IPARM,IPRINT,IHADRON
      common /PDFNAME/namein
      LOGICAL IFIRST 
      data IFIRST,ICTEQ /.TRUE.,6/
      SAVE IFIRST,ICTEQ
      integer i
      integer temp
C============================================================================= 
c F. Intercept this and set it to LHAPDF by default for now
      ICTEQ=-1
c END intercept      
C============================================================================= 
      IHADRON=IHADRONin
      IPARTON=IPARTONin
C============================================================================= 
      PDF=0.0D0
      IF((X.GE.1.0).OR.(X.LE.0.0)) RETURN
C============================================================================= 

      If((IHADRON.eq.1).or.(IHADRON.eq.0)) then
      else
        write(6,*) "  !*** Only setup for PROTON/NEUTRON "
        STOP
      endif

C============================================================================= 
c The notation of F123 is the same as cteq i.e. 1->u, 2->d so we need to convert it

      temp = IPARTON  
      if(IPARTON.eq.-2) temp=-1
      if(IPARTON.eq.-1) temp=-2
      if(IPARTON.eq. 1) temp=2
      if(IPARTON.eq. 2) temp=1
      temp = temp + 6! need to move it into 0-12
      pdf = PDFDISF123LFortranWrapperHook(temp, x, q)/x
C============================================================================= 
C============================================================================= 
      IPRINT=0
      IFIRST=.FALSE.
      RETURN
      END
C============================================================================= 
C============================================================================= 
