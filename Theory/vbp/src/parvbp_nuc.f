c ... 03.05.2007 (IJS/JYY)
c ... parvbp_nuc.f

      SUBROUTINE ParVbp_nuc(Iact, Name, Value, Iret) 
C                                                   
C                                           Simplified set/get parameters
C                                               Iret is a true dummy
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER NAME*(*), Uname*10
      Logical LL, LM

      COMMON
     1 / STCOM_NUC / RTS, Q, AMU, NBM, NTG, LSET1, LSET2, iOrd, iBsn
     1 / LOVBP_NUC / LM, LL, Ischm, Iscal
     1 / F0COM_NUC / XA, XB, FA0(-10:10), FB0(-10:10), IA, IB
     1 / INTCM_NUC / AER, RER, ZRO, ONE
     1 /E866xf/ ixfx

      CALL UPC (NAME, Ln, Uname)

      If     (Iact .Eq. 1) Then
         If     (Uname(1:Ln) .Eq. 'RS') Then
            RTS = Value
         ElseIf (Uname(1:Ln) .Eq. 'Q') Then
            Q = Value
            If (.Not. LL) Amu = Q             
         ElseIf (Uname(1:Ln) .Eq. 'AMU') Then
            Amu = Value
            QMULOG = LOG (Q/AMU)
            LL = ABS(QMULOG) .GE. 1.E-1
         ElseIf (Uname(1:Ln) .Eq. 'RERR') Then
            RER = Value
         ElseIf (Uname(1:Ln) .Eq. 'IHD1') Then
            NBM = Nint(Value)
         ElseIf (Uname(1:Ln) .Eq. 'IHD2') Then
            NTG = Nint(Value)
         ElseIf (Uname(1:Ln) .Eq. 'ISET1') Then
            Lset1= Nint(Value)
c            Call Setpdf(Lset1) ! Is done automatically now
         ElseIf (Uname(1:Ln) .Eq. 'ISET2') Then
            Lset2= Nint(Value)
c            Call Setpdf(Lset2) ! Is done automatically now
         ElseIf (Uname(1:Ln) .Eq. 'IORD') Then
            iOrd= Nint(Value)
         ElseIf (Uname(1:Ln) .Eq. 'IBSN') Then
            iBsn= Nint(Value)
         ElseIf (Uname(1:Ln) .Eq. 'ISCHM') Then
            iSchm= Nint(Value)
            If (iSchm .Eq. 0) Then
               LM = .True.
            ElseIf (iSchm .Eq. 1) Then
               LM = .False.
            Else
               Stop  'ISCHM can only be 0 (msbar) or 1 (DIS).'
            EndIf
         ElseIf (Uname(1:Ln) .Eq. 'ISCAL') Then
            iScal= Nint(Value)
            If (iScal .Eq. 0) Then
C                                       if LL=.False., must set Amu=Q
               LL = .False.
               Amu = Q
            ElseIf (iScal .Eq. 1) Then
C                                     if LL=.True., let Amu be an indep. var.
               LL = .True.
            Else
               Stop  'iScal can only be 0 (Amu=Q) or 1 (Amu.Ne.Q).'
            EndIf
         ElseIf (Uname(1:Ln) .Eq. 'IXFX') Then
            IXFX = Nint(Value)
         Else
            Print '(/A/ 4A7 /8A7)',  
     >      'Legal Vbp parameter names are:','RS','Q','AMU','RERR'
     >      ,'IHD1','IHD2','ISET1','ISET2','ISCHM','ISCAL','IORD'
     #      ,'IBSN','IXFX'
            Stop
         EndIf
      ElseIf (Iact .Eq. 2) Then
         If     (Uname(1:Ln) .Eq. 'RS') Then
            Value = RTS
         ElseIf (Uname(1:Ln) .Eq. 'Q') Then
            Value = Q 
         ElseIf (Uname(1:Ln) .Eq. 'AMU') Then
            Value = Amu 
         ElseIf (Uname(1:Ln) .Eq. 'RERR') Then
            Value = RER 
         ElseIf (Uname(1:Ln) .Eq. 'IHD1') Then
            Value = NBM 
         ElseIf (Uname(1:Ln) .Eq. 'IHD2') Then
            Value = NTG 
         ElseIf (Uname(1:Ln) .Eq. 'ISET1') Then
            Value = Lset1
         ElseIf (Uname(1:Ln) .Eq. 'ISET2') Then
            Value = Lset2
         ElseIf (Uname(1:Ln) .Eq. 'IORD') Then
            Value = iOrd
         ElseIf (Uname(1:Ln) .Eq. 'IBSN') Then
            Value = iBsn
         ElseIf (Uname(1:Ln) .Eq. 'ISCHM') Then
            If (LM) Then
              Value = 0
            Else
              Value = 1
            EndIf
         ElseIf (Uname(1:Ln) .Eq. 'ISCAL') Then
            If (LL) Then
              Value = 1
            Else
              Value = 0
            EndIf
         ElseIf (Uname(1:Ln) .Eq. 'IXFX') Then
            Value = IXFX
         Else
            Print '(/A/ 4A7 /8A7)',  
     >      'Legal Vbp parameter names are:','RS','Q','AMU','RERR'
     >      ,'IHD1','IHD2','ISET1','ISET2','ISCHM','ISCAL','IORD'
     #      ,'IBSN','IXFX'
            Stop
         EndIf

      ElseIf (Iact .Eq. 4) Then
         Print '(/A/ 4A7 /8A7)',  
     >   'Vbp parameter names are:','RS','Q','AMU','RERR'
     >   ,'IHD1','IHD2','ISET1','ISET2','ISCHM','ISCAL','IORD'
     #   ,'IBSN','IXFX'
         Print '(/A/4(1pE11.3)/8I7/)',' Their current values are: '
     >   ,RTS, Q, AMU, RER, NBM, NTG, LSET1, LSET2, iSchm, iScal,iOrd
     #   ,iBsn, IXFX
         Print *, 'The two logical switches LM,LL are:', LM, LL   

      Else
         Print *, 'You used Iset =', Iact
     >      , 'The only legal values for Iact in ParVbp_nuc are 1/2/4.'
         Stop
      Endif

      Return
C                         -----------------
      End
