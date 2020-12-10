C...A somewhat more realistic example: a study of reconnection
C...effects on the W mass at LEP 2.

C-----------------------------------------------------------------

C...Preamble: declarations.
 
C...All real arithmetic in double precision.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C...Three Pythia functions return integers, so need declaring.
      INTEGER PYK,PYCHGE,PYCOMP

C...EXTERNAL statement links PYDATA on most machines.
      EXTERNAL PYDATA

C...Commonblocks.
C...The event record.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
C...Parameters.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
C...Particle properties + some flavour parameters.
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
C...Decay information.
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
C...Selection of hard scattering subprocesses.
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
C...Parameters. 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
C...Supersymmetry parameters.
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)

C...Local arrays and data.
      DIMENSION NEVT(30),IQQ(4),PMWR(4,4),IJJ(4),DISTSV(4,4),
     &NACC(10,5),PDIFF(5),NACC2(10,5),WACC(10,5),WACC2(10,5)
      DATA NEVT/30*0/,NACC/50*0/,WACC/50*0D0/,NACC2/50*0/,
     &WACC2/50*0D0/

C-----------------------------------------------------------------

C...First section: initialization.

C...CM energy. Number of events. Number of cases. pi0 stable.
      ECM=170D0
      NEV=1000
      NRE=7
      MDCY(PYCOMP(111),1)=0

C...Choice of W+W- process.
      MSEL=0
      MSUB(25)=1  

C...Allow only hadronic W decays.
      DO 100 IDC=MDCY(24,2),MDCY(24,2)+MDCY(24,3)-1
        IF(IABS(KFDP(IDC,1)).GE.6) MDME(IDC,1)=
     &  MIN(0,MDME(IDC,1))
 100  CONTINUE

C...No ISR photon radiation.
      MSTP(11)=0       

C...Initialize.
      CALL PYINIT('CMS','e+','e-',ECM)

C...Set up cluster finding for 4 clusters.
      MSTU(41)=1
      MSTU(47)=4
      PARU(44)=8D0

C...Book histograms for all cases.
      DO 110 IRE=1,NRE
      CALL PYBOOK(IRE,'Prob(number of jets)',20,-0.5D0,19.5D0)
      CALL PYBOOK(IRE+10,'Prob(reconnect) as fn n_jets',20,-0.5D0,
     &19.5D0)
      CALL PYBOOK(IRE+20,'<m_W>(true)',100,70D0,90D0)
      CALL PYBOOK(IRE+30,'<m_W>(LUCLUS-true)',100,-10D0,10D0)
      CALL PYBOOK(IRE+40,'<m_W>(LUCLUS min dist1-true)',100,-10D0,10D0)
      CALL PYBOOK(IRE+50,'<m_W>(LUCLUS min dist2-true)',100,-10D0,10D0)
      CALL PYBOOK(IRE+60,'<m_W>(LUCLUS min dist3-true)',100,-10D0,10D0)
      CALL PYBOOK(IRE+70,'charged multiplicity',80,-1D0,159D0)
  110 CONTINUE

C-----------------------------------------------------------------

C...Second section: event loop.

C...Loop through the cases of interest. See manual for description.
      DO 300 IRE=1,NRE
      IF(IRE.EQ.1) MSTP(115)=0
      IF(IRE.EQ.2) MSTP(115)=1
      IF(IRE.EQ.3) MSTP(115)=2
      IF(IRE.EQ.4) MSTP(115)=3
      IF(IRE.EQ.5) MSTP(115)=5
      IF(IRE.EQ.6) MSTP(115)=11
      IF(IRE.EQ.7) MSTP(115)=12

C...Generate events for each case.
      DO 290 IEV=1,NEV
      CALL PYEVNT

C...List first event of each type.
      IF(IEV.EQ.1) CALL PYLIST(1)

C...Read out whether reconnection occured. W+- positions.
      IRR=MSTI(32)
      IW1=7
      IW2=8

C...Do cluster search. Only keep events with four clusters.
      IACC=0
      CALL PYCLUS(NJET)
      CALL PYFILL(IRE,DBLE(NJET),1D0)
      CALL PYFILL(IRE+10,DBLE(NJET),DBLE(IRR))
      IF(NJET.NE.4) GOTO 290
      NEVT(IRE)=NEVT(IRE)+1

C...Cuts so that all clusters are well separated and have good energy.
      DSPMIN=0.5D0
      EMIN=20D0
      ISEP=1
      DO 120 I1=N+1,N+3
      DO 120 I2=I1+1,N+4
      PA1=SQRT(P(I1,1)**2+P(I1,2)**2+P(I1,3)**2)      
      PA2=SQRT(P(I2,1)**2+P(I2,2)**2+P(I2,3)**2)      
      P12=P(I1,1)*P(I2,1)+P(I1,2)*P(I2,2)+P(I1,3)*P(I2,3)
      DSEP=ACOS(MAX(-1D0,MIN(1D0,P12/(PA1*PA2))))
      DISTSV(I1-N,I2-N)=DSEP
      DISTSV(I2-N,I1-N)=DSEP
  120 IF(DSEP.LT.DSPMIN) ISEP=0
      DO 130 I1=N+1,N+4
  130 IF(P(I1,4).LT.EMIN) ISEP=0
      IF(ISEP.EQ.0) GOTO 290
      NEVT(IRE+10)=NEVT(IRE+10)+1

C...Fill true average W mass.
      PMAVG=0.5D0*(P(IW1,5)+P(IW2,5))
      CALL PYFILL(IRE+20,PMAVG,1D0)

C...Relate clusters to partons.
      IQQ(1)=9
      IQQ(2)=10
      IQQ(3)=11
      IQQ(4)=12
      DO 150 I1=N+1,N+4
      DO 150 I2=1,4
      PMWRS=(P(I1,4)+P(IQQ(I2),4))**2-(P(I1,1)+P(IQQ(I2),1))**2-
     &(P(I1,2)+P(IQQ(I2),2))**2-(P(I1,3)+P(IQQ(I2),3))**2
  150 PMWR(I1-N,I2)=SQRT(MAX(0D0,PMWRS))
      PROMIN=1D10
      DO 160 IJ1=1,4
      DO 160 IJ2=1,4
      DO 160 IJ3=1,4
      DO 160 IJ4=1,4
      IF(IJ1.EQ.IJ2.OR.IJ1.EQ.IJ3.OR.IJ1.EQ.IJ4.OR.IJ2.EQ.IJ3.OR.
     &IJ2.EQ.IJ4.OR.IJ3.EQ.IJ4) GOTO 160
      PRONOW=PMWR(IJ1,1)*PMWR(IJ2,2)*PMWR(IJ3,3)*PMWR(IJ4,4)
      IF(PRONOW.LT.PROMIN) THEN
        PROMIN=PRONOW
        IJJ(1)=N+IJ1
        IJJ(2)=N+IJ2
        IJJ(3)=N+IJ3
        IJJ(4)=N+IJ4
      ENDIF
  160 CONTINUE

C...Fill reconstructed average W mass shifts; under idealized conditions.
      PMWRS=(P(IJJ(1),4)+P(IJJ(2),4))**2-(P(IJJ(1),1)+P(IJJ(2),1))**2-
     &(P(IJJ(1),2)+P(IJJ(2),2))**2-(P(IJJ(1),3)+P(IJJ(2),3))**2
      PMWW1=SQRT(MAX(0D0,PMWRS))
      PMWRS=(P(IJJ(3),4)+P(IJJ(4),4))**2-(P(IJJ(3),1)+P(IJJ(4),1))**2-
     &(P(IJJ(3),2)+P(IJJ(4),2))**2-(P(IJJ(3),3)+P(IJJ(4),3))**2
      PMWW2=SQRT(MAX(0D0,PMWRS))
      PMWWA=0.5D0*(PMWW1+PMWW2)
      CALL PYFILL(IRE+30,PMWWA-PMAVG,1D0)   

C...Now use three different semirealistic schemes to pair jets to W's,
C...designed to pick jet combination which gives masses "closest" to 
C...nominal ones: (1) average mass closest to 80 GeV; (2) minimal
C...deviation of W masses from 80 GeV; (3) jets most back-to-back. 
      PMBEST=0D0
      PMDEVI=80D0
      DMAXI=0D0
      DO 180 ICA=1,3
      IF(ICA.EQ.1) THEN
        I1=N+1
        I2=N+2
        I3=N+3
        I4=N+4
      ELSEIF(ICA.EQ.2) THEN
        I1=N+1
        I2=N+3
        I3=N+2
        I4=N+4
      ELSE
        I1=N+1
        I2=N+4
        I3=N+2
        I4=N+3
      ENDIF
      PMWRS1=(P(I1,4)+P(I2,4))**2-(P(I1,1)+P(I2,1))**2-
     &(P(I1,2)+P(I2,2))**2-(P(I1,3)+P(I2,3))**2
      PMWR1=SQRT(MAX(0D0,PMWRS1))
      PMWRS2=(P(I3,4)+P(I4,4))**2-(P(I3,1)+P(I4,1))**2-
     &(P(I3,2)+P(I4,2))**2-(P(I3,3)+P(I4,3))**2
      PMWR2=SQRT(MAX(0D0,PMWRS2))
      PMWRA=0.5D0*(PMWR1+PMWR2)
      IF(ABS(PMWRA-80D0).LT.ABS(PMBEST-80D0)) PMBEST=PMWRA
      PMDEWW=ABS(PMWR1-80D0)+ABS(PMWR2-80D0)
      IF(PMDEWW.LT.PMDEVI) THEN
        PMDEVI=PMDEWW
        PMBESA=PMWRA
      ENDIF
      DTHE=DISTSV(I1-N,I2-N)+DISTSV(I3-N,I4-N)
      IF(DTHE.GT.DMAXI) THEN
        DMAXI=DTHE
        PMBESB=PMWRA
      ENDIF
  180 CONTINUE

C...Fill average mass shifts for three different choices.
      CALL PYFILL(IRE+40,PMBEST-PMAVG,1D0)
      CALL PYFILL(IRE+50,PMBESA-PMAVG,1D0)
      CALL PYFILL(IRE+60,PMBESB-PMAVG,1D0)

C...Fill numerical statistics as alternative to histograms.
      PDIFF(1)=PMWWA-PMAVG
      PDIFF(2)=PMBEST-PMAVG
      PDIFF(3)=PMBESA-PMAVG
      PDIFF(4)=PMBESB-PMAVG
      DO 200 J=1,4
      IF(ABS(PDIFF(J)).LT.10D0) THEN
        NACC(IRE,J)=NACC(IRE,J)+1
        WACC(IRE,J)=WACC(IRE,J)+PDIFF(J)
      ENDIF
      IF(ABS(PDIFF(J)).LT.4D0) THEN
        NACC2(IRE,J)=NACC2(IRE,J)+1
        WACC2(IRE,J)=WACC2(IRE,J)+PDIFF(J)
      ENDIF
  200 CONTINUE     

C...Charged multiplicity.
      CALL PYEDIT(3)
      CALL PYFILL(IRE+70,DBLE(N),1D0)

C...End event and case loops.
  290 CONTINUE
  300 CONTINUE

C-----------------------------------------------------------------

C...Third section: produce output and end.

C...Cross section table.
      CALL PYSTAT(1)

C...Numerical statistics.
      DO 320 IRE=1,NRE
      DO 320 J=1,4
      WACC(IRE,J)=WACC(IRE,J)/NACC(IRE,J)
      WACC2(IRE,J)=WACC2(IRE,J)/NACC2(IRE,J)
  320 CONTINUE
      WRITE(6,'(A32,I7/5X,7I6)') '# events in total and after cuts:',
     &NEV,(NEVT(IRE+10),IRE=1,NRE)
      DO 330 J=1,4
      WRITE(6,'(I5,7I6,7F10.3)') J,(NACC(IRE,J),IRE=1,NRE),
     &(WACC(IRE,J),IRE=1,NRE)
      WRITE(6,'(57X,6F10.3)') (WACC(IRE,J)-WACC(1,J),IRE=2,NRE)
      WRITE(6,'(I5,7I6,7F10.3)') J,(NACC2(IRE,J),IRE=1,NRE),
     &(WACC2(IRE,J),IRE=1,NRE)
      WRITE(6,'(57X,6F10.3)') (WACC2(IRE,J)-WACC2(1,J),IRE=2,NRE)
  330 CONTINUE

C...Rescale and print histograms.
C...Note: in a realistic run, the histrograms would be saved on
C...a file, and shown with curves suitably superimposed.
      FAC=1D0/NEV
      DO 400 IRE=1,NRE      
      CALL PYOPER(IRE+10,'/',IRE,IRE+10,1D0,1D0)
      CALL PYFACT(IRE,FAC)
      CALL PYFACT(IRE+20,5D0*FAC)
      CALL PYFACT(IRE+30,5D0*FAC)
      CALL PYFACT(IRE+40,5D0*FAC)
      CALL PYFACT(IRE+50,5D0*FAC)
      CALL PYFACT(IRE+60,5D0*FAC)
      CALL PYFACT(IRE+70,FAC)
  400 CONTINUE
      CALL PYHIST

      END 
