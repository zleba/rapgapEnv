C...Supersymmetry at a hadron collider.

C-----------------------------------------------------------------

C...Preamble: declarations.
 
C...All real arithmetic in double precision.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C...Three Pythia functions return integers, so need declaring.
      INTEGER PYK,PYCHGE,PYCOMP
C...Parameter statement to help give large particle numbers
C...(left- and righthanded SUSY, technicolor, excited fermions,
C...extra dimensions).
      PARAMETER (KSUSY1=1000000,KSUSY2=2000000,KTECHN=3000000,
     &KEXCIT=4000000,KDIMEN=5000000)

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

C-----------------------------------------------------------------

C...First section: initialization.

C...Select generic SUSY generation.
      MSEL=39

C...Set SUSY parameters in SUGRA scenario.
      IMSS(1)=2
      RMSS(1)=200D0
      RMSS(4)=1D0
      RMSS(5)=10D0
      RMSS(8)=800D0
      RMSS(16)=0D0

C...If interested only in cross sections and resonance decays: 
C...switch off initial and final state radiation, 
C...multiple interactions and hadronization.
C      MSTP(61)=0
C      MSTP(71)=0
      MSTP(81)=0
C      MSTP(111)=0 

C...Initialization for the Tevatron or LHC. 
C      CALL PYINIT('CMS','p','pbar',2000D0)
      CALL PYINIT('CMS','p','p',14000D0)

C...List resonance data: decay channels, widths etc.
      CALL PYSTAT(2)

C...Histograms for mass distributions.      
      CALL PYBOOK(1,'top mass',100,0D0,1000D0)
      CALL PYBOOK(2,'squark_L mass',100,0D0,1000D0)
      CALL PYBOOK(3,'squark_R mass',100,0D0,1000D0)
      CALL PYBOOK(4,'sbottom_1 mass',100,0D0,1000D0)
      CALL PYBOOK(5,'sbottom_2 mass',100,0D0,1000D0)
      CALL PYBOOK(6,'stop_1 mass',100,0D0,1000D0)
      CALL PYBOOK(7,'stop_2 mass',100,0D0,1000D0)
      CALL PYBOOK(8,'slepton mass',100,0D0,1000D0)
      CALL PYBOOK(9,'gluino mass',100,0D0,1000D0)
      CALL PYBOOK(10,'chi0_1 mass',100,0D0,1000D0)
      CALL PYBOOK(11,'chi0_2 mass',100,0D0,1000D0)
      CALL PYBOOK(12,'chi0_3 mass',100,0D0,1000D0)
      CALL PYBOOK(13,'chi0_4 mass',100,0D0,1000D0)
      CALL PYBOOK(14,'chi+_1 mass',100,0D0,1000D0)
      CALL PYBOOK(15,'chi+_2 mass',100,0D0,1000D0)

C-----------------------------------------------------------------

C...Second section: event loop.

C...Loop over the number of events.
      DO 200 IEV=1,10000
        IF(MOD(IEV,1000).EQ.0) WRITE(6,*) 
     &  'Now at event number',IEV

C...Event generation.
        CALL PYEVNT

C...List first few events.
          IF(IEV.LE.2) CALL PYLIST(1)

C...Loop though documentation section of event.
        DO 150 I=1,MSTP(126)
          IF(K(I,1).NE.21) GOTO 150 

C...Fill the masses of interesting (s)particles.
          KFA=IABS(K(I,2))
          IH=0
          IF(KFA.EQ.6) IH=1
          IF(KFA.GE.KSUSY1+1.AND.KFA.LE.KSUSY1+4) IH=2
          IF(KFA.GE.KSUSY2+1.AND.KFA.LE.KSUSY2+4) IH=3
          IF(KFA.EQ.KSUSY1+5) IH=4
          IF(KFA.EQ.KSUSY2+5) IH=5
          IF(KFA.EQ.KSUSY1+6) IH=6
          IF(KFA.EQ.KSUSY2+6) IH=7
          IF(KFA.GE.KSUSY1+11.AND.KFA.LE.KSUSY1+16) IH=8
          IF(KFA.GE.KSUSY2+11.AND.KFA.LE.KSUSY2+16) IH=8
          IF(KFA.EQ.KSUSY1+21) IH=9
          IF(KFA.EQ.KSUSY1+22) IH=10
          IF(KFA.EQ.KSUSY1+23) IH=11
          IF(KFA.EQ.KSUSY1+25) IH=12
          IF(KFA.EQ.KSUSY1+35) IH=13
          IF(KFA.EQ.KSUSY1+24) IH=14
          IF(KFA.EQ.KSUSY1+37) IH=15
          IF(IH.NE.0) CALL PYFILL(IH,P(I,5),1D0)

C...End of documentation and event loops.
  150   CONTINUE  
  200 CONTINUE

C-----------------------------------------------------------------

C...Third section: produce output and end.

C...Cross section table.
      CALL PYSTAT(1)

C...Histograms.
      CALL PYHIST

      END
