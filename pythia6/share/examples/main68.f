C...Sample main program for a heavy Higgs mass spectrum.

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
 
C...Number of events and cm energy. 
      NEV=1000
      ECM=14000D0

C...Processes.
      MSEL=0
      MSUB(102)=1
      MSUB(123)=1
      MSUB(124)=1

C...Higgs mass.
      PMAS(25,1)=600D0

C...Switch off unnecessary aspects: initial and  final state
C...showers, multiple interactions, hadronization.
C...(Optional for faster simulation of the parton-level 
C...processes only.)
      MSTP(61)=0
      MSTP(71)=0
      MSTP(81)=0
      MSTP(111)=0 

C...Initialize.
      CALL PYINIT('CMS','p','p',ECM)

C...List table of resonance decay channels.
C      CALL PYSTAT(2)
 
C...Book histograms.
      CALL PYBOOK(1,'Higgs mass distribution, old',
     &100,0D0,1000D0)
      CALL PYBOOK(2,'Higgs mass distribution, new',
     &100,0D0,1000D0)

C-----------------------------------------------------------------

C...Second section: event loop.

C...Outer loop over Higgs width treatment.
      DO 250 ICA=1,2
      MSTP(49)=ICA-1 

C...Generate events and look at first few.
      DO 200 IEV=1,NEV
        CALL PYEVNT
        IF(IEV.LE.1) CALL PYLIST(1) 

C...Find and fill Higgs mass.
        DO 150 I=7,9
          IF(K(I,2).EQ.25) CALL PYFILL(ICA,P(I,5),1D0) 
  150   CONTINUE

C...End of loops over events and cases.
  200 CONTINUE
  250 CONTINUE

C-----------------------------------------------------------------

C...Third section: produce output and end.

C...Cross section table.
      CALL PYSTAT(1)

C...Histograms.
      CALL PYHIST

      END
