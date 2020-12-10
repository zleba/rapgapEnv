C...Sample main program for leptoquark generation at HERA.

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
C...Process information.
      COMMON/PYINT2/ISET(500),KFPR(500,2),COEF(500,20),ICOL(40,4,2)
C...Supersymmetry parameters.
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)

C...Local arrays.
      DIMENSION NCHN(12),QVEC(4)
      DATA NCHN/12*0/ 

C-----------------------------------------------------------------

C...First section: initialization.

C...Number of events, positron and proton energy.
      NEV=1000
      EPOSI=27.5D0
      EPROT=820D0

C...Process.
      MSEL=0
      MSUB(145)=1

C...Leptoquark mass and type: quark and lepton species.
      PMAS(39,1)=200D0
      KFDP(MDCY(39,2),1)=2
      KFDP(MDCY(39,2),2)=-11

C...Rescale leptoquark Yukawa coupling to suitable value.
      PARU(151)=0.01D0

C...Switch off unnecessary aspects: initial and  final state
C...showers, multiple interactions, hadronization.
C...(Optional for faster simulation of the parton-level 
C...processes only.)
C      MSTP(61)=0
C      MSTP(71)=0
C      MSTP(81)=0
C      MSTP(111)=0 

C...Switch off interference between initial and final partons.
C...(Not there because of long leptoquark lifetime.)
      MSTJ(50)=0
      MSTP(67)=0 

C...Initialize.
      P(1,1)=0D0
      P(1,2)=0D0
      P(1,3)=-EPOSI
      P(2,1)=0D0
      P(2,2)=0D0
      P(2,3)=EPROT
      CALL PYINIT('3mom','e+','p',0D0)

C...List table of resonance decay channels.
      CALL PYSTAT(2)
 
C...Book histograms.
      CALL PYBOOK(1,'Leptoquark mass',100,PMAS(39,1)-10D0,
     &PMAS(39,1)+10D0)
      CALL PYBOOK(2,'x spectrum',100,0D0,1D0) 
      CALL PYBOOK(3,'Q2 distribution',100,0D0,1D5) 

C-----------------------------------------------------------------

C...Second section: event loop.
 
C...Begin event loop.
        DO 200 IEV=1,NEV
          CALL PYEVNT

C...List first few events.
          IF(IEV.LE.5) CALL PYLIST(1)

C...Hard scattering mass from commonblock.
          CALL PYFILL(1,PARI(13),1D0) 

C...Reconstruct Q2 and Bjorken x.
          IIN=1
          IP=2
          IOUT=9
  110     IOUT=IOUT+1
          IF(K(IOUT,2).NE.KFDP(MDCY(39,2),2)) GOTO 110
          DO 120 J=1,4
            QVEC(J)=P(IIN,J)-P(IOUT,J)
  120     CONTINUE
          Q2=QVEC(1)**2+QVEC(2)**2+QVEC(3)**2-QVEC(4)**2
          PQ=P(IP,4)*QVEC(4)-P(IP,1)*QVEC(1)-P(IP,2)*QVEC(2)-
     &    P(IP,3)*QVEC(3)
          X=Q2/(2D0*PQ)
          CALL PYFILL(2,X,1D0)
          CALL PYFILL(3,Q2,1D0)

C...End event loop.
  200   CONTINUE

C-----------------------------------------------------------------

C...Third section: produce output and end.

C...Cross section table.
      CALL PYSTAT(1)

C...Histograms.
      CALL PYHIST

      END
