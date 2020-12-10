C...Run to check that all processes allowed work and to clock them.
C...PYTIME must be interfaced to local clock for this to work.
C...Replace suitably if you want CPU time rather than elapsed one.

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
C...Process information.
      COMMON/PYINT2/ISET(500),KFPR(500,2),COEF(500,20),ICOL(40,4,2)
C...Process names.
      COMMON/PYINT6/PROC(0:500)
      CHARACTER PROC*28
C...Supersymmetry parameters.
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)

C...Local arrays.
      DIMENSION INEP(10),INEE(20),INGAP(10),INGAGA(10),INGAST(10),
     &INEESS(10),INGSTP(10),MSTP14(20),ITIME(6),ITUSED(500),SIGTOT(500)
C...Processes needing other collider types than pp.
      DATA INEP/145,9*0/
      DATA INEE/70,103,146,153,158,343,344,345,346,347,348,9*0/
      DATA INGAP/33,34,35,36,54,80,84,107,2*0/
      DATA INGAGA/58,69,85,108,6*0/
      DATA INGAST/99,137,138,139,140,5*0/
      DATA INEESS/341,342,8*0/
      DATA INGSTP/131,132,133,134,135,136,4*0/
C...Processes needing MSTP(14)=0.
      DATA MSTP14/33,34,35,36,54,58,69,70,80,84,85,103,107,108,153,
     &158,4*0/

C...Default values time consumption.
      DATA ITUSED/500*-1/ 

C-----------------------------------------------------------------

C...First section: common parameters.

C...Number of events per case.
      NEV=100

C...Minimize printout.
      MSTP(122)=0

C...Do not crash if vanishing cross section.
      MSTP(127)=1

C...Set SUSY parameters in SUGRA scenario.
      IMSS(1)=2
      RMSS(1)=200D0
      RMSS(4)=1D0
      RMSS(5)=10D0
      RMSS(8)=800D0
      RMSS(16)=0D0

C-----------------------------------------------------------------

C...Second section: process and event loops.
      CALL PYLIST(0)
 
C...Loop over processes; check which to skip. 
      DO 300 ISUB=1,500
      IF(ISET(ISUB).LT.0.OR.ISET(ISUB).EQ.9) GOTO 300
      WRITE(6,*) 'Now begin process ',ISUB   

C...Switch on desired subprocess and off the rest.
      MSEL=0
      DO 110 J=1,500
 110  MSUB(J)=0
      MSUB(ISUB)=1    

C...Decide on collider type; default is pp.
      ITYPE=1
      DO 120 J=1,10
        IF(ISUB.EQ.INEP(J)) ITYPE=2
        IF(ISUB.EQ.INEE(J)) ITYPE=3
        IF(ISUB.EQ.INEE(J+10)) ITYPE=3
        IF(ISUB.EQ.INGAP(J)) ITYPE=4
        IF(ISUB.EQ.INGAGA(J)) ITYPE=5
        IF(ISUB.EQ.INGAST(J)) ITYPE=6
        IF(ISUB.EQ.INEESS(J)) ITYPE=7
        IF(ISUB.EQ.INGSTP(J)) ITYPE=8
  120 CONTINUE

C...Set MSTP(14)=0 rather than default for some processes..
      MSTP(14)=30
      DO 130 J=1,20
        IF(ISUB.EQ.MSTP14(J)) MSTP(14)=0
  130 CONTINUE
      IF(ISUB.EQ.99) MSTP(14)=26

C...Reset CKIN(3).
      CKIN(3)=0D0

C...Starting time.
      CALL PYTIME(ITIME)
      TOLD=3600D0*ITIME(4)+60D0*ITIME(5)+ITIME(6)
 
C...Initialize for collider type.
C...(These are just convenient choices; could be changed.)
      IF(ITYPE.EQ.1) THEN
        CALL PYINIT('CMS','p','p',14000D0)
      ELSEIF(ITYPE.EQ.2) THEN
        CALL PYINIT('CMS','e+','p',2000D0)
      ELSEIF(ITYPE.EQ.3) THEN
        CALL PYINIT('CMS','e+','e-',2000D0)
      ELSEIF(ITYPE.EQ.4) THEN
        CALL PYINIT('CMS','gamma','p',2000D0)
      ELSEIF(ITYPE.EQ.5) THEN
        CALL PYINIT('CMS','gamma','gamma',2000D0)
      ELSEIF(ITYPE.EQ.6) THEN
        CALL PYINIT('CMS','gamma/e+','gamma/e-',2000D0)
      ELSEIF(ITYPE.EQ.7) THEN
        CALL PYINIT('CMS','e-','e-',200D0)
      ELSEIF(ITYPE.EQ.8) THEN
        CALL PYINIT('CMS','gamma/e-','p',2000D0)
      ENDIF  
 
C...Generate events, unless no cross section found...
      IF(MSTI(53).NE.1) THEN      
        DO 200 IEV=1,NEV
          CALL PYEVNT
  200   CONTINUE
      ELSE
        PARI(1)=-1D0
      ENDIF

C...Finishing time.
      CALL PYTIME(ITIME)
      TNEW=3600D0*ITIME(4)+60D0*ITIME(5)+ITIME(6)

C...Save time and cross sections.
      ITUSED(ISUB)=TNEW-TOLD
      SIGTOT(ISUB)=PARI(1)
      WRITE(6,*) 'Time and cross section: ',ITUSED(ISUB),PARI(1)   

C...End loop over processes.
 300  CONTINUE

C-----------------------------------------------------------------

C...Third section: produce output and end.

C...Loop over processes; write time and cross section.
      DO 400 ISUB=1,500
        IF(ITUSED(ISUB).NE.-1)  WRITE(6,'(2X,I4,3X,A28,I8,1P,E12.4)')
     &  ISUB,PROC(ISUB),ITUSED(ISUB),SIGTOT(ISUB)
  400 CONTINUE   

      END
