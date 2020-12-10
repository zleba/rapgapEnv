C...Production of a single top by s-channel W exchange.
C...The process is e.g. u dbar -> W+ -> t bbar, with t -> b W+.
C...The problem here is that two W's appear in the process,
C...one as the mother of the t and the other as decay product.
C...Therefore the normal forcing of the W decay does not work;
C...it leads to an infinite recursion W -> t -> W -> t -> ....
C...The way out is to use the W', which has default couplings
C...just like the normal W, just a different mass, which then
C...can be changed to agree. The W' is forced to decay to t bbar,
C...while the W can decay freely (or also be forced, if desired).
C...Additionally, it is important to raise CKIN(1), so that the 
C...program does not get stuck in a region of phase space where
C...the cross section vanishes.  

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

C-----------------------------------------------------------------

C...First section: initialization.

C...Main parameters of run: c.m. energy and number of events.
      ECM=2000D0
      NEV=1000

C...Select W' production.
      MSEL=0
      MSUB(142)=1

C...Set W' mass identical to W (widths already are!).
      PMAS(34,1)=PMAS(24,1)

C...Set W'+ decay only to bbar t.
      DO 100 IDC=MDCY(34,2),MDCY(34,2)+MDCY(34,3)-1
        MDME(IDC,1)=MIN(0,MDME(IDC,1))
        IF(KFDP(IDC,1).EQ.-5.AND.KFDP(IDC,2).EQ.6) MDME(IDC,1)=1
  100 CONTINUE

C...Minimum mass cut.
      CKIN(1)=170D0

C...Initialize.
      CALL PYINIT('CMS','p','pbar',ECM)

C...List resonance initialization parameters.
      CALL PYSTAT(2)

C...Histograms.
      CALL PYBOOK(1,'W* mass distribution',100,0D0,1000D0)

C-----------------------------------------------------------------

C...Second section: event loop.
 
C...Begin event loop.
      DO 200 IEV=1,NEV
        CALL PYEVNT

C...List first few events.
        IF(IEV.LE.2) CALL PYLIST(1)

C...Fill mass.
        CALL PYFILL(1,P(7,5),1D0)

C...End event loop.
  200 CONTINUE

C-----------------------------------------------------------------

C...Third section: produce output and end.

C...Cross section table and histograms.
      CALL PYSTAT(1)
      CALL PYHIST

      END
