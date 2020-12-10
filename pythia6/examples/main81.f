C A total cross section of a single W production at LHC 
C with a PDF set used via LHAPDF library

      PROGRAM MAIN81

c-------------------------------------------------------------
c
c Before running this test the environment variable 
c LHAPATH=/path/to/PDFsets must be set 
c
c (See http://hepforge.cedar.ac.uk/lhapdf/manual for details)
c
c-------------------------------------------------------------

      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
c...Three Pythia functions return integers, so need declaring.
      INTEGER PYK,PYCHGE,PYCOMP

c...EXTERNAL statement links PYDATA on most machines.
      EXTERNAL PYDATA

c...Commonblocks.
c...The event record.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
c...Parameters.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
c...Particle properties + some flavour parameters.
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
c...Decay information.
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
c...Selection of hard scattering subprocesses.
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
c...Parameters. 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
c...Supersymmetry parameters.
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)
C...Generation and cross section statistics.
      COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
C...Random number generation stuff
      COMMON/PYDATR/MRPY(6),RRPY(100)
*
      integer IPDFSET(22)
      logical iopened
c... using only central value PDF sets (*.LHgrid only, as
c... calculating the PDFs "in flight" is too CPU consuming)
c... see http://hepforge.cedar.ac.uk/lhapdf/manual5#tth_sEcA
      data IPDFSET/
     &19050,19051,19060,19070,19150,19160,19170,
     &20050,20052,20053,20054,20060,20070,20250,20270,20350,20370,
     &20450,20470,
     &70050,70150,70250        ! 70350 does not work in 5.2.1, removed!
     &/
      data iopened/.FALSE./
c-----------------------------------------------------------------

c...Main parameters of run: c.m. energy and number of events.
      ECM=14.d3
      NEV=1000

c...Selecting single W+- production:
      MSEL=0
      MSUB(2)=1

c...using LHApdf instead of built-in PDFs:
      MSTP(52) = 2

c... PDF set:
      ITEST=1
c... selecting the PDF set
c... (see http://hepforge.cedar.ac.uk/lhapdf/manual, Appendix A):
      MSTP(51) = IPDFSET(ITEST)
c...Initialize for LHC 
      CALL PYINIT('CMS','p','p',ECM)
c...Event loop:
      DO 02 IEV=1,NEV
02      CALL PYEVNT ()
c...Printing out the cross section table:
      call PYSTAT(1)
c...Printing the total cross section and its uncerainty to testi.dat:
      write(*,100) ITEST, XSEC(0,3),XSEC(0,3)/sqrt(dble(NGEN(0,3)))

100   format('iset, xsec, staterrorxsec ',I3,2E15.6)
      END
