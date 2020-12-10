C...Sample program for jet or minumum bias events 
C...in ep or e+e- by gamma*p or gamma*gamma* interactions,
C...internally convoluted with gamma* flux in e.

C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER PYK,PYCHGE,PYCOMP

C...EXTERNAL statement links PYDATA on most machines.
      EXTERNAL PYDATA

C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      
C...Select machine: 1 = HERA, 2 = LEP2.
      ISEL=1

C...Number of events.
      NEV=1000

C...Three (main) options of event classes:
C...= 1 : jet events only (set pTmin in CKIN(3)).
C...= 2 : jet plus low-pT events (pTmin irrelevant).
C...= 3 : jet plus low-pT plus elastic and diffractive events.  
      MOPT=3
      IF(MOPT.EQ.1) THEN
        CKIN(3)=4D0
      ELSEIF(MOPT.EQ.2) THEN   
        CKIN(3)=0D0
      ELSE
        CKIN(3)=0D0
        MSEL=2    
      ENDIF    

C...Set minimal and maximal W of gamma*gamma* system.
      CKIN(77)=50D0
      CKIN(78)=200D0

C...Set minimum photon virtuality.
c      CKIN(65)=0D0
c      CKIN(67)=0D0    
C...Set maximum photon virtuality.
      CKIN(66)=4D0
      CKIN(68)=4D0  

C...Modify contribution of resolved longitudinal photons.
c      MSTP(17)=4
c      PARP(165)=1D0 

C...Simplify simulation.
c      MSTP(61)=0 
c      MSTP(71)=0 
c      MSTP(81)=0 
c      MSTP(91)=0 
c      MSTP(111)=0 

C...Initialize for HERA.
      IF(ISEL.EQ.1) THEN
        DO 100 J=1,3
          P(1,J)=0D0
          P(2,J)=0D0
 100    CONTINUE
        P(1,3)=-27.5D0
        P(2,3)=820D0
C...Use minimal W to set minimum photon energy fractions.
        ECM=SQRT(4D0*ABS(P(1,3)*P(2,3))) 
        CKIN(61)=(CKIN(77)/ECM)**2  
        CKIN(73)=CKIN(61)
C...Check that not too low pTmin.
        IF(MOPT.EQ.1) THEN
          PTMINS=PARP(81)*(ECM/PARP(89))**PARP(90) 
          CKIN(3)=MAX(CKIN(3),1.01D0*PTMINS)
        ENDIF
        CALL PYINIT('3MOM','GAMMA/E-','P',0D0)

C...Initialize for LEP2.
      ELSE
        ECM=200D0
C...Use minimal W to set minimum photon energy fractions. 
        CKIN(61)=(CKIN(77)/ECM)**2  
        CKIN(63)=CKIN(61)
        CKIN(73)=CKIN(61)
        CKIN(75)=CKIN(61)
C...Check that not too low pTmin.
        IF(MOPT.EQ.1) THEN
          PTMINS=PARP(81)*(ECM/PARP(89))**PARP(90) 
          CKIN(3)=MAX(CKIN(3),1.01D0*PTMINS)
        ENDIF
        CALL PYINIT('CMS','GAMMA/E+','GAMMA/E-',ECM)
      ENDIF 

C...Book histograms.
      CALL PYBOOK(1,'x = energy fraction of photons in lepton',
     &100,0D0,1D0)
      CALL PYBOOK(2,'W distribution',100,0D0,200D0)
      CALL PYBOOK(3,'log10(Q2)',100,-15D0,10D0)
      CALL PYBOOK(4,'x = energy fraction of partons in photon',
     &100,0D0,1D0)
      CALL PYBOOK(5,'pT of hard scattering',100,0D0,20D0)
      CALL PYBOOK(6,'number of hard interactions in event',
     &20,-0.5D0,19.5D0)

C...Event loop. List first few events.
      DO 200 IEV=1,NEV
        CALL PYEVNT
        IF(IEV.LE.2) CALL PYLIST(1)

C...Analyze event.
        ISUB=MSTI(1)
        X1G=PARI(103)
        X2G=PARI(104)
        W=PARI(11)
        Q21=PARI(105)
        Q22=PARI(106)
C...Shift a direct photon away from x_p = 1 for better histogramming.
        X1P=MIN(0.999D0,PARI(33))
        X2P=MIN(0.999D0,PARI(34))
C...pT=0 for low-pT, but nonvanishing for elastic/diffractive events.
        PT=PARI(17)
        NMUL=MSTI(31)

C...Fill histograms. End event loop.
        CALL PYFILL(1,X1G,1D0)
        IF(ISEL.EQ.2) CALL PYFILL(1,X2G,1D0)
        CALL PYFILL(2,W,1D0)
        CALL PYFILL(3,LOG10(Q21),1D0)
        IF(ISEL.EQ.2) CALL PYFILL(3,LOG10(Q22),1D0)
        CALL PYFILL(4,X1P,1D0)
        IF(ISEL.EQ.2) CALL PYFILL(4,X2P,1D0)
        CALL PYFILL(5,PT,1D0)
        CALL PYFILL(6,DBLE(NMUL),1D0)
  200 CONTINUE

C...Final statistics and histograms.
      CALL PYSTAT(1)
      CALL PYHIST 

      END
