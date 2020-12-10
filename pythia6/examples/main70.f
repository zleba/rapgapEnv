C...Long example how to interface user-defined processes to PYTHIA
C...based on the Les Houches commonblock agreement. Generates events 
C...of several different kinds, to test the new code under varied 
C...conditions, but kinematics selection and cross sections are
C...completely unphysical.

C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER PYK,PYCHGE,PYCOMP

C...EXTERNAL statement links PYDATA on most machines.
      EXTERNAL PYDATA

C...User process event common block.
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500) 
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP  
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,IDUP(MAXNUP),
     &ISTUP(MAXNUP),MOTHUP(2,MAXNUP),ICOLUP(2,MAXNUP),PUP(5,MAXNUP),
     &VTIMUP(MAXNUP),SPINUP(MAXNUP)  
 
C...Standard PYTHIA commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)

C...Extra commonblock to transfer run info.
      COMMON/PRIV/MODE,NLIM  

C...Local arrays.
      DIMENSION IPRID(100)

C...Switch process mode; agrees with IDWTUP code (+-1,+-2,+-3,+-4).
      MODE=-2

C...Simulate limited suppy of events of Wbb kind.
      NLIM=10000

C...Maximum number of events to generate.
      NEV=1000

C...Set pi0 stable to trim event listings.
      MDCY(PYCOMP(111),1)=0

C...Expanded event listing (required for histogramming).
      MSTP(125)=2

C...Histograms.
      CALL PYBOOK(1,'charged multiplicity',100,-1D0,399D0)
      CALL PYBOOK(2,'starting virtualities, 2 -> 2',100,0D0,200D0)
      CALL PYBOOK(3,'starting virtualities, 2 -> 3',100,0D0,200D0)
      CALL PYBOOK(4,'starting virtualities, 2 -> 4',100,0D0,200D0)
      CALL PYBOOK(5,'starting virtualities, 2 -> 5',100,0D0,200D0)
      CALL PYBOOK(6,'starting virtualities, 2 -> 6',100,0D0,200D0)
      CALL PYBOOK(7,'starting virtualities, 2 -> 7',100,0D0,200D0)

C...Initialize with external process.
      CALL PYINIT('USER',' ',' ',0D0)
      NACC=0
      NPRID=0

C...Event loop; generate event; check it was obtained or quit.
      DO 130 IEV=1,NEV
        CALL PYEVNT
        IF(MSTI(51).EQ.1) GOTO 140  
        NACC=NACC+1

C...List one event of each new type.
        ISUB=MSTI(1)
        IAGR=0
        DO 100 I=1,NPRID
          IF(ISUB.EQ.IPRID(I)) IAGR=I
  100   CONTINUE
        IF(IAGR.EQ.0) THEN
          NPRID=NPRID+1
          IPRID(NPRID)=ISUB
          CALL PYLIST(7)
          CALL PYLIST(2)
        ENDIF
 
C...Analyze events. End of event loop.
        IF(IDPRUP.GE.212.AND.IDPRUP.LE.217) THEN
          NGLU=MOD(IDPRUP,10)
          ICM=0
          DO 110 I=1,N
            IF(K(I,2).EQ.94) ICM=I
  110     CONTINUE
          IF(ICM.NE.0) THEN
            DO 120 I=1,NGLU
              CALL PYFILL(NGLU,P(ICM+I,5),XWGTUP)
  120       CONTINUE
          ENDIF
        ENDIF  
        CALL PYEDIT(3)
        CALL PYFILL(1,DBLE(N),1D0)
  130 CONTINUE

C...Statistics and histograms.
  140 CALL PYSTAT(1)
      CALL PYHIST

      END
   
C*********************************************************************
 
C...UPINIT
C...Routine to be called by user to set up user-defined processes.
C...Code below only intended as example, without any claim of realism.
C...Especially it shows what info needs to be put in HEPRUP.
 
      SUBROUTINE UPINIT
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)

C...User process initialization commonblock.
      INTEGER MAXPUP
      PARAMETER (MAXPUP=100)
      INTEGER IDBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,LPRUP
      DOUBLE PRECISION EBMUP,XSECUP,XERRUP,XMAXUP
      COMMON/HEPRUP/IDBMUP(2),EBMUP(2),PDFGUP(2),PDFSUP(2),
     &IDWTUP,NPRUP,XSECUP(MAXPUP),XERRUP(MAXPUP),XMAXUP(MAXPUP),
     &LPRUP(MAXPUP)
      SAVE /HEPRUP/

C...Extra commonblock to transfer run info.
      COMMON/PRIV/MODE,NLIM   
      SAVE/PRIV/

C....Pythia commonblock - needed for setting PDF's; see below.
C      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
C      SAVE /PYPARS/      

C...Set incoming beams: Tevatron Run II.
      IDBMUP(1)=2212
      IDBMUP(2)=-2212
      EBMUP(1)=1000D0
      EBMUP(2)=1000D0

C...Set PDF's of incoming beams: CTEQ 5L.
C...Note that Pythia will not look at PDFGUP and PDFSUP.  
      PDFGUP(1)=4
      PDFSUP(1)=46
      PDFGUP(2)=PDFGUP(1)
      PDFSUP(2)=PDFSUP(1)
      
C...If you want Pythia to use PDFLIB, you have to set it by hand.
C...(You also have to ensure that the dummy routines
C...PDFSET, STRUCTM and STRUCTP in Pythia are not linked.)      
C      MSTP(52)=2
C      MSTP(51)=1000*PDFGUP(1)+PDFSUP(1)

C...Decide on weighting strategy: unweighted on input.
      IDWTUP=MODE

C...Number of external processes. 
      NPRUP=9

C...Set up q qbar -> t tbar.
      XSECUP(1)=0.6D0
      XMAXUP(1)=0.8D0
      LPRUP(1)=661

C...Set up g g -> t tbar.
      XSECUP(2)=0.06D0
      XMAXUP(2)=0.15D0
      LPRUP(2)=662

C...Set up u dbar -> W+ b bbar.
      XSECUP(3)=0.5D0
      XMAXUP(3)=0.5D0
      LPRUP(3)=2455

C...Set up g g -> 2 to 7 gluons.
      DO 100 IPR=4,9
        XSECUP(IPR)=0.1D0
        XMAXUP(IPR)=1D0
        LPRUP(IPR)=208+IPR
  100 CONTINUE

C...Negative weights for some processes.
      MODEA=IABS(MODE)
      DO 110 I=1,9
        IF(MODE.LT.0.AND.MOD(I,2).EQ.0) THEN
          XSECUP(I)=-XSECUP(I)
          XMAXUP(I)=-XMAXUP(I)
        ENDIF

C...In some modes XSECUP or XMAXUP need not be given.
        IF(MODEA.EQ.1.OR.MODEA.EQ.4) XSECUP(I)=0D0
        IF(MODEA.EQ.3.OR.MODEA.EQ.4) XMAXUP(I)=0D0
  110 CONTINUE
        

      RETURN
      END
 
C*********************************************************************
 
C...UPEVNT
C...Sample routine to generate events of various kinds.
C...Not intended to be realistic, but rather to show in closed
C...and understandable form what such a routine might look like.
C...Especially it shows what info needs to be put in HEPEUP.
 
      SUBROUTINE UPEVNT
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)

C...Extra commonblock to transfer run info.
      COMMON/PRIV/MODE,NLIM   
      SAVE/PRIV/

C...User process event common block.
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500) 
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP  
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,IDUP(MAXNUP),
     &ISTUP(MAXNUP),MOTHUP(2,MAXNUP),ICOLUP(2,MAXNUP),PUP(5,MAXNUP),
     &VTIMUP(MAXNUP),SPINUP(MAXNUP)  
      SAVE /HEPEUP/ 

C...If PYTHIA is supposed to select event type, do not modify this choice.
      IF(IABS(MODE).LE.2) THEN

C...Else free hands to mix; here evenly.
      ELSE  
        IPYR=1+MIN(8,INT(9D0*PYR(0)))
        IF(IPYR.LE.2) THEN
          IDPRUP=660+IPYR
        ELSEIF(IPYR.EQ.3) THEN
          IDPRUP=2455
        ELSE
          IDPRUP=208+IPYR
        ENDIF 
      ENDIF

C...Call the respective routine to generate event.
      IF(IDPRUP.EQ.661.OR.IDPRUP.EQ.662) THEN
        CALL MYTTB
      ELSEIF(IDPRUP.EQ.2455) THEN
        CALL MYWBB
      ELSEIF(IDPRUP.GE.212.AND.IDPRUP.LE.217) THEN
        CALL MYGLU
      ELSE
        WRITE(*,*) 'Fatal error! Unknown process',IDPRUP 
        STOP  
      ENDIF

C...Ensure proper normalization of weights for MODE=+-3
      IF(MODE.EQ.3) THEN
        XWGTUP=1D0
      ELSEIF(MODE.EQ.-3) THEN
        XWGTUP=SIGN(1D0,XWGTUP)
      ENDIF   

      RETURN
      END 
 
C*********************************************************************
 
C...MYTTB
C...Sample routine to generate q qbar or g + g -> t tbar events.
C...Not intended to be realistic
 
      SUBROUTINE MYTTB
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      PARAMETER (PI=3.141592653589793D0)

C...Extra commonblock to transfer run info.
      COMMON/PRIV/MODE,NLIM   
      SAVE/PRIV/

C...User process event common block.
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500) 
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP  
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,IDUP(MAXNUP),
     &ISTUP(MAXNUP),MOTHUP(2,MAXNUP),ICOLUP(2,MAXNUP),PUP(5,MAXNUP),
     &VTIMUP(MAXNUP),SPINUP(MAXNUP)  
      SAVE /HEPEUP/ 

C...PYTHIA commonblock.  
      COMMON/PYINT1/MINT(400),VINT(400)

C...CM energy of system.
      ECM=VINT(1)

C...Input (to be provided from common block eventually).
C...Top mass and Breit-Wigner parameters.
      PMTOP=175D0
      PWTOP=1.5D0
      PMTOPL=150D0
      PMTOPU=200D0
C....W mass and Breit-Wigner parameters. 
      PMW=80.4D0
      PWW=2.1D0
      PMWL=60D0
      PMWU=100D0
C....b mass.
      PMB=4.5D0 

C...Zero some arrays in common blocks to simplify filling.
      NUP=12
      DO 100 I=1,NUP
        MOTHUP(1,I)=0
        MOTHUP(2,I)=0
        ICOLUP(1,I)=0
        ICOLUP(2,I)=0
        SPINUP(I)=9D0
        PUP(1,I)=0D0    
        PUP(2,I)=0D0    
        PUP(5,I)=0D0
        VTIMUP(I)=1D0
  100 CONTINUE

C...Generate top and antitop masses according to Breit-Wigners.
      ATL=ATAN((PMTOPL**2-PMTOP**2)/(PMTOP*PWTOP))
      ATU=ATAN((PMTOPU**2-PMTOP**2)/(PMTOP*PWTOP))
      PMBW=PMTOP**2+PMTOP*PWTOP*TAN(ATL+PYR(0)*(ATU-ATL))
      PMT1=MIN(PMTOPU,MAX(PMTOPL,SQRT(MAX(0D0,PMBW))))
      PMBW=PMTOP**2+PMTOP*PWTOP*TAN(ATL+PYR(0)*(ATU-ATL))
      PMT2=MIN(PMTOPU,MAX(PMTOPL,SQRT(MAX(0D0,PMBW))))
          
C...Generate flavour, tau and y of q qbar or g g  pair. 
  110 CONTINUE
      IFL=MAX(1,MIN(3,INT(1D0+2.2D0*PYR(0))))
      IF(PYR(0).LT.0.2D0) IFL=-IFL
      IF(IDPRUP.EQ.662) IFL=21 
      TAUMIN=(PMT1+PMT2)**2/ECM**2
      TAU=TAUMIN**PYR(0)  
      YMAX=-0.5D0*LOG(TAU)
      Y=YMAX*(2D0*PYR(0)-1D0)
      X1=SQRT(TAU)*EXP(Y)
      X2=SQRT(TAU)*EXP(-Y)
      
C...Generate scattering angles; translate to Mandelstam variables.
      COSTHE=2D0*PYR(0)-1D0
      SHAT=TAU*ECM**2
      BETA34=SQRT(MAX(0D0,(1D0-PMT1**2/SHAT-PMT2**2/SHAT)**2-
     &4D0*(PMT1**2/SHAT)*(PMT2**2/SHAT)))
      THAT=-0.5D0*(SHAT-PMT1**2-PMT2**2-SHAT*BETA34*COSTHE)
      UHAT=-0.5D0*(SHAT-PMT1**2-PMT2**2+SHAT*BETA34*COSTHE)
      PT2HAT=0.25D0*SHAT*BETA34**2*(1D0-COSTHE**2)

C...Here is excellent place to evaluate parton densities and
C...matrix elements. This is currently not done, except primitively.
C...If you want unweighted events, goto 100 in case of failure,
C...else carry along weight as part of event weight.
      IF(IFL.NE.21) THEN   
        WT=(1D0-X1)**3*(1D0-X2)**6
      ELSE
        WT=(1D0-X1)**4*(1D0-X2)**7
      ENDIF
      IF(WT.LT.PYR(0)) GOTO 110 

C...Set up flavour and history of q + qbar or g + g -> t + tbar. 
      IDUP(1)=IFL
      IDUP(2)=-IFL
      IF(IFL.EQ.21) IDUP(2)=IFL
      IDUP(3)=6
      IDUP(4)=-6 
      ISTUP(1)=-1
      ISTUP(2)=-1
      ISTUP(3)=2
      ISTUP(4)=2
      MOTHUP(1,3)=1
      MOTHUP(2,3)=2
      MOTHUP(1,4)=1
      MOTHUP(2,4)=2

C...Set up colour of g +g or q + qbar -> t + tbar. 
      IF(IFL.EQ.21) THEN
        ICOLUP(1,1)=101
        ICOLUP(2,2)=102
        ICOLUP(2,1)=109
        ICOLUP(1,2)=109
      ELSEIF(IFL.GT.0) THEN
        ICOLUP(1,1)=101
        ICOLUP(2,2)=102
      ELSE
        ICOLUP(2,1)=102
        ICOLUP(1,2)=101
      ENDIF         
      ICOLUP(1,3)=101
      ICOLUP(2,4)=102

C...Set up kinematics of q + qbar -> t + tbar. 
      PHI=2D0*PI*PYR(0)
      PUP(4,1)=0.5D0*X1*ECM
      PUP(3,1)=PUP(4,1)
      PUP(4,2)=0.5D0*X2*ECM
      PUP(3,2)=-PUP(4,2)
      PUP(1,3)=SQRT(PT2HAT)*COS(PHI)  
      PUP(2,3)=SQRT(PT2HAT)*SIN(PHI)  
      PUP(3,3)=0.5D0*SQRT(SHAT)*BETA34*COSTHE
      PUP(4,3)=0.5D0*(SHAT+PMT1**2-PMT2**2)/SQRT(SHAT)  
      PUP(5,3)=PMT1
      PUP(1,4)=-PUP(1,3)  
      PUP(2,4)=-PUP(2,3)  
      PUP(3,4)=-PUP(3,3)  
      PUP(4,4)=0.5D0*(SHAT+PMT2**2-PMT1**2)/SQRT(SHAT)  
      PUP(5,4)=PMT2
      CALL BSTUP(3,4,0D0,0D0,(X1-X2)/(X1+X2))

C...Generate W+- masses according to Breit-Wigners.
      AWL=ATAN((PMWL**2-PMW**2)/(PMW*PWW))
      AWU=ATAN((PMWU**2-PMW**2)/(PMW*PWW))
      PMBW=PMW**2+PMW*PWW*TAN(AWL+PYR(0)*(AWU-AWL))
      PMW1=MIN(PMWU,MAX(PMWL,SQRT(MAX(0D0,PMBW))))
      PMBW=PMW**2+PMW*PWW*TAN(AWL+PYR(0)*(AWU-AWL))
      PMW2=MIN(PMWU,MAX(PMWL,SQRT(MAX(0D0,PMBW))))

C...Set up flavour, history  and colour of two t -> b + W decays.  
      IDUP(5)=5
      IDUP(6)=24
      IDUP(7)=-5
      IDUP(8)=-24
      ISTUP(5)=1
      ISTUP(6)=2
      ISTUP(7)=1
      ISTUP(8)=2
      MOTHUP(1,5)=3
      MOTHUP(1,6)=3
      MOTHUP(1,7)=4
      MOTHUP(1,8)=4
      ICOLUP(1,5)=101
      ICOLUP(2,7)=102

C...Set up flavour, history  and colour of two W -> u dbar decays.
      IDUP(9)=2
      IDUP(10)=-1 
      IDUP(11)=-2
      IDUP(12)=1 
      ISTUP(9)=1
      ISTUP(10)=1
      ISTUP(11)=1
      ISTUP(12)=1
      MOTHUP(1,9)=6
      MOTHUP(1,10)=6
      MOTHUP(1,11)=8
      MOTHUP(1,12)=8
      ICOLUP(1,9)=103
      ICOLUP(2,10)=103
      ICOLUP(2,11)=104
      ICOLUP(1,12)=104

C...Construct top decay kinematics isotropically in angle and boost. 
  120 COSTHE=2D0*PYR(0)-1D0
      PHI=2D0*PI*PYR(0)
      PABS=0.5D0*SQRT((PMT1**2-PMW1**2-PMB**2)**2-4D0*PMW1**2*PMB**2)/
     &PMT1
      PUP(1,5)=PABS*SQRT(1D0-COSTHE**2)*COS(PHI)
      PUP(2,5)=PABS*SQRT(1D0-COSTHE**2)*SIN(PHI)
      PUP(3,5)=PABS*COSTHE
      PUP(4,5)=0.5D0*(PMT1**2+PMB**2-PMW1**2)/PMT1
      PUP(5,5)=PMB
      PUP(1,6)=-PUP(1,5)
      PUP(2,6)=-PUP(2,5)
      PUP(3,6)=-PUP(3,5)
      PUP(4,6)=0.5D0*(PMT1**2+PMW1**2-PMB**2)/PMT1
      PUP(5,6)=PMW1
      CALL BSTUP(5,6,PUP(1,3)/PUP(4,3),PUP(2,3)/PUP(4,3),
     &PUP(3,3)/PUP(4,3)) 
      COSTHE=2D0*PYR(0)-1D0
      PHI=2D0*PI*PYR(0)
      PABS=0.5D0*SQRT((PMT2**2-PMW2**2-PMB**2)**2-4D0*PMW2**2*PMB**2)/
     &PMT2
      PUP(1,7)=PABS*SQRT(1D0-COSTHE**2)*COS(PHI)
      PUP(2,7)=PABS*SQRT(1D0-COSTHE**2)*SIN(PHI)
      PUP(3,7)=PABS*COSTHE
      PUP(4,7)=0.5D0*(PMT2**2+PMB**2-PMW2**2)/PMT2
      PUP(5,7)=PMB
      PUP(1,8)=-PUP(1,7)
      PUP(2,8)=-PUP(2,7)
      PUP(3,8)=-PUP(3,7)
      PUP(4,8)=0.5D0*(PMT2**2+PMW2**2-PMB**2)/PMT2
      PUP(5,8)=PMW2
      CALL BSTUP(7,8,PUP(1,4)/PUP(4,4),PUP(2,4)/PUP(4,4),
     &PUP(3,4)/PUP(4,4)) 

C...Construct W decay kinematics isotropically in angle and boost. 
      COSTHE=2D0*PYR(0)-1D0
      PHI=2D0*PI*PYR(0)
      PUP(1,9)=0.5D0*PMW1*SQRT(1D0-COSTHE**2)*COS(PHI)
      PUP(2,9)=0.5D0*PMW1*SQRT(1D0-COSTHE**2)*SIN(PHI)
      PUP(3,9)=0.5D0*PMW1*COSTHE
      PUP(4,9)=0.5D0*PMW1
      PUP(1,10)=-PUP(1,9)
      PUP(2,10)=-PUP(2,9)
      PUP(3,10)=-PUP(3,9)
      PUP(4,10)=PUP(4,9)
      CALL BSTUP(9,10,PUP(1,6)/PUP(4,6),PUP(2,6)/PUP(4,6),
     &PUP(3,6)/PUP(4,6)) 
      COSTHE=2D0*PYR(0)-1D0
      PHI=2D0*PI*PYR(0)
      PUP(1,11)=0.5D0*PMW2*SQRT(1D0-COSTHE**2)*COS(PHI)
      PUP(2,11)=0.5D0*PMW2*SQRT(1D0-COSTHE**2)*SIN(PHI)
      PUP(3,11)=0.5D0*PMW2*COSTHE
      PUP(4,11)=0.5D0*PMW2
      PUP(1,12)=-PUP(1,11)
      PUP(2,12)=-PUP(2,11)
      PUP(3,12)=-PUP(3,11)
      PUP(4,12)=PUP(4,11)
      CALL BSTUP(11,12,PUP(1,8)/PUP(4,8),PUP(2,8)/PUP(4,8),
     &PUP(3,8)/PUP(4,8)) 

C...Now all decay kinematics is known. Here is excellent place to
C...insert weighting of the decay angle correlations.  
C...If you want unweighted events, goto 120 in case of failure.
      IF(IDPRUP.EQ.661) THEN
        XWGTUP=1D0
      ELSE
        XWGTUP=0.2D0   
        IF(MODE.LT.0) XWGTUP=-0.2D0
      ENDIF  

C...Some other compulsory quantities.
      SCALUP=PMTOP

      RETURN
      END
 
C*********************************************************************
 
C...MYWBB
C...Sample routine to generate u dbar -> W+ b bbar events.
C...Not intended to be realistic
 
      SUBROUTINE MYWBB
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      PARAMETER (PI=3.141592653589793D0)

C...Extra commonblock to transfer run info.
      COMMON/PRIV/MODE,NLIM   
      SAVE/PRIV/

C...User process event common block.
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500) 
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP  
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,IDUP(MAXNUP),
     &ISTUP(MAXNUP),MOTHUP(2,MAXNUP),ICOLUP(2,MAXNUP),PUP(5,MAXNUP),
     &VTIMUP(MAXNUP),SPINUP(MAXNUP)  
      SAVE /HEPEUP/   

C...PYTHIA commonblock.  
      COMMON/PYINT1/MINT(400),VINT(400)

C...Event number counter.
      DATA ILIM/0/

C...CM energy of system.
      ECM=VINT(1)

C...Simulate limited supply of events.
      ILIM=ILIM+1
      IF(ILIM.GT.NLIM) THEN
        WRITE(*,*) 'No more Wbb events!'
        NUP=0  
        RETURN
      ENDIF

C...Input (to be provided from common block eventually).
C....W mass and Breit-Wigner parameters. 
      PMW=80.4D0
      PWW=2.1D0
      PMWL=60D0
      PMWU=100D0
C....b mass.
      PMB=4.5D0 
          
C...Zero some arrays in common blocks to simplify filling.
      NUP=5
      DO 100 I=1,NUP
        MOTHUP(1,I)=0
        MOTHUP(2,I)=0
        ICOLUP(1,I)=0
        ICOLUP(2,I)=0
        SPINUP(I)=9D0
        PUP(1,I)=0D0    
        PUP(2,I)=0D0    
        PUP(5,I)=0D0
        VTIMUP(I)=0D0
  100 CONTINUE

C...Set up flavour and history of u dbar -> W+ b bbar. 
      IDUP(1)=2
      IDUP(2)=-1
      IDUP(3)=24
      IDUP(4)=5 
      IDUP(5)=-5 
      ISTUP(1)=-1
      ISTUP(2)=-1
      ISTUP(3)=1
      ISTUP(4)=1
      ISTUP(5)=1
      MOTHUP(1,3)=1
      MOTHUP(2,3)=2
      MOTHUP(1,4)=1
      MOTHUP(2,4)=2
      MOTHUP(1,5)=1
      MOTHUP(2,5)=2

C...Set up colour of u dbar -> W+ b bbar.  
      ICOLUP(1,1)=101
      ICOLUP(1,4)=101
      ICOLUP(2,2)=102
      ICOLUP(2,5)=102

C...Set up kinematics of b and bbar.
      PTBMX=50D0
      DO 120 I=4,5
        PUP(5,I)=PMB
        RR=PYR(0)
        PTB=PMB*PTBMX*SQRT(RR/(PMB**2+(1D0-RR)*PTBMX**2))
        PHIB=2D0*PI*PYR(0)
        PUP(1,I)=PTB*COS(PHIB) 
        PUP(2,I)=PTB*SIN(PHIB) 
        PUP(3,I)=100D0*(2D0*PYR(0)-1D0)
        PUP(4,I)=SQRT(PUP(1,I)**2+PUP(2,I)**2+PUP(3,I)**2+PUP(5,I)**2)
  120 CONTINUE

C...Generate W+ mass according to Breit-Wigner.
      AWL=ATAN((PMWL**2-PMW**2)/(PMW*PWW))
      AWU=ATAN((PMWU**2-PMW**2)/(PMW*PWW))
      PMBW=PMW**2+PMW*PWW*TAN(AWL+PYR(0)*(AWU-AWL))
      PUP(5,3)=MIN(PMWU,MAX(PMWL,SQRT(MAX(0D0,PMBW))))

C...Set up kinematics of W+.
      PUP(1,3)=-PUP(1,4)-PUP(1,5)   
      PUP(2,3)=-PUP(2,4)-PUP(2,5)   
      PUP(3,3)=100D0*(2D0*PYR(0)-1D0)
      PUP(4,3)=SQRT(PUP(1,3)**2+PUP(2,3)**2+PUP(3,3)**2+PUP(5,3)**2)

C...Set up kinematics of incoming u dbar.
      ESUM=PUP(4,3)+PUP(4,4)+PUP(4,5)
      PZSUM=PUP(3,3)+PUP(3,4)+PUP(3,5)
      PUP(4,1)=0.5D0*(ESUM+PZSUM)
      PUP(4,2)=0.5D0*(ESUM-PZSUM)
      PUP(3,1)=PUP(4,1) 
      PUP(3,2)=-PUP(4,2) 

C...Now all decay kinematics is known. Fictitious weight.
      XWGTUP=PYR(0)

C...Some other compulsory quantities.
      SCALUP=PMW

      RETURN
      END
 
C*********************************************************************
 
C...MYGLU
C...Sample routine to generate g g -> 2 to 7 gluons.
C...Not intended to be realistic.
 
      SUBROUTINE MYGLU
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      PARAMETER (PI=3.141592653589793D0)

C...Extra commonblock to transfer run info.
      COMMON/PRIV/MODE,NLIM   
      SAVE/PRIV/

C...User process event common block.
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500) 
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP  
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,IDUP(MAXNUP),
     &ISTUP(MAXNUP),MOTHUP(2,MAXNUP),ICOLUP(2,MAXNUP),PUP(5,MAXNUP),
     &VTIMUP(MAXNUP),SPINUP(MAXNUP)  
      SAVE /HEPEUP/ 

C...PYTHIA commonblock.  
      COMMON/PYINT1/MINT(400),VINT(400)

C...Local arrays.
      DIMENSION PSUM(4)

C...Number of gluons.
      NGLU=MOD(IDPRUP,10)

C...CM energy of system.
      ECM=VINT(1)
          
C...Zero some arrays in common blocks to simplify filling.
      NUP=NGLU+2
      DO 110 I=1,NUP          
        IDUP(I)=21
        MOTHUP(1,I)=0
        MOTHUP(2,I)=0
        ICOLUP(1,I)=0
        ICOLUP(2,I)=0
        SPINUP(I)=9D0
        PUP(1,I)=0D0    
        PUP(2,I)=0D0    
        PUP(5,I)=0D0
        VTIMUP(I)=0D0
  110 CONTINUE

C...Pick scattering subsystem energy. 
      EHAT=MAX(10D0,0.9D0*ECM*PYR(0)**2)

C...Set up mothers.
      ISTUP(1)=-1
      ISTUP(2)=-1
      PUP(4,1)=0.5D0*EHAT
      PUP(4,2)=0.5D0*EHAT
      PUP(3,1)=PUP(4,1)
      PUP(3,2)=-PUP(4,2)

C...Pick outgoing momenta for daughters, and set mother pointers.
      DO 120 J=1,4
        PSUM(J)=0D0
  120 CONTINUE  
      DO 140 I=3,NUP
        ISTUP(I)=1
        MOTHUP(1,I)=1
        MOTHUP(2,I)=2 
        PABS=-LOG(PYR(0)*PYR(1))
        CTHE=2D0*PYR(0)-1D0
        STHE=SQRT(MAX(0D0,1D0-CTHE**2))
        PHI=2D0*PI*PYR(0)
        PUP(1,I)=PABS*STHE*COS(PHI)
        PUP(2,I)=PABS*STHE*SIN(PHI)
        PUP(3,I)=PABS*CTHE
        PUP(4,I)=PABS
        DO 130 J=1,4
          PSUM(J)=PSUM(J)+PUP(J,I)
  130   CONTINUE  
  140 CONTINUE 

C...Boost to rest frame and rescale.
      CALL BSTUP(3,NUP,-PSUM(1)/PSUM(4),-PSUM(2)/PSUM(4),
     &-PSUM(3)/PSUM(4))
      ESUM=0D0
      DO 150 I=3,NUP
        ESUM=ESUM+PUP(4,I)
  150 CONTINUE
      FAC=EHAT/ESUM
      DO 170 I=3,NUP
        DO 160 J=1,4
          PUP(J,I)=FAC*PUP(J,I)
  160   CONTINUE  
  170 CONTINUE 
   
C...Pick rapidity and do boost.
      YMAX=LOG(ECM/EHAT)
  180 Y=(2D0*PYR(0)-1D0)*YMAX
      IF((1-ABS(Y)/YMAX)**2.LT.PYR(0)) GOTO 180
      BETAZ=(EXP(2D0*Y)-1D0)/(EXP(2D0*Y)+1D0) 
      CALL BSTUP(1,NUP,0D0,0D0,BETAZ)      

C...Now all decay kinematics is known. Here is excellent place to
C...insert weighting of the decay angle correlations.  
      XWGTUP=PYR(0)**NGLU

C...Simulate negative weights.
      IF(MODE.LT.0.AND.MOD(NGLU,2).EQ.0) XWGTUP=-XWGTUP 

C...Set colour flow.
      ICOLUP(2,1)=501
      ICOLUP(1,2)=501
      ICOLUP(1,1)=502
      ICOLUP(1,3)=502
      ICOLUP(2,2)=503
      ICOLUP(2,NUP)=503
      DO 190 I=3,NUP-1
        ICOLUP(2,I)=501+I
        ICOLUP(1,I+1)=501+I
  190 CONTINUE            

C...Some other compulsory quantities.
      SCALUP=EHAT/NGLU

      RETURN
      END
 
C*********************************************************************
 
C...BSTUP
C...Performs boosts on user-process entries.
 
      SUBROUTINE BSTUP(IMIN,IMAX,BETAX,BETAY,BETAZ)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)

C...User process event common block.
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500) 
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP  
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,IDUP(MAXNUP),
     &ISTUP(MAXNUP),MOTHUP(2,MAXNUP),ICOLUP(2,MAXNUP),PUP(5,MAXNUP),
     &VTIMUP(MAXNUP),SPINUP(MAXNUP)  
      SAVE /HEPEUP/   
 
C...Boost, typically from rest to momentum/energy=beta.
      GAMMA=1D0/SQRT(1D0-BETAX**2-BETAY**2-BETAZ**2)
      DO 100 I=IMIN,IMAX
        BETAP=BETAX*PUP(1,I)+BETAY*PUP(2,I)+BETAZ*PUP(3,I)
        GABEP=GAMMA*(GAMMA*BETAP/(1D0+GAMMA)+PUP(4,I))
        PUP(1,I)=PUP(1,I)+GABEP*BETAX
        PUP(2,I)=PUP(2,I)+GABEP*BETAY
        PUP(3,I)=PUP(3,I)+GABEP*BETAZ
        PUP(4,I)=GAMMA*(PUP(4,I)+BETAP)
  100 CONTINUE
 
      RETURN
      END
