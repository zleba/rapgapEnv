C...Sample main program and subroutines illustrating how to set up
C...an external process in Pythia, including cross section calculation.
C...Example: q + g -> q + gamma at the Tevatron. This process 
C...already exists in standard Pythia, som comparisons are possible
C...between the two implementations.
C...Reminder: since dummy copies of the UPINIT and UPEVNT routines 
C...are included in the Pythia standard distribution (to avoid 
C...problems with unresolved external references) it is necessary 
C...to remove (the linking of) this dummy to avoid clashes with the 
C...routines below.

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
C...Parameters. 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)

C...The user's own transfer of information.
      COMMON/MYCOMM/ECM,PTMIN,INIT
      SAVE/MYCOMM/

C...Counter for number of generated events of each type.
      DIMENSION NCOUNT(2)
      DATA NCOUNT/2*0/

C-----------------------------------------------------------------

C...First section: initialization.
 
C...Number of events and cm energy. 
      NEV=1000
      ECM=2000D0

C...Define a pTmin cutoff for the process.
      PTMIN=10D0

C...Switch off unnecessary aspects: initial and  final state
C...showers, multiple interactions, hadronization.
C...(Optional for faster simulation of the parton-level 
C...processes only.)
C      MSTP(61)=0
C      MSTP(71)=0
C      MSTP(81)=0
C      MSTP(111)=0 

C...Initialize.
      CALL PYINIT('USER',' ',' ',0D0)
 
C...Book histograms.
      CALL PYBOOK(1,'pT spectrum internal',100,0D0,100D0)
      CALL PYBOOK(2,'pT spectrum external',100,0D0,100D0)
      CALL PYBOOK(3,'n_ch distribution internal',100,-1D0,199D0)
      CALL PYBOOK(4,'n_ch distribution external',100,-1D0,199D0)

C-----------------------------------------------------------------

C...Second section: event loop.

C...Generate events and look at first few.
      DO 200 IEV=1,NEV
        CALL PYEVNT
        ISUB=MSTI(1)
        ICASE=1
        IF(ISUB.NE.29) ICASE=2
        NCOUNT(ICASE)=NCOUNT(ICASE)+1
        IF(NCOUNT(ICASE).LE.1) THEN
          WRITE(6,*) ' Following event is subprocess',ISUB 
          IF(ICASE.EQ.2) CALL PYLIST(7)
          CALL PYLIST(1)
        ENDIF

C...Read out generated pT distribution and charged multiplicity.
        PTGEN=PARI(17)
        CALL PYFILL(ICASE,PTGEN,1D0)
        CALL PYEDIT(3)
        CALL PYFILL(ICASE+2,DBLE(N),1D0)

C...End of loop over events.
  200 CONTINUE

C-----------------------------------------------------------------

C...Third section: produce output and end.

C...Cross section table.
      CALL PYSTAT(1)

C...Histograms.
      CALL PYHIST

      END
 
C*********************************************************************

      SUBROUTINE UPINIT
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)

C...Selection of hard scattering subprocesses.
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
 
C...User process initialization commonblock.
      INTEGER MAXPUP
      PARAMETER (MAXPUP=100)
      INTEGER IDBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,LPRUP
      DOUBLE PRECISION EBMUP,XSECUP,XERRUP,XMAXUP
      COMMON/HEPRUP/IDBMUP(2),EBMUP(2),PDFGUP(2),PDFSUP(2),
     &IDWTUP,NPRUP,XSECUP(MAXPUP),XERRUP(MAXPUP),XMAXUP(MAXPUP),
     &LPRUP(MAXPUP)
      SAVE /HEPRUP/
 
C...User process event common block.
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,IDUP(MAXNUP),
     &ISTUP(MAXNUP),MOTHUP(2,MAXNUP),ICOLUP(2,MAXNUP),PUP(5,MAXNUP),
     &VTIMUP(MAXNUP),SPINUP(MAXNUP)
      SAVE /HEPEUP/

C...The user's own transfer of information.
      COMMON/MYCOMM/ECM,PTMIN,INIT
      SAVE/MYCOMM/

C...Set up incoming beams.
      IDBMUP(1)= 2212
      IDBMUP(2)= -2212
      EBMUP(1) = 0.5D0*ECM 
      EBMUP(2) = 0.5D0*ECM

C...Set up the external process.
      IDWTUP   = 1
      NPRUP    = 1
      LPRUP(1) = 1001     

C...Find some reasonable maximum for external process.
C...(This is done by a cheating call to UPEVNT with special switch.)
      INIT=1
      CALL UPEVNT
      INIT=0
      XMAXUP(1) = XWGTUP  

C...Also switch on same internal process, for comparison.
      MSUB(29)=1
      CKIN(3)=PTMIN
 
      RETURN
      END
 
C*********************************************************************

C...This routine samples the phase space, evaluates the process cross 
C...and sets the colour topology and parton shower scales. 
 
      SUBROUTINE UPEVNT
 
C...All real arithmetic in double precision.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
 
C...User process event common block.
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,IDUP(MAXNUP),
     &ISTUP(MAXNUP),MOTHUP(2,MAXNUP),ICOLUP(2,MAXNUP),PUP(5,MAXNUP),
     &VTIMUP(MAXNUP),SPINUP(MAXNUP)
      SAVE /HEPEUP/

C...The user's own transfer of information.
      COMMON/MYCOMM/ECM,PTMIN,INIT
      SAVE/MYCOMM/

C...Local arrays and parameters.
      DIMENSION XPP(-25:25),XPPBAR(-25:25),TERM(20)
      DATA PI/3.141592653589793D0/
      DATA CONV/0.3894D9/

C-----------------------------------------------------------------

C...Default return value.
      XWGTUP=0D0

C...Default flag for event rejection.
      IREJ=0

C-----------------------------------------------------------------

C...The first step is to decide how the phase space should be sampled.
C...The total cross section for a 2 -> 2 process can be written
C...   sigma = Integral d(x_1) Integral d(x_2) Integral d(t-hat) 
C...           f_1(x_1,Q^2) f_2(x_2,Q^2) d(sigma-hat)/d(t-hat)
C...(see manual section 7.2), so a simple recipe would be to pick
C...x_1, x_2 and t_hat uniformly. However, the matrix element and 
C...parton distributions are peaked, so this would be inefficient.
C...In particlar, the matrix element is proportional to 1/u-hat,
C...so note that d(t-hat)= d(u-hat), and therefore rewrite as
C...   sigma = Integral d(x_1)/x_1 Integral d(x_2)/x_2 
C...           Integral d(u-hat)/u-hat  x_1 f_1(x_1,Q^2) 
C...           x_2 f_2(x_2,Q^2) u-hat d(sigma-hat)/d(t-hat).
C...(More fancy methods are used in Pythia to improve efficiency,
C...but here the principles are more important.)

C...The pTmin cut translates into cuts on x_1, x_2 and u-hat.
      X1MIN=4D0*PTMIN**2/ECM**2  
      X2MIN=X1MIN
      UHTMIN=PTMIN**2 
      UHTMAX=ECM**2      

      IF(INIT.EQ.0) THEN
C...Normally select a phase space point at random from 
C...d(x_1)/x_1 d(x_2)/x_2 d(u-hat)/u-hat.
        X1=X1MIN**PYR(0)
        X2=X2MIN**PYR(0)
        UHAT=-UHTMIN*(UHTMAX/UHTMIN)**PYR(0)

      ELSE
C...For special initialization call pick lowest point as 
C...likely to represent maximum cross section.
        X1=2D0*PTMIN/ECM
        X2=X1
        UHAT=-2D0*PTMIN**2
      ENDIF

C...Derive other kinematical quantities.
      SHAT=X1*X2*ECM**2
      THAT=-SHAT-UHAT
      PT2=THAT*UHAT/SHAT

C...By decoupling the x_1, x_2 and u-hat variables we overestimate
C...the phase space region. We therefore have to reject events 
C...outside allowed phase space by putting the cross section = 0.
C...However, to keep bookkeeping correct, one still needs to set
C...variables as if there is an event, e.g. picking the same
C...point as for the initialization call.
      IF(INIT.EQ.0.AND.PT2.LT.PTMIN**2) THEN
        IREJ=1
        X1=2D0*PTMIN/ECM
        X2=X1
        UHAT=-2D0*PTMIN**2
        SHAT=X1*X2*ECM**2
        THAT=-SHAT-UHAT
        PT2=THAT*UHAT/SHAT
      ENDIF    

C...The selection strategy for the phase space point relies on
C...a phase space volume V inside which points are sampled. The
C...expression sigma = Integral dV d(sigma)/dV =
C...                 = < d(sigma)/dV > * Integral dV
C...(with < ... > denoting (Monte Carlo sampled) mean value over V)
C...implies we need Integral dV = 
C...Integral d(x_1)/x_1 Integral d(x_2)/x_2 Integral d(u-hat)/u-hat.
      PHSPV=(-LOG(X1MIN))*(-LOG(X2MIN))*LOG(UHTMAX/UHTMIN)

C-----------------------------------------------------------------

C...The second step is to evaluate the cross section in the 
C...phase space point selected above.

C...Pick Q2 scale (which involves some arbitrariness) 
C...and evaluate alpha_em and alpha_s.
      Q2=PT2
      AEM=PYALEM(Q2)
      ALPS=PYALPS(Q2) 

C...Evaluate matrix element d(sigma-hat)/d(t-hat), except for
C...a factor e_q^2 that varies over the quark flavours.
      DSIGDT=(PI*AEM*ALPS/3D0)*(SHAT**2+UHAT**2)/(-SHAT**3*UHAT)

C...Correction factor u-hat from our choice of phase space sampling.
      DSIGDT=(-UHAT)*DSIGDT

C...Now need sum of charge-square-weighted quark*gluon distributions.
C...The relative size of these terms is also needed for the 
C...Monte Carlo selection of one specific flavour in initial state,
C...and therefore the terms are stored in an array.

C...Parton distributions (multiplied by x) of p and pbar.
      CALL PYPDFU(2212,X1,Q2,XPP)
      CALL PYPDFU(-2212,X2,Q2,XPPBAR)

C...Loop over quark flavours (up to b).
      SUM=0D0
      DO 100 IQ=1,5
C...Squared quark charge.
        EQ2=1D0/9D0
        IF(IQ.EQ.2.OR.IQ.EQ.4) EQ2=4D0/9D0
C...q from p, g from pbar.
        TERM(4*IQ-3)=EQ2*XPP(IQ)*XPPBAR(21)
C...qbar from p, g from pbar.
        TERM(4*IQ-2)=EQ2*XPP(-IQ)*XPPBAR(21)
C...g from p, q from pbar.
        TERM(4*IQ-1)=EQ2*XPP(21)*XPPBAR(IQ)
C...g from p, qbar from pbar.
        TERM(4*IQ)=EQ2*XPP(21)*XPPBAR(-IQ)
        SUM=SUM+TERM(4*IQ-3)+TERM(4*IQ-2)+TERM(4*IQ-1)+TERM(4*IQ)
  100 CONTINUE  

C...Now possible to define the SIGEV return value.
C...The CONV factor converts from GeV**(-2) to pb.
      IF(IREJ.EQ.0) XWGTUP=CONV*PHSPV*DSIGDT*SUM

C...Done when initializing for maximum cross section.
      IF(INIT.EQ.1) RETURN

C...Store scale choice etc.
      SCALUP=SQRT(Q2)
      AQEDUP=AEM
      AQCDUP=ALPS

C-----------------------------------------------------------------

C...The third step is to set up the partonic process that is selected.

C...Define number of partons - two incoming and two outgoing.
      NUP=4

C...Select one of the possible flavour terms at random.
      RANSUM=PYR(0)*SUM
      ISEL=0
  110 ISEL=ISEL+1
      RANSUM=RANSUM-TERM(ISEL)
      IF(ISEL.LT.20.AND.RANSUM.GT.0D0) GOTO 110

C...Read off quark flavour and quark side from selection.
      KFLQ=(ISEL+3)/4
      IF(ISEL-4*KFLQ.EQ.-2.OR.ISEL-4*KFLQ.EQ.0) KFLQ=-KFLQ
      ISIDE=1
      IF(ISEL-4*KFLQ.GE.-1) ISIDE=2

C...Flavour codes for entries. Note that definition of t-hat
C...means quark either is 1 and 3 or 2 and 4.
      IDUP(ISIDE)=KFLQ
      IDUP(3-ISIDE)=21
      IDUP(ISIDE+2)=KFLQ
      IDUP(5-ISIDE)=22

C...Status codes.
      ISTUP(1)=-1
      ISTUP(2)=-1
      ISTUP(3)=1
      ISTUP(4)=1

C...Mother codes.
      MOTHUP(1,1)=0
      MOTHUP(2,1)=0
      MOTHUP(1,2)=0
      MOTHUP(2,2)=0
      MOTHUP(1,3)=1
      MOTHUP(2,3)=2
      MOTHUP(1,4)=1
      MOTHUP(2,4)=2

C...Colour flow.
      IF(KFLQ.GT.0) THEN
C...Colour stretched from initial quark to gluon to final quark.
        ICOLUP(1,ISIDE)  =501
        ICOLUP(2,3-ISIDE)=501
        ICOLUP(1,3-ISIDE)=502
        ICOLUP(1,ISIDE+2)=502
        ICOLUP(2,ISIDE)  =0
        ICOLUP(2,ISIDE+2)=0
        ICOLUP(1,5-ISIDE)=0
        ICOLUP(2,5-ISIDE)=0
      ELSE
C...Ditto but with colour <--> anticolour.
        ICOLUP(2,ISIDE)  =501
        ICOLUP(1,3-ISIDE)=501
        ICOLUP(2,3-ISIDE)=502
        ICOLUP(2,ISIDE+2)=502
        ICOLUP(1,ISIDE)  =0
        ICOLUP(1,ISIDE+2)=0
        ICOLUP(1,5-ISIDE)=0
        ICOLUP(2,5-ISIDE)=0
      ENDIF  

C...Reset momenta to zero.
      DO 130 I=1,4
        DO 120 J=1,5
          PUP(J,I)=0D0
  120   CONTINUE
  130 CONTINUE

C...Masses of final state entries; initial assumed massless.
      PUP(5,3)=PYMASS(IDUP(3))
      PUP(5,4)=PYMASS(IDUP(4))

C...Four-momenta of the incoming partons simple.
      PUP(4,1)=0.5D0*X1*ECM
      PUP(3,1)=PUP(4,1)
      PUP(4,2)=0.5D0*X2*ECM
      PUP(3,2)=-PUP(4,2)

C...Energies and absolute momentum of the outgoing partons in 
C...the subsystem frame.
      RTSHAT=SQRT(X1*X2)*ECM
      PABS=0.5D0*SQRT(MAX(0D0,(RTSHAT**2-PUP(5,3)**2-
     &PUP(5,4)**2)**2-4D0*PUP(5,3)**2*PUP(5,4)**2))/RTSHAT
      PE3=0.5D0*(RTSHAT**2+PUP(5,3)**2-PUP(5,4)**2)/RTSHAT
      PE4=RTSHAT-PE3

C...Subsystem scattering angle defined neglecting quark mass.
      COSTHE=(THAT-UHAT)/SHAT
      SINTHE=SQRT(MAX(0D0,1D0-COSTHE**2))

C...Azimuthal angle at random.
      PHI=2D0*PI*PYR(0)

C...Momenta of outgoing partons in the subsystem frame.
      PUP(1,3)=PABS*SINTHE*COS(PHI)
      PUP(2,3)=PABS*SINTHE*SIN(PHI)
      PZ3=PABS*COSTHE
      PUP(1,4)=-PUP(1,3)
      PUP(2,4)=-PUP(2,3)
      PZ4=-PZ3

C...Longitudinal boost of outgoing partons to cm frame.
      BETA=(X1-X2)/(X1+X2)
      GAMMA=0.5D0*(X1+X2)/SQRT(X1*X2)
      PUP(3,3)=GAMMA*(PZ3+BETA*PE3)
      PUP(4,3)=GAMMA*(PE3+BETA*PZ3)
      PUP(3,4)=GAMMA*(PZ4+BETA*PE4)
      PUP(4,4)=GAMMA*(PE4+BETA*PZ4)

C...Done.
      RETURN
      END
