C...Sample main program for technicolour production at the Tevatron.

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

C...Local arrays.
      DIMENSION NCHN(12)
      DATA NCHN/12*0/ 

C-----------------------------------------------------------------

C...First section: initialization.
 
C...Number of events and cm energy. 
      NEV=1000
      ECM=2000D0

C...Possibility to set masses freely:
C...pi_tech0
      PMAS(PYCOMP(KTECHN+111),1)=100D0
C...pi_tech+-.
      PMAS(PYCOMP(KTECHN+211),1)=100D0
C...rho_tech0
      PMAS(PYCOMP(KTECHN+113),1)=200D0
C...rho_tech+-
      PMAS(PYCOMP(KTECHN+213),1)=200D0
C...omega_tech0
      PMAS(PYCOMP(KTECHN+223),1)=200D0

C...Parameters of technicolour scenario.
C...sin(chi).
c      PARP(141)=1D0/3D0
C...F_T.
c      PARP(142)=82D0
C...Q_U (Q_D = Q_U - 1 always implied). 
c      PARP(143)=1.0D0
C...N_TC.
c      PARP(144)=4D0 
C...M_T.
c      PARP(145)=200D0   

C...Technicolour processes.
      MSEL=0
C...rho_tech0.
      MSUB(191)=1
C...rho_tech+-. 
      MSUB(192)=1
C...omega_tech.
      MSUB(193)=1
C...rho_tech/omega_tech with full interference for one flavour.
C...Supersedes processes 191 and 193 for fermion pair production;
C...to include both at the same time would be doublecounting.
c      MSUB(194)=1
C...Pick flavour of above process.
c      KFPR(194,1)=13 

C...Kinematical cuts in mass - absolutely necessary for process
C...194; optional for the others.
      CKIN(1)=150D0
c      CKIN(2)=210D0

C...Force decay to muon pair only of rho_T and omega_T.
C...(Of interest when comparing with process 194).
c      KCRHOT=PYCOMP(KTECHN+113)
c      DO 100 IDC=MDCY(KCRHOT,2),MDCY(KCRHOT,2)+MDCY(KCRHOT,3)-1
c      IF(IABS(KFDP(IDC,1)).NE.13) MDME(IDC,1)=MIN(0,MDME(IDC,1))
c 100  CONTINUE
c      KCOMET=PYCOMP(KTECHN+223)
c      DO 110 IDC=MDCY(KCOMET,2),MDCY(KCOMET,2)+MDCY(KCOMET,3)-1
c      IF(IABS(KFDP(IDC,1)).NE.13) MDME(IDC,1)=MIN(0,MDME(IDC,1))
c 110  CONTINUE

C...Switch off unnecessary aspects: initial and  final state
C...showers, multiple interactions, hadronization.
C...(Optional for faster simulation of the parton-level 
C...processes only.)
c      MSTP(61)=0
c      MSTP(71)=0
c      MSTP(81)=0
c      MSTP(111)=0 

C...Initialize.
      CALL PYINIT('CMS','p','pbar',ECM)

C...List table of resonance decay channels.
      CALL PYSTAT(2)
 
C...Book histograms.
      CALL PYBOOK(1,'mass distribution rho_tech0',100,0D0,250D0)
      CALL PYBOOK(2,'mass distribution rho_tech+-',100,0D0,250D0)
      CALL PYBOOK(3,'mass distribution omega_tech0',100,0D0,250D0)
      CALL PYBOOK(4,'mass distribution pi_tech0',100,0D0,150D0)
      CALL PYBOOK(5,'mass distribution pi_tech+-',100,0D0,150D0)
      CALL PYBOOK(6,'mass distribution Z0',100,0D0,150D0)
      CALL PYBOOK(7,'mass distribution W+-',100,0D0,150D0)
      CALL PYBOOK(8,'mass distribution (rho+omega)_tech0',
     &100,0D0,250D0)
      CALL PYBOOK(11,'cos(theta*) proc 191',50,-1D0,1D0) 
      CALL PYBOOK(12,'cos(theta*) proc 192',50,-1D0,1D0) 
      CALL PYBOOK(13,'cos(theta*) proc 193',50,-1D0,1D0) 
      CALL PYBOOK(14,'cos(theta*) proc 194',50,-1D0,1D0) 

C-----------------------------------------------------------------

C...Second section: event loop.

C...Generate events and look at first few.
      DO 200 IEV=1,NEV
        CALL PYEVNT
        IF(IEV.LE.2) CALL PYLIST(1) 
        ISUB=MSTI(1)

C...Fill masses of resonances.
        DO 120 I=7,9
          KFA=IABS(K(I,2))
          IKF=0
          IF(KFA.EQ.KTECHN+113) IKF=1
          IF(KFA.EQ.KTECHN+213) IKF=2
          IF(KFA.EQ.KTECHN+223) IKF=3
          IF(KFA.EQ.KTECHN+111) IKF=4
          IF(KFA.EQ.KTECHN+211) IKF=5
          IF(KFA.EQ.23) IKF=6
          IF(KFA.EQ.24) IKF=7
          IF(IKF.NE.0) CALL PYFILL(IKF,P(I,5),1D0)
  120   CONTINUE
        IF(ISUB.EQ.194) CALL PYFILL(8,PARI(13),1D0)

C...Fill decay angle of primary resonances in rest frame.
        CTHE=PARI(41)
        CALL PYFILL(ISUB-180,CTHE,1D0)

C...Fill relative probability of decay channels.
        KDA1=MIN(IABS(K(8,2)),IABS(K(9,2)))
        KDA2=MAX(IABS(K(8,2)),IABS(K(9,2)))
        IF(ISUB.EQ.191) THEN
          IF(KDA2.EQ.24) THEN
            NCHN(1)=NCHN(1)+1
          ELSEIF(KDA1.EQ.24) THEN
            NCHN(2)=NCHN(2)+1
          ELSEIF(KDA1.EQ.KTECHN+211) THEN
            NCHN(3)=NCHN(3)+1
          ELSE
            NCHN(4)=NCHN(4)+1
          ENDIF
        ELSEIF(ISUB.EQ.192) THEN
          IF(KDA2.EQ.24) THEN 
            NCHN(5)=NCHN(5)+1
          ELSEIF(KDA1.EQ.24) THEN
            NCHN(6)=NCHN(6)+1
          ELSEIF(KDA1.EQ.23) THEN
            NCHN(7)=NCHN(7)+1
          ELSEIF(KDA1.EQ.KTECHN+111) THEN     
            NCHN(8)=NCHN(8)+1
          ELSE
            NCHN(9)=NCHN(9)+1
          ENDIF
        ELSEIF(ISUB.EQ.193) THEN
          IF(KDA1.EQ.22) THEN
            NCHN(10)=NCHN(10)+1
          ELSE
            NCHN(11)=NCHN(11)+1
          ENDIF
        ELSE
          NCHN(12)=NCHN(12)+1 
        ENDIF         

C...End of loop over events.
  200 CONTINUE

C-----------------------------------------------------------------

C...Third section: produce output and end.

C...Cross section table.
      CALL PYSTAT(1)

C...Print cross sections by resonance decay channel. 
      FAC=PARI(2)*1D9
      WRITE(6,'(A,I8,A,F10.3,A)') ' # rho_T0   -> W+ W-        =',
     &NCHN(1),' =',FAC*NCHN(1),' pb'
      WRITE(6,'(A,I8,A,F10.3,A)') ' # rho_T0   -> W+- pi_T-+   =', 
     &NCHN(2),' =',FAC*NCHN(2),' pb'
      WRITE(6,'(A,I8,A,F10.3,A)') ' # rho_T0   -> pi-T+ pi_T-  =',
     &NCHN(3),' =',FAC*NCHN(3),' pb'
      WRITE(6,'(A,I8,A,F10.3,A)') ' # rho_T0   -> f fbar       =',
     &NCHN(4),' =',FAC*NCHN(4),' pb'
      WRITE(6,'(A,I8,A,F10.3,A)') ' # rho_T+-  -> W+- Z0       =', 
     &NCHN(5),' =',FAC*NCHN(5),' pb'
      WRITE(6,'(A,I8,A,F10.3,A)') ' # rho_T+-  -> W+- pi_T0    =', 
     &NCHN(6),' =',FAC*NCHN(6),' pb'
      WRITE(6,'(A,I8,A,F10.3,A)') ' # rho_T+-  -> pi_T+- Z0    =',
     &NCHN(7),' =',FAC*NCHN(7),' pb'
      WRITE(6,'(A,I8,A,F10.3,A)') ' # rho_T+-  -> pi_T+- pi_T0 =', 
     &NCHN(8),' =',FAC*NCHN(8),' pb'
      WRITE(6,'(A,I8,A,F10.3,A)') ' # rho_T+-  -> f fbar''      =', 
     &NCHN(9),' =',FAC*NCHN(9),' pb'
      WRITE(6,'(A,I8,A,F10.3,A)') ' # omega_T0 -> gamma pi_T0  =', 
     &NCHN(10),' =',FAC*NCHN(10),' pb'
      WRITE(6,'(A,I8,A,F10.3,A)') ' # omega_T0 -> f fbar       =', 
     &NCHN(11),' =',FAC*NCHN(11),' pb'
      WRITE(6,'(A,I8,A,F10.3,A)') ' anything else              =', 
     &NCHN(12),' =',FAC*NCHN(12),' pb'

C...Histograms.
      CALL PYHIST

      END
