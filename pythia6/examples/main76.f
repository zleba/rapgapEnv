      PROGRAM MAIN76
C...Example program to generate stop1 pair production events at LHC, using a 
C...SuSy spectrum calculated by an external RGE code (in this case SoftSusy),
C...input via the SuSy Les Houches Accord, 
C...                                 see http://home.fnal.gov/~skands/slha.
C...DEFAULT INPUT:
C...Spectrum of couplings and masses for SPS1a from SoftSusy. You can of
C...course replace this by any other SLHA spectrum you want to use. Decays 
C...and cross sections are calculated by Pythia. If you wish to use input 
C...from an SLHA file with a decay list, the relevant line below is 
C...included but commented out. Simply uncomment it.
C...NOTE: in its default form, this program assumes the file 'softsusy.spc' 
C...to exist in the current directory (and contain an SLHA spectrum).

C--------------- PREAMBLE: COMMON BLOCK DECLARATIONS ETC -------------
C...All real arithmetic done in double precision.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C...The PYTHIA event record:
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
C...PYTHIA MSSM and subprocess common blocks
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
C...EXTERNAL statement links PYDATA on most machines.
      EXTERNAL PYDATA

C-------------------------- PYTHIA SETUP -----------------------------
C - Number of events to generate.
      NEV=5000

C - Subprocesses: q qbar -> ~t_1 ~t_1bar and g g -> ~t_1 ~t_1bar.
      MSEL=0
      MSUB(261)=1
      MSUB(264)=1

C - CM energy
      ECM=14000D0

C - Set SuSy parameters, using SLHA input.
C...Switch on SUSY MSSM input from an SLHA file.
      IMSS(1)=11
C...Open the SLHA file and tell Pythia where (on which LUN) to find it:
      LUNSPC=22
      OPEN(LUNSPC,file='softsusy.spc',status='old')
      IMSS(21)=LUNSPC
C...If the file also contains an SLHA decay table that you want to use,
C...uncomment the next line: 
C      IMSS(22)=LUNSPC
      
C - Initialize PYTHIA
      CALL PYINIT('CMS','p','p',ECM)

C - Close the SLHA file (only used during initialization).
      CLOSE(LUNSPC)

C------------------------- GENERATE EVENTS ---------------------------
      DO 100 IEV=1,NEV
        CALL PYEVNT
C...Make a print of the event record for the first event.
        IF (IEV.EQ.1) THEN 
          CALL PYLIST(2)
        ENDIF

 100  CONTINUE

C----------------------------- FINALIZE ------------------------------
C - Print info on cross sections and errors/warnings
      CALL PYSTAT(1)

      END
