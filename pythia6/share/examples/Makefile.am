
bin_PROGRAMS = main61 main62 main63 main64 main65 main66 main67 main68 main69 main70 main71 main72 main76 main77


AM_LDFLAGS =

INCLUDES = -I$(PYTHIA6_DIR)/include
AM_FFLAGS = $(INCLUDES)

main61_LDFLAGS = $(PYTHIA6_DIR)/lib/pydata.o -L$(PYTHIA6_DIR)/lib -lpythia6 -lpythia6_dummy -lpythia6_pdfdummy

main62_LDFLAGS = $(main61_LDFLAGS)
main63_LDFLAGS = $(main61_LDFLAGS)
main64_LDFLAGS = $(main61_LDFLAGS)
main65_LDFLAGS = $(main61_LDFLAGS)
main66_LDFLAGS = $(main61_LDFLAGS)
main67_LDFLAGS = $(main61_LDFLAGS)
main68_LDFLAGS = $(main61_LDFLAGS)
main69_LDFLAGS = $(main61_LDFLAGS)
main70_LDFLAGS = $(main61_LDFLAGS)
main71_LDFLAGS = $(main61_LDFLAGS)
main72_LDFLAGS = $(main61_LDFLAGS)
main76_LDFLAGS = $(main61_LDFLAGS)
main77_LDFLAGS = $(main61_LDFLAGS)

main61_SOURCES = main61.f  
main62_SOURCES = main62.f  
main63_SOURCES = main63.f  
main64_SOURCES = main64.f  
main65_SOURCES = main65.f  
main66_SOURCES = main66.f  
main67_SOURCES = main67.f  
main68_SOURCES = main68.f  
main69_SOURCES = main69.f  
main70_SOURCES = main70.f  
main71_SOURCES = main71.f  
main72_SOURCES = main72.f  
main76_SOURCES = main76.f  
main77_SOURCES = main77.f



if HEPMC2

bin_PROGRAMS += main78 main79
INCLUDES += -I$(HEPMC_DIR)/include
main78_LDFLAGS = $(main61_LDFLAGS) -L$(HEPMC_DIR)/lib -lHepMC -lHepMCfio
main79_LDFLAGS = $(main78_LDFLAGS)
main78_SOURCES = main78.cc
main79_SOURCES = main79.cc

else # HepMC < 2 is used, need CLHEP for this:

if HEPMC

if CLHEP

bin_PROGRAMS += main78 main79
INCLUDES += -I$(HEPMC_DIR)/include -I$(CLHEP_DIR)/include
main78_LDFLAGS = $(main61_LDFLAGS) -L$(HEPMC_DIR)/lib -lHepMC -lHepMCfio -L$(CLHEP_DIR)/lib -lCLHEP
main79_LDFLAGS = $(main78_LDFLAGS)
main78_SOURCES = main78.cc
main79_SOURCES = main79.cc

else

main78 main79: main7879-msg

main7879-msg:
	@echo "ERROR: run ./configure --with-clhep=/path/to/CLHEP/ to use HepMC < 2.00.00">&2  ; exit 1

endif

else # no HepMC at all:

main78 main79: main7879-nohepmc-msg

main7879-nohepmc-msg:
	@echo "INFO: run ./configure --with-hepmc=/path/to/HepMC/ to compile $@" ; exit 1

endif

endif



if LHAPDF

bin_PROGRAMS += main81
main81_LDFLAGS = $(PYTHIA6_DIR)/lib/pydata.o -L$(PYTHIA6_DIR)/lib -lpythia6 -lpythia6_dummy -L$(LHAPDF_DIR)/lib -lLHAPDF
main81_SOURCES = main81.f

else
main81: main81-msg

main81-msg:
	@echo "Can't compile $@ : run ./configure --with-lhapdf=/path/to/LHAPDF/" ; exit 1
endif

dist-hook: 
	rm -rf `find $(distdir) -name '.svn'` `find $(distdir) -type 'd' -name 'CVS'`
