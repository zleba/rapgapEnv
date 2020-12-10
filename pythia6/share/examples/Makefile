
_main77-main77.f_inc_files = aldata.f 
_main77-main77.f_inc_paths =  $(PYTHIA6_DIR)/include  $(HEPMC_DIR)/include  $(HEPMC_DIR)/include $(CLHEP_DIR)/include 
_main77-main77.f_deps := $(shell for __f in $(_main77-main77.f_inc_files) ; do for __incdir in $(_main77-main77.f_inc_paths) ; do if test -s $$__incdir/$$__f ; then echo $$__incdir/$$__f ; break ; fi ; done ; done)
main77-main77.f.lo: $(_main77-main77.f_deps)
main77.lo: $(_main77-main77.f_deps)
