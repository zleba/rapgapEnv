PYTHIA6_DUMMY_SOURCES   =  \
fhhiggscorr.f \
fhsetflags.f \
fhsetpara.f \
pyevwt.f \
pykcut.f \
pytaud.f \
pytime.f \
ssmssm.f \
sugra.f \
upevnt.f \
upinit.f \
upveto.F \
visaje.f 

# No automatic dependency detection yet:
#
upveto.lo: $(top_srcdir)/include/hepevt.inc

lib_LTLIBRARIES = libpythia6_dummy.la
INCLUDES = -I$(top_srcdir)/include
libpythia6_dummy_la_SOURCES = $(PYTHIA6_DUMMY_SOURCES)

