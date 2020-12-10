# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([pythia6-examples],[421],[],[pythia6-examples])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE

AC_LANG(Fortran 77)
AC_LANG(C++)

			     
PYTHIA6_DIR=no

if test -s paths.sh ; then
    eval `cat paths.sh`
fi


AC_ARG_WITH([pythia6],
	      [AS_HELP_STRING([--with-pythia6],
			      [path to pythia6 installation,
                   is taken from ./paths.sh by default]
			     )
	      ],
	      [],
	      [	with_pythia6=${PYTHIA6_DIR} ]
	   )

if test "x$with_pythia6" = "xno" ; then
	echo  "WARNING: path to Pythia 6 is not specified. Trying to guess ..."
	echo -n "Checking ../lib/ ... "
	if test -d ../lib ; then
		if test "x`ls -1 ../lib/libpythia6.*`" != "x" ; then
			echo "found: " `ls ../lib/libpythia6.*`
			with_pythia6=`pwd`/..
		else
			echo "libpythia6.* NOT found"
		fi
	else
		echo " NOT found. Please specify path to Pythia 6 installation by --with-pythia6=/path/to/pythia6"
		exit 1
	fi
fi
		

AS_IF([test "x$with_pythia6" != "xno" && test "x`ls $with_pythia6/lib/libpythia6.*`" != "x" ], [AC_SUBST(PYTHIA6_DIR,[$with_pythia6])],
      [echo "ERROR: please specify correct path to Pythia 6 installation: --with-pythia6=/path/to/pythia6/.">&2 ; exit 1 ]
     )


AC_ARG_WITH([hepmc],
	    [AS_HELP_STRING([--with-hepmc],
			     [path to HepMC installation @<:@default=no@:>@ .....
                  On CERN AFS:
                  /afs/cern.ch/sw/lcg/external/HepMC/version/platform
                 ]
			    )
	    ],
	    [],
	    [with_hepmc=no]
       )

HEPMC_VERSION=
HEPMC_VERSION_MAJOR=

if test "x$with_hepmc" != "xno" ; then
	echo -n "Checking for valid HepMC installation at ${with_hepmc} ..."
	if test "x` ls $with_hepmc/lib/libHepMC.* 2>/dev/null `" != "x" ; then
		echo "libHepMC is found" ;
		echo -n "Checking for $with_hepmc/include/ ... " 
		if test -d $with_hepmc/include ; then
			echo "found"
		else
			echo " NOT found. Check HepMC installation">&2
			exit 1
		fi
		echo -n "Checking for HepMC version ... "
		HEPMC_VERSION=[`grep -r -E "^[ ]*#define [ ]*VERSION " ${with_hepmc}/include | sed 's,^.*VERSION [ \"]*,,;s,[ \"]*$,,'`]
		HEPMC_VERSION_MAJOR=[`echo $HEPMC_VERSION | sed 's,^\([0-9][0-9]*\)\..*$,\1,'`]
		echo "$HEPMC_VERSION (major: $HEPMC_VERSION_MAJOR)"
	else
		echo " NOT found. Check HepMC installation" >&2
		exit 1
	fi
fi

AM_CONDITIONAL([HEPMC], [test "x$with_hepmc" != "xno" && test "x$HEPMC_VERSION_MAJOR" = "x1"])
AM_CONDITIONAL([HEPMC2], [test "x$with_hepmc" != "xno" && test "x$HEPMC_VERSION_MAJOR" = "x2"])

AS_IF([test "x$with_hepmc" != "xno"], [AC_SUBST(HEPMC_DIR,[$with_hepmc])],
      [echo "WARNING: path to HepMC installation is not given, test using HepMC won't be compiled."]
     )



AC_ARG_WITH([clhep],
	    [AS_HELP_STRING([--with-clhep],
			     [path to CLHEP installation @<:@default=no@:>@ ....
                  On CERN AFS:
                  /afs/cern.ch/sw/lcg/external/clhep/version/platform
                 ]
			    )
	    ],
	    [],
	    [with_clhep=no]
	   )

	
AM_CONDITIONAL([CLHEP], [test "x$with_clhep" != "xno" && test "x`ls ${with_clhep}/lib/libCLHEP.* 2>/dev/null `" != "x"])
AS_IF([test "x$with_clhep" != "xno"], [AC_SUBST(CLHEP_DIR,[$with_clhep])],
      [if test "x$HEPMC_VERSION_MAJOR" = "x1" ; then \
       		echo "ERROR: correct path to CLHEP installation is not given. Please run ./configure --with-clhep=/path/to/CLHEP/ or use HepMC >= 2.00.00">&2 ; \
		exit 1 ; \
       fi \ 
      ]
     )
		    
AC_ARG_WITH([lhapdf],
	    [AS_HELP_STRING([--with-lhapdf],
			     [path to LHAPDF installation @<:@default=no@:>@ ...
                  On CERN AFS:
                  /afs/cern.ch/sw/lcg/external/MCGenerators/lhapdf/version/platform
                 ]
			    )
	    ],
	    [],
	    [with_lhapdf=no]
	   )

AM_CONDITIONAL([LHAPDF], [test "x$with_lhapdf" != "xno"])

AS_IF([test "x$with_lhapdf" != "xno" && test "x`ls $with_lhapdf/lib/libLHAPDF.*`" != "x" ], [AC_SUBST(LHAPDF_DIR,[$with_lhapdf])],
      [test "x$with_lhapdf" = "xno" && echo "WARNING: examples with LHAPDF won't be compiled. Use --with-lhapdf=/path/to/lhapdf/." || echo "ERROR: can't find $with_lhapdf/lib/libLHAPDF.*">&2 ]
     )


AC_ARG_WITH([lhapdfsets],
            [AS_HELP_STRING([--with-lhapdfsets],
              [path to LHAPDF data (PDFsets). To be used for non-standard
               LHAPDF installations, for example in absence
               of script lhapdf-config @<:@default=no@:>@
              ]
                           )
            ],
            [],
            [with_lhapdfsets=no]
           )

rm -f ld_config.sh
rm -f ld_config.csh
(test "x$with_pythia6" != "xno" && echo "export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:${with_pythia6}/lib ; export DYLD_LIBRARY_PATH=\${DYLD_LIBRARY_PATH}:${with_pythia6}/lib ; "
 test "x$with_hepmc"   != "xno" && echo "export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:${with_hepmc}/lib   ; export DYLD_LIBRARY_PATH=\${DYLD_LIBRARY_PATH}:${with_hepmc}/lib   ; "
 test "x$with_clhep"   != "xno" && echo "export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:${with_clhep}/lib   ; export DYLD_LIBRARY_PATH=\${DYLD_LIBRARY_PATH}:${with_clhep}/lib   ; "
 test "x$with_lhapdf"  != "xno" && echo "export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:${with_lhapdf}/lib  ; export DYLD_LIBRARY_PATH=\${DYLD_LIBRARY_PATH}:${with_lhapdf}/lib  ; "
) > ld_config.sh

echo "#!/bin/csh" > ld_config.csh
echo 'if( ! $?LD_LIBRARY_PATH ) then'                 >> ld_config.csh
echo "  setenv LD_LIBRARY_PATH ${with_pythia6}/lib"   >> ld_config.csh
echo "else"                                           >> ld_config.csh
echo "  setenv LD_LIBRARY_PATH" '${LD_LIBRARY_PATH}'":${with_pythia6}/lib" >> ld_config.csh
echo "endif"                                          >> ld_config.csh
echo 'if( ! $?DYLD_LIBRARY_PATH ) then'               >> ld_config.csh
echo "  setenv DYLD_LIBRARY_PATH ${with_pythia6}/lib" >> ld_config.csh
echo "else"                                           >> ld_config.csh
echo "  setenv DYLD_LIBRARY_PATH" '${DYLD_LIBRARY_PATH}'":${with_pythia6}/lib" >> ld_config.csh
echo "endif"                                          >> ld_config.csh
test "x$with_hepmc"   != "xno" && echo "setenv LD_LIBRARY_PATH \${LD_LIBRARY_PATH}:${with_hepmc}/lib  ; setenv DYLD_LIBRARY_PATH \${DYLD_LIBRARY_PATH}:${with_hepmc}/lib  ; " >> ld_config.csh
test "x$with_clhep"   != "xno" && echo "setenv LD_LIBRARY_PATH \${LD_LIBRARY_PATH}:${with_clhep}/lib  ; setenv DYLD_LIBRARY_PATH \${DYLD_LIBRARY_PATH}:${with_clhep}/lib  ; " >> ld_config.csh
test "x$with_lhapdf"  != "xno" && echo "setenv LD_LIBRARY_PATH \${LD_LIBRARY_PATH}:${with_lhapdf}/lib ; setenv DYLD_LIBRARY_PATH \${DYLD_LIBRARY_PATH}:${with_lhapdf}/lib ; " >> ld_config.csh

if test "x$with_lhapdf" != "xno" ; then
    LHAPATH=
    if test "x$with_lhapdfsets" != "xno" ; then
       LHAPATH=$with_lhapdfsets
    else
       echo -n "Checking for lhapdf-config ..."
       if test "x`ls $with_lhapdf/bin/lhapdf-config 2>/dev/null`" != "x"  ; then
         echo " found."
         LHAPATH=`$with_lhapdf/bin/lhapdf-config --pdfsets-path 2>/dev/null`
         if test "x`ls $LHAPATH/MRST2001lo.LHgrid 2>/dev/null`" != "x" ; then
           echo " PDF sets found"
         else
           echo "but PDF sets not found"
           LHAPATH=
         fi
       fi
       if test "x$LHAPATH" = "x" ; then
        echo " PDF sets location not found with lhapdf-config. Trying to guess PDFsets location ..."
        LHAPATH=$with_lhapdf/PDFsets
        if test "x`ls $LHAPATH/MRST2001lo.LHgrid 2>/dev/null`" != "x" ; then
          echo " found"
        else
          LHAPATH=
        fi
        if test "x$LHAPATH" = "x" ; then
          LHAPATH=$with_lhapdf/../share/PDFsets
          if test "x`ls $LHAPATH/MRST2001lo.LHgrid 2>/dev/null`" != "x" ; then
            echo " found"
          else
            echo "PDFsets location not found. examples with LHAPDF won't work"
            LHAPATH=
          fi
        fi
      fi
    fi
    echo -n "Checking for PDF sets ... "
    if test "x`ls $LHAPATH/MRST2001lo.LHgrid 2>/dev/null`" != "x" ; then
       echo " found"
       echo "export LHAPATH=${LHAPATH} ;" >> ld_config.sh
       echo "setenv LHAPATH ${LHAPATH}" >> ld_config.csh
    else
       echo "NOT found. LHAPATH is not set correctly: examples with LHAPDF won't work"
    fi
fi

AC_CONFIG_SRCDIR([main61.f])
AC_CONFIG_HEADERS([config/config.h])

## Checks for programs.
AC_PROG_CXX
AC_PROG_F77

AC_PROG_LIBTOOL

AC_PROG_INSTALL
AC_PROG_LN_S


AC_CONFIG_FILES([Makefile])

if uname | grep -E "^Darwin" >&/dev/null ; then
        LDFLAGS="$LDFLAGS -Wl,-flat_namespace"
fi

AC_OUTPUT

