# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([pythia6],[427.2],[],[pythia6])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([foreign])
AC_LANG(Fortran 77)
AC_LANG(C)

AC_PREFIX_DEFAULT([`pwd`])

			     
AC_CONFIG_SRCDIR([src/pyevnt.f])
AC_CONFIG_HEADERS([config/config.h])

## Checks for programs.
AC_PROG_CC

# here is a trick to preffer
#   gfortran for gcc4 and higher
#   and g77 for gcc3 and lower
flist="gfortran g77 f77"
if [[[ "$(gcc -dumpversion | cut -d . -f 1)" < "4" ]]] ; then
  flist="g77 f77 gfortran"
fi

AC_PROG_F77([$flist])

AC_PROG_LIBTOOL

AC_PROG_INSTALL
AC_PROG_LN_S


AM_CONDITIONAL([DO_INSTALL],[test "x$prefix" != "xNONE"])

AC_CONFIG_FILES([Makefile src/Makefile dummy/Makefile pdfdummy/Makefile rootinterface/Makefile include/Makefile])

if uname | grep -E "^Darwin" >&/dev/null ; then
        LDFLAGS="$LDFLAGS -Wl,-flat_namespace"
fi

AC_OUTPUT

rm -f examples/paths.sh
echo "PYTHIA6_DIR=$prefix" > examples/paths.sh

#---- Fortran crude dependency tracking: ---

AC_DEFUN([AC_FORTRAN_DEPENDENCIES],
[[
echo ""
echo " *** Crude GENSER hack for automatic dependency tracking of Fortran source files. Is there any regular solution? *** "
echo ""
for _mkfile in `find . -name 'Makefile.am' ` ; do 
	_dir=`echo $_mkfile | sed 's,^\(.*\)Makefile.am[ ]*$,\1,'`
	_topdir=`pwd`
	cd $_dir
	echo " In $_dir: "
	_include_paths=`
	 cat Makefile.am \
	 | sed 's,$,@@ENDLINE@@,' | tr '\n' '@' | sed 's,\\\\@@ENDLINE@@@, ,g' | sed 's,@@ENDLINE@@@,#,g' \
	 | tr '#' '\n' \
	 | grep -E "=[^=]*\-I[^=]*$" \
	 |sed 's,^.*=[ ]*,,; s,^-I, ,; s, -I, ,g;'
	`
	cat Makefile.am \
	 | sed 's,$,@@ENDLINE@@,' | tr '\n' '@' | sed 's,\\@@ENDLINE@@@, ,g' | sed 's,@@ENDLINE@@@,#,g' | tr '#' '\n' \
	 | grep _SOURCES \
	 | while read str ; do
	  	_prefix_SOURCES=`echo $str | sed 's,_SOURCES[ ]*=.*$,,'`
		_source_files=`echo $str | sed 's,^.*[ ]*=[ ]*,,' | tr ' \t' '\n\n' | grep -E "\.[Ff][^\.]*$"`
		for _f in $_source_files ; do
			_includes=`
		 	 sed '/^[ \t]*[#]\{0,1\}[Ii][Nn][Cc][Ll][Uu][Dd][Ee][ \t]*/!d' < $_f \
			 | sed "s/^[ \t]*[#]\{0,1\}[Ii][Nn][Cc][Ll][Uu][Dd][Ee][ \t]*[\'\"]\([^\'\"][^\'\"]*\)[\'\"].*$/\1/" \
			 | sort | uniq
			`
		 	if test "x$_includes" != "x" ; then
				_file_basename=`echo $_f | sed 's,^.*/\([^/][^/]*\)\.[Ff][^\,]*$,\1,'`
				_file_prefix=${_prefix_SOURCES}-${_file_basename}
				_obj_file=${_file_prefix}.lo
				_obj_file1=`echo ${_file_basename} | sed 's,\.[Ff][^\.]*$,,'`.lo
				echo "   $_f: $_includes" | tr '\n' ' ' ; echo
		 		(
				  echo ""
				  echo "_${_file_prefix}_inc_files = ${_includes}" | tr '\n' ' ' ; echo
				  echo "_${_file_prefix}_inc_paths = ${_include_paths}" | tr '\n' ' ' ; echo
				  echo "_${_file_prefix}_deps := \$(shell for __f in \$(_${_file_prefix}_inc_files) ; do for __incdir in \$(_${_file_prefix}_inc_paths) ; do if test -s \$\$__incdir/\$\$__f ; then echo \$\$__incdir/\$\$__f ; break ; fi ; done ; done)"
				  echo "${_obj_file}: \$(_${_file_prefix}_deps)"
				  echo "${_obj_file1}: \$(_${_file_prefix}_deps)"
				) >> Makefile
			fi
		done
	   done
	cd $_topdir
done
echo
]]
)

AC_FORTRAN_DEPENDENCIES
