SUBDIRS = src dummy pdfdummy rootinterface include

#SUBDIRS += examples

if DO_INSTALL

_dirs_to_copy  = examples

#TODO: check if examples/ contains *.o, etc.:
#
install-data-local:
	@echo "Copying " $(_dirs_to_copy) " to " $(datadir) " ... " ; \
	mkdir -p $(datadir) && \
	for _d in $(_dirs_to_copy) ; do \
		if ! test -d $$_d ; then \
			continue ; \
		fi ; \
		if test -d $(datadir)/$$_d ; then \
			rm -fr $(datadir)/$$_d ;  \
		fi ; \
		cp -r $$_d $(datadir)/$$_d ; \
	done

endif

install-exec-local:
	@echo "Creating symbolic links for archive libraries ..." ; \
	mkdir -p $(libdir)/archive/ && \
	cd $(libdir)/archive/ && \
	for _l in `grep old_library ../*.la 2>/dev/null | sed "s,^.*=',,;s,'[ ]*$$,,;"`; do \
		rm -f $$_l && ln -sf ../$$_l . ; \
	done
	
dist-hook: 
	rm -rf `find $(distdir) -name '.svn'` `find $(distdir) -type 'd' -name 'CVS'`
