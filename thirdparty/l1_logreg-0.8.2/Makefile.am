
ACLOCAL_AMFLAGS = -I config

AUTOMAKE_OPTIONS = gnu

SUBDIRS = src_c examples doc util_m

dist_noinst_SCRIPTS = bootstrap test_script

# binary only in distribution
dist-bin:
	rm -rf /tmp/dist
	make install prefix=/tmp/dist
	mv /tmp/dist/bin /tmp/dist/@PACKAGE@-@VERSION@-@host@
	tar -C /tmp/dist/ -cvf @PACKAGE@-@VERSION@-@host@.tar @PACKAGE@-@VERSION@-@host@
	gzip -9 @PACKAGE@-@VERSION@-@host@.tar

# binary and examples in distribution
dist-bin-deprecated:
	rm -rf /tmp/dist
	make install prefix=/tmp/dist/@PACKAGE@-@VERSION@-@host@
	cp -r dist_bin/* /tmp/dist/@PACKAGE@-@VERSION@-@host@
	tar -C /tmp/dist/ -cvf @PACKAGE@-@VERSION@-@host@.tar @PACKAGE@-@VERSION@-@host@
	gzip -9 @PACKAGE@-@VERSION@-@host@.tar

# make webpage directory
web:
	make doc
	cp -r doc/html web
	mkdir web/download
	cp @PACKAGE@-@VERSION@*.tar.gz web/download/
	rm -rf web/.svn

doc:
if DOXYGEN_OK
	doxygen doc/document.conf
else
	@echo
	@echo "doxygen is required to build documentations"
	@echo
endif

.PHONY: doc


