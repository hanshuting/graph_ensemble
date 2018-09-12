
AUTOMAKE_OPTIONS=gnu

bin_PROGRAMS = l1_logreg_train l1_logreg_classify l1_logreg_regpath
l1_logreg_train_SOURCES = train.c ${common_SOURCES}
l1_logreg_classify_SOURCES = classify.c ${common_SOURCES}
l1_logreg_regpath_SOURCES = regpath.c ${common_SOURCES}


common_SOURCES = l1_logreg.c util.c pcg.c mmio.c dmatrix.c \
                 l1_logreg.h util.h pcg.h mmio.h dmatrix.h \
		 def.c def.h pkgdef.h blas.h lapack.h

man_MANS = 
EXTRA_DIST = 
