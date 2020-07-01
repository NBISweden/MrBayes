EXTRA_DIST=	\
	runtests.sh \
	test1.nex \
	test2.nex

CLEANFILES=	\
	test1.nex.out \
	test2.nex.out

check: $(top_builddir)/src/mb $(EXTRA_DIST)
	./runtests.sh
