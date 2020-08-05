TEST_EXTENSIONS=	.nex

NEX_LOG_DRIVER=	\
	env AM_TAP_AWK='$(AWK)' $(SHELL) \
	$(top_srcdir)/autotools/tap-driver.sh

NEX_LOG_COMPILER=	./runtests.sh

TESTS=	test1.nex \
	test2.nex

EXTRA_DIST=	\
	$(TESTS) \
	runtests.sh
