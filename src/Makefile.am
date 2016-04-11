EXTRA_DIST=	COPYING INSTALL

bin_PROGRAMS=	mb

mb_SOURCES=	\
	bayes.c \
	bayes.h \
	best.c \
	best.h \
	command.c \
	command.h \
	likelihood.c \
	likelihood.h \
	mbbeagle.c \
	mbbeagle.h \
	mcmc.c \
	mcmc.h \
	model.c \
	model.h \
	proposal.c \
	proposal.h \
	sumpt.c \
	sumpt.h \
	utils.c \
	utils.h

mb_CFLAGS=	$(SIMD_FLAGS) $(CPUEXT_FLAGS) -DHOST=$(host) -DHOSTCPU=$(host_cpu)

mb_LDADD=	$(BEAGLE_LIBS)
mb_CFLAGS+=	$(BEAGLE_CFLAGS)
