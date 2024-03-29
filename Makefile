SHELL=/bin/sh


########################################################
# Makefile flags and defintions

MAKEFILE     = Makefile
MAKEFILES    = $(MAKEFILE) Makefile.site Makefile.base $(if $(LIBONLY),,Makefile.test)

include Makefile.site
include Makefile.base
include Makefile.setup
ifdef LIBONLY
else
include Makefile.test
endif


# Default shell commands
RM ?= /bin/rm

# Use C++11 standard, flags differ by compiler
# -MMD generates a dependecy list for each file as a side effect
ifeq ($(CXXCOMPNAME),gnu)
   CXXFLAGS_STD = -std=c++11
   DEPFLAG = -MMD
else ifeq ($(CXXCOMPNAME), pgi)
   CXXFLAGS_STD = -std=c++11
   DEPFLAG = -MMD
else ifeq ($(CXXCOMPNAME), ibm)
   CXXFLAGS_STD = -std=c++11
   DEPFLAG = -MMD
else ifeq ($(CXXCOMPNAME), llvm)
   CXXFLAGS_STD = -std=c++11
   DEPFLAG = -MMD
else
   $(info $(CXXCOMPNAME) compiler not yet supported.)
endif


# Combine all compiler and linker flags
ifeq ($(DEBUG),true)
CXXFLAGS = $(CXXFLAGS_STD) $(CXXFLAGS_DEBUG) -I$(BUILDDIR) $(CXXFLAGS_BASE) \
           $(CXXFLAGS_TEST_DEBUG) $(CXXFLAGS_AMREX)
else
CXXFLAGS = $(CXXFLAGS_STD) $(CXXFLAGS_PROD) -I$(BUILDDIR) $(CXXFLAGS_BASE) \
           $(CXXFLAGS_TEST_PROD) $(CXXFLAGS_AMREX)
endif
LDFLAGS  = -L$(LIB_BITTREE) -lbittree $(LDFLAGS_TEST) $(LDFLAGS_STD)


# Add code coverage flags
ifeq ($(CODECOVERAGE), true)
CXXFLAGS += $(CXXFLAGS_COV)
LDFLAGS  += $(LDFLAGS_COV)
endif


# List of sources, objects, and dependencies
C_SRCS    = $(SRCS_BASE) $(SRCS_TEST)
C_OBJS    = $(addsuffix .o, $(basename $(notdir $(C_SRCS))))
DEPS      = $(C_OBJS:.o=.d)

OBJS_TEST = $(addsuffix .o, $(basename $(notdir $(SRCS_TEST))))
OBJS_BASE = $(addsuffix .o, $(basename $(notdir $(SRCS_BASE))))

# Use vpath as suggested here: http://make.mad-scientist.net/papers/multi-architecture-builds/#single
vpath %.cpp $(sort $(dir $(C_SRCS)))


##########################################################
# Makefile commands:

.PHONY: default all clean library test install
default: $(if $(LIBONLY), libbittree.a, $(BINARYNAME))
all:     $(if $(LIBONLY), libbittree.a, $(BINARYNAME))
library: libbittree.a
ifdef LIBONLY
test:
else
test: $(BINARYNAME)
	./$(BINARYNAME)
endif

# If code coverage is being build into the test, remove any previous gcda files to avoid conflict.
$(BINARYNAME): $(OBJS_TEST) $(MAKEFILES) libbittree.a
ifeq ($(CODECOVERAGE), true)
	$(RM) -f *.gcda
endif
	$(CXXCOMP) -o $(BINARYNAME) $(OBJS_TEST) $(LDFLAGS)

%.o: %.cpp $(MAKEFILES)
	$(CXXCOMP) -c $(DEPFLAG) $(CXXFLAGS) -o $@ $<

libbittree.a: $(OBJS_BASE) $(MAKEFILES)
	ar -rcs $@ $(OBJS_BASE)

install:
ifdef LIBONLY
	$(info Creating prefix: $(LIB_BITTREE_PREFIX))
	@mkdir -p $(LIB_BITTREE_PREFIX)
	@mkdir -p $(LIB_BITTREE_PREFIX)/include
	@mkdir -p $(LIB_BITTREE_PREFIX)/lib

	$(info Installing library...)
	@cp libbittree.a $(LIB_BITTREE_PREFIX)/lib

	$(info Installing headers...)
	@cp setup.log $(LIB_BITTREE_PREFIX)
	@cp Bittree_constants.h $(LIB_BITTREE_PREFIX)/include
	@for filename in $(HEADERS_BASE); do \
	    cp $$filename $(LIB_BITTREE_PREFIX)/include; \
	    done
	$(info Success!)
endif

# Clean removes all intermediate files
clean:
	$(RM) -f *.o
	$(RM) -f *.d
	$(RM) -f *.a
ifeq ($(CODECOVERAGE), true)
	$(RM) -f *.gcno
	$(RM) -f *.gcda
endif
	$(RM) -f lcov_temp.info


.PHONY: coverage
coverage:
ifeq ($(CODECOVERAGE), true)
	$(LCOV) -o lcov_temp.info -c -d . -b $(BASEDIR) --no-external
	$(GENHTML)  -o Coverage_Report lcov_temp.info
else
	$(info Include --coverage in your setup line to enable code coverage.)
endif


# Include dependencies generated by compiler
-include $(DEPS)

