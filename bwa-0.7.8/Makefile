TARGET  := bwa

SOURCES := \
bamlite.c \
bntseq.c \
bwa.c \
bwamem.c \
bwamem_pair.c \
bwape.c \
bwase.c \
bwaseqio.c \
bwtaln.c \
bwt.c \
bwtgap.c \
bwt_gen.c \
bwtindex.c \
bwt_lite.c \
bwtsw2_aux.c \
bwtsw2_chain.c \
bwtsw2_core.c \
bwtsw2_main.c \
bwtsw2_pair.c \
fastmap.c \
is.c \
kopen.c \
kstring.c \
ksw.c \
kthread_batch.c \
kthread.c \
main.c \
malloc_wrap.c \
pemerge.c \
QSufSort.c \
utils.c

SOURCE_CPP := qa_top.cc
AAL_SRCDIR := /curr/haoyc/aalsdk_splrm-4.1.0/aalsdk_splrm-4.1.0
AAL_BUILDDIR := /curr/haoyc/aalsdk_splrm-4.1.0/aalsdk_splrm-4.1.0/mybuild

COBJECTS := $(patsubst %.c, %.o, $(SOURCES))
CCOBJECTS := $(patsubst %.cc, %.o, $(SOURCE_CPP))
LIBPATH := .

AAL_INC=\
-I$(AAL_SRCDIR)/include \
-I$(AAL_BUILDDIR)

AAL_LDADD=\
$(AAL_BUILDDIR)/aas/OSAL/libOSAL.la \
$(AAL_BUILDDIR)/aas/AASLib/libAAS.la \
$(AAL_BUILDDIR)/CCILib/libCCI.la \
$(AAL_BUILDDIR)/clp/libaalclp.la

CXX = gcc
CXXFLAGS    = -g -Wall -Wno-unused-function -O2 $(AAL_INC)
WRAP_MALLOC = -DUSE_MALLOC_WRAPPERS
LDFLAGS     = -DHAVE_PTHREAD $(WRAP_MALLOC) -lm -lz -lpthread -lstdc++

all: $(TARGET)

$(TARGET): $(COBJECTS) $(CCOBJECTS)
	libtool --mode=link g++ $(CFLAGS) *.o -o $(TARGET) $(LDFLAGS) $(AAL_LDADD)

$(COBJECTS): %.o: %.c
	$(CXX) -c $(CXXFLAGS) -o $@ $<

$(CCOBJECTS): %.o: %.cc
	gcc -c -x c++ $(CXXFLAGS) -o $@ $<

clean:
	rm -f *.o $(TARGET)

