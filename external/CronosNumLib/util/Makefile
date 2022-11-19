include /home/dproksch/master_thesis/CronosCelerityBlock/external/CronosNumLib/makeinclude/include.mak

OFLAGS  = $(X_OFLAGS) -DNDEBUG
CXXFLAGS = $(OFLAGS) $(X_INC)

all: libutil_mt.a

ParameterFileReader.o: ParameterFileReader.C ParameterFileReader.H

CronosOstream.o: CronosOstream.C CronosOstream.H

util.o: util.C util.H

libutil_mt.a: util.o ParameterFileReader.o CronosOstream.o

install: install_lib install_headers

install_lib: libutil_mt.a 
	$(INSTALL) -m $(MOD_DATA) $^ $(X_LIB_DIR)

install_headers: util.H complex.H ParameterFileReader.H CronosOstream.H
	$(INSTALL) -m $(MOD_DATA) $^ $(X_INC_DIR)

clean:
	$(X_CLEAN)

