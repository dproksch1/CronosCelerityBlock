X_LIB_ROOT_DIR = /home/kissmrbu/src/CronosGit/CronosNumLib
CRONOS_ROOT_DIR = /home/kissmrbu/src/CronosGit/CronosCode.master-beta

include $(X_LIB_ROOT_DIR)/makeinclude/include.mak

PATH_GENERIC = $(CRONOS_ROOT_DIR)/cronos/generic


OFLAGS = $(X_OFLAGS) -Wall -fexceptions 

CXXFLAGS = $(OFLAGS) -I.. -I. $(X_INC) -I$(PATH_GENERIC)


LDFLAGS = $(X_LIB)

LOADLIBES = -lutil_mt -lmatrix_mt -lgsl -lgslcblas -lm\
	$(X_COMPLEX_LIB) $(X_HDF5_LIBS) 


all: TestKeplerCart TestKeplerCyl orbitTest

orbitTest: orbitTest.C ../Orbit.H
	gcc $< -lm -lstdc++ -o $@



OBJS_KEPLER_BASE = KeplerOrbit.o normalisation.o units.o CException.o
OBJS_KEPLER_CART = $(OBJS_KEPLER_BASE) testKeplerCart.o
OBJS_KEPLER_CYL = $(OBJS_KEPLER_BASE) testKeplerCyl.o

%.o: %.C
	$(CXX) $(CXXFLAGS) $(CXXINC) -c $^

KeplerOrbit.o: ../KeplerOrbit.C
	$(CXX) $(CXXFLAGS) $(CXXINC) -c $^

normalisation.o: $(PATH_GENERIC)/normalisation.C
	$(CXX) $(CXXFLAGS) $(CXXINC) -c $^

units.o: $(PATH_GENERIC)/units.C
	$(CXX) $(CXXFLAGS) $(CXXINC) -c $^

CException.o: $(PATH_GENERIC)/CException.C
	$(CXX) $(CXXFLAGS) $(CXXINC) -c $^


TestKeplerCart: $(OBJS_KEPLER_CART)
	$(CXX) $^ $(LDFLAGS) $(LOADLIBES) $(CXXLIBS) -o $@

TestKeplerCyl: $(OBJS_KEPLER_CYL)
	$(CXX) $^ $(LDFLAGS) $(LOADLIBES) $(CXXLIBS) -o $@

clean:
	rm -rf *.o TestKeplerCart
