include /home/dproksch/Documents/master/master_thesis/repos/CronosNumLib/makeinclude/include.mak

OFLAGS  = $(X_OFLAGS) $(CFLAGS) -DNDEBUG

CXXFLAGS = $(OFLAGS) -I. $(X_INC) $(X_NO_AUTO_INCLUDE)

MATRIX_SOURCES = matrix.C matrix_1d_inc.C matrix_2d_inc.C matrix_3d_inc.C \
	         num_matrix_1d_inc.C num_matrix_2d_inc.C num_matrix_3d_inc.C \

MATRIX_HEADERS = matrix.H matrix_1d_inc.H matrix_2d_inc.H matrix_3d_inc.H \
	         matrix_4d_inc.H matrix_5d_inc.H \
		 num_matrix_1d_inc.H num_matrix_2d_inc.H num_matrix_3d_inc.H \
		 num_matrix_4d_inc.H num_matrix_5d_inc.H

BOUNDARY_SOURCES = boundary.C num_boundary_1d_inc.C num_boundary_2d_inc.C \
	           num_boundary_3d_inc.C

BOUNDARY_HEADERS = boundary.H num_boundary_1d_inc.H num_boundary_2d_inc.H \
	           num_boundary_3d_inc.H

HEADERS = array.H fftmatrix.H \
	$(MATRIX_HEADERS) $(BOUNDARY_HEADERS)

all: libmatrix_mt.a

array.o: array.C array.H

matrix.o: $(MATRIX_SOURCES) $(MATRIX_HEADERS)

boundary.o: $(BOUNDARY_SOURCES) $(BOUNDARY_HEADERS)

libmatrix_mt.a: array.o matrix.o boundary.o

install: install_libmatrix_mt.a install_libmatrix_mt_headers

install_libmatrix_mt.a: libmatrix_mt.a 
	$(INSTALL) -m $(MOD_DATA) $^ $(X_LIB_DIR)

install_libmatrix_mt_headers: $(HEADERS)
	$(INSTALL) -m $(MOD_DATA) $^ $(X_INC_DIR)

clean:
	$(X_CLEAN) 


