APPLICATION_CFLAGS  = -DNDEBUG -D 'FN_XBODY="xdmf_body.partxmf"'
APPLICATION_LDFLAGS =
APPLICATION_INCLUDE = 
APPLICATION_CXX     = H5Convert.C XdmfGen.C attrvstr.c

include ../make/path_include
include $(CRONOS_ROOT_DIR)/cronos/make/Makefile


# Set default rule
.DEFAULT_GOAL := converter_all

.PHONY: converter_all
converter_all: Converter XdmfGen

BUILD_DataReader = $(BUILD_DIR)/cronos/DataReader

# The recipe to compile %.o files is provided by the base makefile under $(CRONOS_ROOT_DIR)/cronos/make/Makefile
# The base makefile will translate every %.C file provided by APPLICATION_CXX to its %.o file
# They are going to be compiled under $(BUILD_DataReader), e.g. ./Linux-amd64/cronos/DataReader
# We just have to specify the linking recipes

Attr.test: $(BUILD_DIR)/Attr.test

$(BUILD_DIR)/Attr.test: $(BUILD_DataReader)/attrvstr.o
	$(CXX) $(CXXFLAGS_FINAL) $< -o $@

Converter: $(BUILD_DIR)/Converter

$(BUILD_DIR)/Converter: $(BUILD_DataReader)/H5Convert.o $(BUILD_DIR)/cronos/generic/Hdf5File_cbase.o $(BUILD_DIR)/cronos/generic/buildinfo.o
	$(CXX) $^ $(CXXFLAGS_FINAL) $(LDFLAGS) $(LOADLIBES) -o $@

XdmfGen: $(BUILD_DIR)/XdmfGen

$(BUILD_DIR)/XdmfGen: $(BUILD_DataReader)/XdmfGen.o
	$(CXX) $^ $(CXXFLAGS_FINAL) $(LDFLAGS) $(LOADLIBES) -o $@
