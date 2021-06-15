include $(X_LIB_ROOT_DIR)/makeinclude/include.mak

# Global structure
# use $(shell pwd) instead of $(PWD), the latter is wrong in some cases
# e.g. referring to the calling script instead of the make location
BUILD_DIR = $(shell pwd)/$(X_OSTYPE)


# --------------------------------------------------------------------------------------
#	GLOBAL BUILD INFOS
# --------------------------------------------------------------------------------------

# Getting infos
CRONOS_GIT_VERSION 	:= $(shell cd $(CRONOS_ROOT_DIR) && git --no-pager describe --dirty --always --tags)
CRONOS_GIT_COMMIT  	:= $(shell cd $(CRONOS_ROOT_DIR) && git rev-parse --verify HEAD)
APP_GIT_REPO		:= $(shell git remote get-url --all origin)
APP_GIT_VERSION 	:= $(shell git --no-pager describe --dirty --always --tags)
APP_GIT_COMMIT  	:= $(shell git rev-parse --verify HEAD)
BUILD_DATE  		:= $(shell date --iso=hours)
GLOBAL_INFOS := \
	-DCRONOS_ROOT_DIR=\"$(CRONOS_ROOT_DIR)\" \
	-DCRONOS_GIT_VERSION=\"$(CRONOS_GIT_VERSION)\" \
	-DCRONOS_GIT_COMMIT=\"$(CRONOS_GIT_COMMIT)\" \
	-DAPP_GIT_REPO=\"$(APP_GIT_REPO)\" \
	-DAPP_GIT_VERSION=\"$(APP_GIT_VERSION)\" \
	-DAPP_GIT_COMMIT=\"$(APP_GIT_COMMIT)\" \
	-DBUILD_DATE=\"$(BUILD_DATE)\"

# --------------------------------------------------------------------------------------
#	BUILD FLAGS
# --------------------------------------------------------------------------------------

# General compilation and linking flags
# The variables APPLICATION_CFLAGS, APPLICATION_LDFLAGS, APPLICATION_INCLUDE, APPLICATION_CXX
# are intended to be supplied by the user as needed

OFLAGS = $(X_OFLAGS) $(APPLICATION_OFLAGS) $(APPLICATION_CFLAGS) -Wall -fexceptions
APPLICATION_INCLUDE += ./
_APPLICATION_INCLUDE = $(foreach d, $(APPLICATION_INCLUDE), -I$d)

# If one of those flags change a recompilation of the whole code is triggered
CXXFLAGS_RECOMPILE = \
	$(OFLAGS) \
	$(X_INC) \
	-I$(CRONOS_ROOT_DIR)/cronos/generic \
	-I$(CRONOS_ROOT_DIR)/cronos/InterfaceFunctions \
	-I$(CRONOS_ROOT_DIR)/cronos/RiemannSolvers \
	-I$(CRONOS_ROOT_DIR)/cronos/RiemannSolvers/MHD \
	-I$(CRONOS_ROOT_DIR)/cronos/RiemannSolvers/HD \
	-I$(CRONOS_ROOT_DIR)/cronos/Limiter \
	-I$(CRONOS_ROOT_DIR)/cronos/Geom \
	-I$(CRONOS_ROOT_DIR)/cronos/Util \
	$(X_NO_AUTO_INCLUDE)

# Final compilation flags
CXXFLAGS_FINAL =\
	$(CXXFLAGS_RECOMPILE) \
	$(_APPLICATION_INCLUDE) \

# Linking flags
LDFLAGS = $(X_LIB) $(APPLICATION_LDFLAGS)
LOADLIBES = -lutil_mt -lmatrix_mt -lm -lgsl -lgslcblas $(X_COMPLEX_LIB) $(X_HDF5_LIBS)


# --------------------------------------------------------------------------------------
#	SOURCE FILES
# --------------------------------------------------------------------------------------

# Source files that need to be compiled
PROJ_CXX += \
	$(CRONOS_ROOT_DIR)/cronos/Limiter/limiter.C\
	$(CRONOS_ROOT_DIR)/cronos/generic/Hdf5File_cbase.C\
	$(CRONOS_ROOT_DIR)/cronos/generic/buildinfo.C\
	$(CRONOS_ROOT_DIR)/cronos/generic/movie.C\
	$(CRONOS_ROOT_DIR)/cronos/generic/CException.C\
	$(CRONOS_ROOT_DIR)/cronos/generic/gridParam.C\
	$(CRONOS_ROOT_DIR)/cronos/generic/eos.C\
	$(CRONOS_ROOT_DIR)/cronos/generic/normalisation.C\
	$(CRONOS_ROOT_DIR)/cronos/generic/units.C\
	$(CRONOS_ROOT_DIR)/cronos/generic/gridgen.C\
	$(CRONOS_ROOT_DIR)/cronos/generic/data.C\
	$(CRONOS_ROOT_DIR)/cronos/generic/multifluid.C\
	$(CRONOS_ROOT_DIR)/cronos/generic/specific.C\
	$(CRONOS_ROOT_DIR)/cronos/generic/grid.C\
	$(CRONOS_ROOT_DIR)/cronos/generic/gridfunc.C\
	$(CRONOS_ROOT_DIR)/cronos/generic/timestepping.C\
	$(CRONOS_ROOT_DIR)/cronos/generic/randgen.C\
	$(CRONOS_ROOT_DIR)/cronos/generic/proj.C\
	$(CRONOS_ROOT_DIR)/cronos/generic/user_equations.C\
	$(CRONOS_ROOT_DIR)/cronos/generic/fluid.C\
	$(CRONOS_ROOT_DIR)/cronos/InterfaceFunctions/problem.C\
	$(CRONOS_ROOT_DIR)/cronos/InterfaceFunctions/singlestep.C\
	$(CRONOS_ROOT_DIR)/cronos/Geom/sources.C\
	$(CRONOS_ROOT_DIR)/cronos/RiemannSolvers/reconst_const.C\
	$(CRONOS_ROOT_DIR)/cronos/RiemannSolvers/transformations.C\
	$(CRONOS_ROOT_DIR)/cronos/RiemannSolvers/save_data.C\
	$(CRONOS_ROOT_DIR)/cronos/RiemannSolvers/reconst_2nd.C\
	$(CRONOS_ROOT_DIR)/cronos/RiemannSolvers/fieldLists.C\
	$(CRONOS_ROOT_DIR)/cronos/RiemannSolvers/reconst.C\
	$(CRONOS_ROOT_DIR)/cronos/RiemannSolvers/fields_local.C\
	$(CRONOS_ROOT_DIR)/cronos/RiemannSolvers/HD/hllc_hd.C\
	$(CRONOS_ROOT_DIR)/cronos/RiemannSolvers/HD/PhysFluxesHD.C\
	$(CRONOS_ROOT_DIR)/cronos/RiemannSolvers/HD/hll_hd.C\
	$(CRONOS_ROOT_DIR)/cronos/RiemannSolvers/HD/vChar_hd.C\
	$(CRONOS_ROOT_DIR)/cronos/RiemannSolvers/RiemannSolver.C\
	$(CRONOS_ROOT_DIR)/cronos/RiemannSolvers/MHD/DissipationMHD.C\
	$(CRONOS_ROOT_DIR)/cronos/RiemannSolvers/MHD/ConstrainedTransport.C\
	$(CRONOS_ROOT_DIR)/cronos/RiemannSolvers/MHD/CTStone.C\
	$(CRONOS_ROOT_DIR)/cronos/RiemannSolvers/MHD/CTLondrilloDelZanna.C\
	$(CRONOS_ROOT_DIR)/cronos/RiemannSolvers/MHD/vChar.C\
	$(CRONOS_ROOT_DIR)/cronos/RiemannSolvers/MHD/hll.C\
	$(CRONOS_ROOT_DIR)/cronos/RiemannSolvers/MHD/hlld.C\
	$(CRONOS_ROOT_DIR)/cronos/RiemannSolvers/MHD/PhysFluxesMHD.C\
	$(CRONOS_ROOT_DIR)/cronos/RiemannSolvers/MHD/bc_emf.C\
	$(CRONOS_ROOT_DIR)/cronos/RiemannSolvers/cfl.C\
	$(CRONOS_ROOT_DIR)/cronos/RiemannSolvers/changes.C\
	$(CRONOS_ROOT_DIR)/cronos/RiemannSolvers/reconst_WENO.C
#	$(CRONOS_ROOT_DIR)/cronos/RiemannSolvers/singlestep_multifluid.C\


# Users can add files in the setup Makefile by just defining/appending to APPLICATION_CXX
APPLICATION_CXX +=	$(realpath ./modules.C)

# Use local files if supplied
# Otherwise fall back on default ones
ifneq ($(findstring solver.C, $(wildcard ./*.C)), )
APPLICATION_CXX += $(realpath ./solver.C)
else
APPLICATION_CXX += $(CRONOS_ROOT_DIR)/cronos/InterfaceFunctions/solver.C
endif

ifneq ($(findstring user_pde.C, $(wildcard ./*.C)), )
APPLICATION_CXX += $(realpath ./user_pde.C)
else
APPLICATION_CXX += $(CRONOS_ROOT_DIR)/cronos/InterfaceFunctions/user_pde.C
endif

# Convert to absolute paths
# This allows the user to use relative paths as well (in the setup Makefile)
APPLICATION_CXX_REAL = $(foreach cxx, $(APPLICATION_CXX), $(realpath $(cxx)))

PROJ_CXX += $(APPLICATION_CXX_REAL)

# Create a joined list of all relevant source files (sort to remove duplicates)
# This is relevant (for the future) if more programs (that share parts of the cronos source files)
# should be compiled with this Makefile.
ALL_CXX = $(sort $(PROJ_CXX))


# --------------------------------------------------------------------------------------
#	OBJECT FILES AND DEPENDENCIES
# --------------------------------------------------------------------------------------

# Definition of helper functions

# Divide the source files in two groups (cronos/application)
# This is needed since the build paths will be different for the respective .o files
SORT_CXX	= $(filter $(CRONOS_ROOT_DIR)/%.C, $(1)) $(filter-out $(CRONOS_ROOT_DIR)/%.C, $(1))

# Generate a list of object files (.o) from a given list of source (.C) files.
CXX2OBJ		= $(patsubst $(CRONOS_ROOT_DIR)/%.C, $(BUILD_DIR)/%.o, $(filter $(CRONOS_ROOT_DIR)/%.C, $(1))) \
				$(patsubst %.C, $(BUILD_DIR)/application/%.o, $(notdir $(filter-out $(CRONOS_ROOT_DIR)/%.C, $(1))))

# Generate a list of mpi object files (_mpi.o) from a given list of normal object (.o) files.
OBJ2MPIOBJ		= $(patsubst %.o, %_mpi.o, $1)

# Generate a list of dependency files (.d) from a given list of normal object (.o) files.
OBJ2DEP			= $(patsubst %.o, %.d, $1)

# end function definition

# Generate the actual files (objects and dependencies)
ALL_CXX_SORTED 	= $(call SORT_CXX, $(ALL_CXX))
ALL_OBJS 		= $(call CXX2OBJ, $(ALL_CXX_SORTED))
ALL_MPI_OBJS	= $(call OBJ2MPIOBJ, $(ALL_OBJS))

ALL_DEPS		= $(call OBJ2DEP, $(ALL_OBJS))
ALL_MPI_DEPS	= $(call OBJ2DEP, $(ALL_MPI_OBJS))

# Include all .d files
-include $(ALL_DEPS)
-include $(ALL_MPI_DEPS)

# Generate object sublists needed for the high-level recipes
PROJ_OBJS			= $(call CXX2OBJ, $(PROJ_CXX))
PROJ_MPI_OBJS		= $(call OBJ2MPIOBJ, $(PROJ_OBJS))


# --------------------------------------------------------------------------------------
#	HIGH LEVEL RECIPES
# --------------------------------------------------------------------------------------

.DEFAULT_GOAL := all

.PHONY: all
all: proj proj_MPI

.PHONY: proj
proj: $(BUILD_DIR)/proj

.PHONY: proj_MPI
proj_MPI: $(BUILD_DIR)/proj_MPI

# Actual target of the binary - depends on all .o files.
$(BUILD_DIR)/proj: $(PROJ_OBJS) $(BUILD_DIR)/compiler_flags
	$(CXX) $(PROJ_OBJS) $(LDFLAGS) $(LOADLIBES) $(FFTLIBES) -o $@

$(BUILD_DIR)/proj_MPI: $(PROJ_MPI_OBJS) $(BUILD_DIR)/compiler_flags
	$(MPICXX) $(PROJ_MPI_OBJS) $(LDFLAGS) $(LOADLIBES) -D parallel -o $@

# Cleaning up
.PHONY: clean
clean:
	rm -rf $(BUILD_DIR)


# --------------------------------------------------------------------------------------
#	OBJECT RECIPIES
# --------------------------------------------------------------------------------------

# Build target for every single object file.
# The potential dependency on header files is covered
# by calling `-include $(ALL_DEPS)`.
# The -MMD flags additionaly creates a .d file with
# the same name as the .o file.

# Define the build rule template for the serial and mpi case
define SERIAL_BUILD_RULE_TEMPLATE
$(1): $(2) $$(BUILD_DIR)/compiler_flags
	@mkdir -p $$(@D)
	$$(CXX) $$(CXXFLAGS_FINAL) -MMD -c $$< -o $$@

endef
define MPI_BUILD_RULE_TEMPLATE
$(1): $(2) $$(BUILD_DIR)/compiler_flags
	@mkdir -p $$(@D)
	$$(MPICXX) $$(CXXFLAGS_FINAL) -D parallel -MMD -c $$< -o $$@

endef

# Now we need to deduce build rules for all object files
# To do this we need a joined iteration over both the object file and its respective cxx source file.
# Since there are no joined or indexed iterations we need to be creative
# Hack: create a joined list formated as (obj):(cxx) and iterate over that list
SERIAL_LIST_PRE = $(join $(ALL_OBJS:%=%:),     $(ALL_CXX_SORTED))
MPI_LIST_PRE    = $(join $(ALL_MPI_OBJS:%=%:), $(ALL_CXX_SORTED))

# Remove the buildinfo.o since we will define a manual recipe for this object later
SERIAL_LIST = $(filter-out %buildinfo.C, $(SERIAL_LIST_PRE))
MPI_LIST = $(filter-out %buildinfo.C, $(MPI_LIST_PRE))


# Define some helper functions to get the individual entries back
GET_OBJ = $(word 1,$(subst :, ,$1))
GET_CXX = $(word 2,$(subst :, ,$1))

# Finally apply the templates
$(foreach joined, $(SERIAL_LIST), $(eval $(call SERIAL_BUILD_RULE_TEMPLATE, $(call GET_OBJ, $(joined)), $(call GET_CXX, $(joined)))))
$(foreach joined, $(MPI_LIST),    $(eval $(call MPI_BUILD_RULE_TEMPLATE,    $(call GET_OBJ, $(joined)), $(call GET_CXX, $(joined)))))


# --------------------------------------------------------------------------------------
#	COMPILER FLAGS CHECK
# --------------------------------------------------------------------------------------

# Write the compiler flags to a file
$(shell mkdir -p $(BUILD_DIR) && echo '$(CXXFLAGS_RECOMPILE)' | cmp -s - $(BUILD_DIR)/compiler_flags || echo '$(CXXFLAGS_RECOMPILE)' > $(BUILD_DIR)/compiler_flags)

# Define a target that will be outdated if the compiler flags have changed
# This target will be used to trigger recompilation in this case
$(BUILD_DIR)/compiler_flags:

# Write build infos to a file
$(shell mkdir -p $(BUILD_DIR) && echo '$(GLOBAL_INFOS)' | cmp -s - $(BUILD_DIR)/build_infos || echo '$(GLOBAL_INFOS)' > $(BUILD_DIR)/build_infos)

# Define a target that will be outdated if the build infos have changed
# This target will be used to trigger recompilation in this case
$(BUILD_DIR)/build_infos:

# If the build infos have changed we only want to recompile the buildinfo getters (in buildinfo.C)
# We therefore override the previously deduced build recipe (in the previous section) for this source file

$(BUILD_DIR)/cronos/generic/buildinfo.o: $(CRONOS_ROOT_DIR)/cronos/generic/buildinfo.C $(BUILD_DIR)/compiler_flags $(BUILD_DIR)/build_infos
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS_FINAL) $(GLOBAL_INFOS) -MMD -c $< -o $@

$(BUILD_DIR)/cronos/generic/buildinfo_mpi.o: $(CRONOS_ROOT_DIR)/cronos/generic/buildinfo.C $(BUILD_DIR)/compiler_flags $(BUILD_DIR)/build_infos
	@mkdir -p $(@D)
	$(MPICXX) $(CXXFLAGS_FINAL) $(GLOBAL_INFOS) -D parallel -MMD -c $< -o $@


# --------------------------------------------------------------------------------------
#	LEGACY SYNTAX
# --------------------------------------------------------------------------------------

# Enable the legacy syntax to compile individual .o files
# Define rule template for the legacy syntax
define LEGACY_SYNTAX_RULE_TEMPLATE
.PHONY: $(notdir $(1))
$(notdir $(1)): $(1)
endef

# Evaluate the template for every object
$(foreach obj,$(PROJ_OBJS),$(eval $(call LEGACY_SYNTAX_RULE_TEMPLATE,$(obj))))
$(foreach obj,$(PROJ_MPI_OBJS),$(eval $(call LEGACY_SYNTAX_RULE_TEMPLATE,$(obj))))


# Disable implicit rules so canonical targets will work.
.SUFFIXES:
# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =
