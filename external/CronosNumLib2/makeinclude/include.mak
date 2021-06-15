include /home/dproksch/Documents/master/master_thesis/repos/CronosNumLib/makeinclude/ostype.mak
include /home/dproksch/Documents/master/master_thesis/repos/CronosNumLib/makeinclude/$(X_OSTYPE)

VPATH		= ..

X_LIB_ROOT_DIR	= /home/dproksch/Documents/master/master_thesis/repos/CronosNumLib
X_INC_DIR	= $(X_LIB_ROOT_DIR)/include
X_LIB_DIR	= $(X_LIB_ROOT_DIR)/lib/$(X_OSTYPE)
X_INC		= -I$(X_INC_DIR)
X_LIB		= -L$(X_LIB_DIR)

MOD_DATA	= 644
MOD_EXEC	= 755

%.a : 
	$(AR) $(ARFLAGS) $@ $^

LINK.o = $(CXX) $(LDFLAGS) $(TARGET_ARCH)
