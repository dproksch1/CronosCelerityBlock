include $(HOME)/src/makeinclude/ostype.mak
include $(HOME)/src/makeinclude/$(X_OSTYPE)

VPATH		= ..

X_ROOT_DIR	= $(HOME)
X_INC_DIR	= $(X_ROOT_DIR)/include
X_LIB_DIR	= $(X_ROOT_DIR)/lib/$(X_OSTYPE)
X_INC		= -I$(X_ROOT_DIR)/include -I/usr/lib/openmpi/include
X_LIB		= -L$(X_ROOT_DIR)/lib/$(X_OSTYPE)

MOD_DATA	= 644
MOD_EXEC	= 755

%.a : 
	$(AR) $(ARFLAGS) $@ $^

LINK.o = $(CXX) $(LDFLAGS) $(TARGET_ARCH)
