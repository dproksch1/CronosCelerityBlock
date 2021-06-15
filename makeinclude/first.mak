include $(HOME)/src/makeinclude/ostype.mak

.SUFFIXES:

%:
	@test -d $(X_OSTYPE) || mkdir $(X_OSTYPE)
	@cd $(X_OSTYPE); $(MAKE) -f../Makefile $@

all:

