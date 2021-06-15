include /home/dproksch/Documents/master/master_thesis/repos/CronosNumLib/makeinclude/ostype.mak

.SUFFIXES:

%:
	@test -d $(X_OSTYPE) || mkdir $(X_OSTYPE)
	@cd $(X_OSTYPE); $(MAKE) -f../Makefile $@

all:

