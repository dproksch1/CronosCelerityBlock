include /home/dproksch/master_thesis/CronosCelerityBlock/external/CronosNumLib/makeinclude/ostype.mak

.SUFFIXES:

%:
	@test -d $(X_OSTYPE) || mkdir $(X_OSTYPE)
	@cd $(X_OSTYPE); $(MAKE) -f../Makefile $@

all:

