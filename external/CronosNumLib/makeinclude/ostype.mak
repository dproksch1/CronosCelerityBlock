X_OSTYPE := $(shell uname)-$(shell uname -m | sed -e s/x86_64/amd64/ -e s/.86/386/  -e 's/\/800//')
