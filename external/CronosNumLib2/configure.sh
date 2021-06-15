#!/bin/bash

LOCAL_PATH=${PWD} 

echo ${LOCAL_PATH}
NUMLIB_ROOT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo ${NUMLIB_ROOT_DIR}

UTIL_DIR="$NUMLIB_ROOT_DIR/util"
UTIL_MAKE_LOW="$UTIL_DIR/makefile"
UTIL_MAKE_UP="$UTIL_DIR/Makefile"

echo ${UTIL_DIR}

MAKEINCLUDE_DIR="$NUMLIB_ROOT_DIR/makeinclude"

MAKE_LOW_STRING="include $MAKEINCLUDE_DIR/first.mak"
MAKE_UP_STRING="include $MAKEINCLUDE_DIR/include.mak"

# util dir:
echo " Treating $UTIL_DIR "

# write to makefile (with backup)
if [ ! -f "$UTIL_MAKE_LOW.back" ]; then
    cp $UTIL_MAKE_LOW "$UTIL_MAKE_LOW.back"
fi
echo "$MAKE_LOW_STRING" > $UTIL_MAKE_LOW
chmod u+x $UTIL_MAKE_LOW

# write to Makefile (with backup)
if [ ! -f "$UTIL_MAKE_UP.back" ]; then
    cp $UTIL_MAKE_UP "$UTIL_MAKE_UP.back"
fi

BUFF="$UTIL_MAKE_UP.dump"
sed "1 c $MAKE_UP_STRING" <$UTIL_MAKE_UP >$BUFF
mv $BUFF $UTIL_MAKE_UP
chmod u+x $UTIL_MAKE_UP



MATRIX_DIR="$NUMLIB_ROOT_DIR/Matrix"
MATRIX_MAKE_LOW="$MATRIX_DIR/makefile"
MATRIX_MAKE_UP="$MATRIX_DIR/Makefile"



# Matrix dir:
echo " Treating $MATRIX_DIR "

# write to makefile (with backup)
if [ ! -f "$MATRIX_MAKE_LOW.back" ]; then
    cp $MATRIX_MAKE_LOW "$MATRIX_MAKE_LOW.back"
fi
echo "$MAKE_LOW_STRING" > $MATRIX_MAKE_LOW
chmod u+x $MATRIX_MAKE_LOW

# write to Makefile (with backup)
if [ ! -d "$MATRIX_MAKE_UP.back" ]; then
    cp $MATRIX_MAKE_UP "$MATRIX_MAKE_UP.back"
fi

BUFF="$MATRIX_MAKE_UP.dump"
sed "1 c $MAKE_UP_STRING" <$MATRIX_MAKE_UP >$BUFF
mv $BUFF $MATRIX_MAKE_UP
chmod u+x $MATRIX_MAKE_UP


# Now address the makeinclude dir:
FIRST_MAK="$MAKEINCLUDE_DIR/first.mak"
INCLUDE_MAK="$MAKEINCLUDE_DIR/include.mak"

echo " Treating $MAKEINCLUDE_DIR "

# Make a copy and  edit file
OSTYPE_STRING="include $MAKEINCLUDE_DIR/ostype.mak"
OSTYPE_STRING2="include $MAKEINCLUDE_DIR/\$(X_OSTYPE)"

# treat first.mak
if [ ! -f "$FIRST_MAK.back" ]; then
    cp $FIRST_MAK "$FIRST_MAK.back"
fi

BUFF="$FIRST_MAK.dump"
sed "1 c $OSTYPE_STRING" <$FIRST_MAK >$BUFF
mv $BUFF $FIRST_MAK
chmod u+x $FIRST_MAK


#treat include.mak
if [ ! -f "$INCLUDE_MAK.back" ]; then
    cp $INCLUDE_MAK "$INCLUDE_MAK.back"
fi

NUMLIB_ROOT_DIR_STRING="X_LIB_ROOT_DIR	= $NUMLIB_ROOT_DIR"

BUFF="$INCLUDE_MAK.dump"
#sed "s:$(HOME):$basedir:" <$include_mak >$buff
sed -e "1 c $OSTYPE_STRING" -e "2 c $OSTYPE_STRING2" \
    -e "6 c $NUMLIB_ROOT_DIR_STRING" <$INCLUDE_MAK >$BUFF
mv $BUFF $INCLUDE_MAK
chmod u+x $INCLUDE_MAK


# Get OS type to make corresponding directories (if not available)
X_OSTYPE=$(uname)-$(uname -m | sed -e s/x86_64/amd64/ -e s/.86/386/  -e 's/\/800//')
echo $X_OSTYPE

echo "Making directories"
NUMLIB_LIB_DIR="$NUMLIB_ROOT_DIR/lib"
NUMLIB_LIB_SUBDIR="$NUMLIB_ROOT_DIR/lib/$X_OSTYPE"
NUMLIB_INC_DIR="$NUMLIB_ROOT_DIR/include"

if [ ! -d "$NUMLIB_LIB_DIR" ]; then
	echo "making lib dir"
	mkdir $NUMLIB_LIB_DIR
fi

if [ ! -d "$NUMLIB_LIB_SUBDIR" ]; then
	echo "making lib sub dir"
	mkdir $NUMLIB_LIB_SUBDIR
fi

if [ ! -d "$NUMLIB_INC_DIR" ]; then
	echo "making include dir"
	mkdir $NUMLIB_INC_DIR
fi

export NUMLIB_ROOT_DIR="$NUMLIB_ROOT_DIR"

echo " Finished "

