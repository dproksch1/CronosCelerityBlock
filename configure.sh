#!/bin/bash

usage ()
{
    echo "NAME"
    echo "      $(basename "$0") - Script to configure Cronos "
    echo ""
    echo "SYNOPSIS"
    echo "      $(basename "$0") [CronosNumLib_Dir]"
    echo ""
    echo "DESCRIPTION"
    echo "      This script has to be executed before CronosCode can be compiled."
    echo "      It sets up the relevant paths in ./cronos/make/path_include,"
    echo "      i.e. linking the CronosNumLib."
    echo ""
    echo "OPTIONS"
    echo "      [CronosNumLib_Dir]"
    echo "              Path to the ConosNumLib root directory."
    echo "              This has to be set if neither the environment variable"
    echo "              NUMLIB_ROOT_DIR nor NUMLIB_PATH are set."
}

do_setup ()
{
    # Get rid of tildes
    MY_NUMLIB_PATH="${MY_NUMLIB_PATH/#\~/$HOME}"

    # convert path if relative
    MY_NUMLIB_PATH=$(realpath $MY_NUMLIB_PATH)

    # exit if path is not a directory
    if ! [ -d "$MY_NUMLIB_PATH" ]; then
        echo "Invalid Path: $MY_NUMLIB_PATH"
        echo ""
        usage
        exit 1
    fi

    echo "Using path: $MY_NUMLIB_PATH"

    CRONOS_ROOT_DIR=$(dirname $0)
    CRONOS_ROOT_DIR=$(realpath $CRONOS_ROOT_DIR)
    PATH_INCLUDE_FILE=$CRONOS_ROOT_DIR/cronos/make/path_include
    echo "Setting up: $PATH_INCLUDE_FILE"

    echo "X_LIB_ROOT_DIR = $MY_NUMLIB_PATH" > $PATH_INCLUDE_FILE
    echo "CRONOS_ROOT_DIR = $CRONOS_ROOT_DIR" >> $PATH_INCLUDE_FILE
}

# check if the legacy syntax is used
for i do
    if [[ $i == --prefix* ]]; then
        MY_NUMLIB_PATH=${i#*=}
        echo " -- Using legacy syntax -- "
        do_setup
        exit 0
    fi
done

# Use argument if specified
if [ $# -eq 1 ]
then
    MY_NUMLIB_PATH=$1
    do_setup
    exit 0
fi

# Use environment variables if set
if ! [ -z "$NUMLIB_ROOT_DIR" ] ; then
    echo "Using environment variable: NUMLIB_ROOT_DIR"
    MY_NUMLIB_PATH="$NUMLIB_ROOT_DIR"
    do_setup
    exit 0
fi

if ! [ -z "$NUMLIB_PATH" ] ; then
    echo "Using environment variable: NUMLIB_PATH"
    MY_NUMLIB_PATH="$NUMLIB_PATH"
    do_setup
    exit 0
fi

# if we did nothing so far
usage
exit 1



