#!/bin/bash

## Make a distribution of the DADA2 nifH pipeline.
## Creates a temp directory, copies needed files there, and tar's them.
## This script must be consistent with INSTALL.txt:
##   -- The path to the distributed tgz's
##   -- The extracted files must be within directory "DADA2_nifH_pipeline"
##      because the path setting step assumes this.
##
## Some of the scripts needed by the DADA2 pipeline are general utilities that
## will be used in other pipelines, e.g. extractFasta.pl. When making the tgz,
## copy over such utilities to the scripts dir.  Installation of the DADA2
## pipeline sets the path to include scripts.
##

STAMP=$(echo $(date) | awk -F' ' '{print $6 $2 $3}')
OUTTGZ="Shared/DADA2_nifH_pipeline.$STAMP.tgz"
if [ -f "$OUTTGZ" ] ; then
    read -p "$OUTTGZ already exists. Want to clobber it? " ans
    if [ "$ans" != 'y' ] && [ "$ans" != 'Y' ] ; then
        echo "Aborted."
        exit -1
    fi
fi

UDIR="$HOME/ToolsShared/Utils"
UTILS="$UDIR/extractFasta.pl"
NEEDED_STUFF="run_DADA2_pipeline.sh Installation/ Example/ scripts/"

## The files in the tgz must be below "DADA2_nifH_pipeline"
echo "Tarring up $OUTTGZ ..."
mkdir -p Shared/DADA2_nifH_pipeline/scripts && \
  cp -Rf `echo $UTILS` Shared/DADA2_nifH_pipeline/scripts/ && \
  cp -Rf `echo "$NEEDED_STUFF"` Shared/DADA2_nifH_pipeline/ && \
  tar cfz "$OUTTGZ" -C Shared DADA2_nifH_pipeline
if [ "$?" -ne 0 ] ; then
    echo "Something went wrong. Cleaning up..."
    rm -f "$OUTTGZ"  # silent if doesn't exist
fi
## Clean up even if tar fails.
rm -rf Shared/DADA2_nifH_pipeline
echo "Done!"
