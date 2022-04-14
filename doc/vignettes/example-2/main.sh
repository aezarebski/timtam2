#!/usr/bin/env sh

CHAIN_LENGTH=100000
BEAST_JAR=../../../../beast2/build/dist/beast.jar
BEAST_CLASS=beast.app.beastapp.BeastMain

EXTRA_ARGS="-overwrite -seed 1 -D chainLength=$CHAIN_LENGTH"
BEAST_CMD="java -cp $BEAST_JAR $BEAST_CLASS $EXTRA_ARGS"

$BEAST_CMD xml/bdsky-serial.xml
mv ex2-bdsky-serial.log out/

$BEAST_CMD xml/timtam.xml
mv ex2-timtam.log out/

Rscript src/make-posterior-plots.R
