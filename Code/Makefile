# Makefile for Writing Make Files Example
 
# *****************************************************
# Variables to control Makefile operation
 
CC = g++
CFLAGS = -Wall -g
 
# ****************************************************
# Targets needed to bring the executable up to date
 
CF2: addAWorm.o check.o defs.o growTheWorm.o initialise.o main.o measure.o next.o nspFunctions.o qfFunctions.o viewer.o
	$(CC) $(CFLAGS) -o CF2 main.o addAWorm.o check.o defs.o growTheWorm.o initialise.o measure.o next.o nspFunctions.o qfFunctions.o viewer.o
 
# The main.o target can be written more simply
 
main.o: addAWorm.o defs.o initialise.o measure.o qfFunctions.o viewer.h
	$(CC) $(CFLAGS) -c main.cpp
 
addAWorm.o: addAWorm.h defs.h next.h check.h growTheWorm.h

check.o: check.h

defs.o: defs.h

growTheWorm.o: defs.h growTheWorm.h next.h

initialise.o: defs.h initialise.h nspFunctions.h

measure.o: defs.h measure.h next.h nspFunctions.h

next.o: defs.h next.h

nspFunctions.o: defs.h nspFunctions.h

qfFunctions.o: defs.h qfFunctions.h addAWorm.h measure.h

viewer.o: defs.h viewer.h
