# This file is part of fftautocorr
# 
# Copyright (C) 2020 CareF
# Author: CareF
# Licensed under a 3-clause BSD style license - see LICENSE.md
#
# The library is intended to be used at source code level. FFTW and pocketfft
# dependencies is introduced only for testing purposes.
#
# The Makefile is presented as a demo of compilation and as infrastructure
# for testing purposes. Compilation settings may vary depending on the platform.

CFLAGS += -O3 -Wall

TESTFLAGS += -lfftw3
ifeq ($(OS),Windows_NT)
	echo Win not tested!
else 
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Darwin)
		TESTFLAGS += -L/usr/local/lib -I/usr/local/include 
	else
		# Assuming linux
		LDLIBS += -lm
	endif
endif


.PHONY : test clean
.DEFAULT_GOAL := fftautocorr.o
.SECONDARY : fftautocorr.o

test: test.out
	./test.out

test.out : test.c fftautocorr.o  pocketfft.o
	$(CC) $(CFLAGS) $^ $(LDLIBS) $(TESTFLAGS) $(LDFLAGS) -o $@

fftautocorr.o : fftautocorr.c fftautocorr.h factortable.h
	$(CC) $(CFLAGS) -c $< -o $@

pocketfft.o : pocketfft/pocketfft.c pocketfft/pocketfft.h
	$(CC) $(CFLAGS) -c $< -o $@

factortable.h : tablegen.py
	./tablegen.py $@

pocketfft/% : 
	git clone https://gitlab.mpcdf.mpg.de/mtr/pocketfft.git

.PHONY : clean
clean :
	@$(RM) *.o test.out
