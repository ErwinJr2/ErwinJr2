CFLAGS += -Ofast -Wall -fPIC
# CFLAGS += -D _DEBUG
# CC = icc
# CFLAGS = -march=native -Ofast -gcc-name=gcc-6 -Wall
OMPFLAGS = -fopenmp

ifeq ($(OS),Windows_NT)
	CFLAGS += -D WIN32
	EXT := dll
else
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Darwin)
		CFLAGS += -D MACOS
		EXT := dylib
		OMPFLAGS += -I/usr/local/opt/llvm/includ
		LDFLAGS += -L/usr/local/opt/llvm/lib
	else
		# Assuming linux
		EXT := so
	endif
endif

.PHONY : all
.DEFAULT : all
.SECONDARY :


all: 1DSchrodinger.$(EXT)
MP:  1DSchrodinger_MP.$(EXT)

1DSchrodinger_MP.$(EXT) : 1DSchrodinger_MP.o band.o fftautocorr.o
	$(CC) -shared -fPIC $(OMPFLAGS) $^ -o $@ -lm $(LDFLAGS)

1DSchrodinger.$(EXT) : 1DSchrodinger.o band.o fftautocorr.o
	$(CC) -shared -fPIC $^ -o $@ -lm

1DSchrodinger.o : 1DSchrodinger.c science.h band.h fftautocorr/fftautocorr.h
	$(CC) $(CFLAGS) -c $< -o $@

1DSchrodinger_MP.o : 1DSchrodinger.c science.h band.h fftautocorr/fftautocorr.h
	$(CC) $(CFLAGS) -D __MP $(OMPFLAGS) -c $< -o $@

band.o : band.c band.h
	$(CC) $(CFLAGS) -c $< -o $@

fftautocorr.o : fftautocorr/fftautocorr.c fftautocorr/fftautocorr.h fftautocorr/factortable.h
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY : clean
clean :
	@$(RM) {1DSchrodinger,1DSchrodinger_MP}.{so,o,dll,dylib} band.o fftautocorr.o

