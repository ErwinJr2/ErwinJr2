CFLAGS += -Ofast -Wall -fPIC #-D _DEBUG
# CC = icc
# CFLAGS = -march=native -Ofast -gcc-name=gcc-6 -Wall

ifeq ($(OS),Windows_NT)
    CFLAGS += -D WIN32
    EXT := dll
else 
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Darwin)
        EXT := dylib
    else
        # Assuming linux
        EXT := so
    endif
endif

.PHONY : all
.DEFAULT : all
.SECONDARY : 


all: 1DSchrodinger.$(EXT) 1DThermal.$(EXT) 1DMaxwell.$(EXT)
MP:  1DSchrodinger_MP.$(EXT)

1DSchrodinger_MP.$(EXT) : 1DSchrodinger_MP.o band.o
	$(CC) -shared -fPIC -fopenmp $^ -o $@ -lm

1DSchrodinger.$(EXT) : 1DSchrodinger.o band.o
	$(CC) -shared -fPIC $^ -o $@ -lm

%.$(EXT) : %.o
	$(CC) -shared -fPIC $^ -o $@ -lm

1DSchrodinger.o : 1DSchrodinger.c science.h band.h
	$(CC) $(CFLAGS) -c $< -o $@

1DSchrodinger_MP.o : 1DSchrodinger.c science.h band.h
	$(CC) $(CFLAGS) -fopenmp -D __MP -c $< -o $@

%.o : %.c science.h 
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY : clean
clean :
	@$(RM) {1DSchrodinger,1DThermal,1DMaxwell,1DSchrodinger_MP}.{so,o,dll,dylib} band.o
