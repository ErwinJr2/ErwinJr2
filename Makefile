# CC = icc
# CFLAGS = -march=native -Ofast -gcc-name=gcc-6 -Wall -fPIC
CC = gcc
CFLAGS = -Ofast -Wall -fPIC

.PHONY : all
.DEFAULT : all
.SECONDARY : 

all: 1DSchrodinger.so 1DThermal.so 1DMaxwell.so

1DSchrodinger_MP.so : 1DSchrodinger_MP.o band.o
	$(CC) -shared -fPIC -fopenmp $^ -o $@ -lm
	mv 1DSchrodinger_MP.so 1DSchrodinger.so

1DSchrodinger.so : 1DSchrodinger.o band.o
	$(CC) -shared -fPIC $^ -o $@ -lm

%.so : %.o
	$(CC) -shared -fPIC $^ -o $@ -lm

1DSchrodinger.o : 1DSchrodinger.c science.h band.h
	$(CC) $(CFLAGS) -c $< -o $@

1DSchrodinger_MP.o : 1DSchrodinger.c science.h band.h
	$(CC) $(CFLAGS) -fopenmp -D __MP -c $< -o $@

%.o : %.c science.h 
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY : clean
clean :
	@$(RM) {1DSchrodinger,1DThermal,1DMaxwell}.{so,o} band.o 1DSchrodinger_MP.o
