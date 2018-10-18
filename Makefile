# CC = icc
# CFLAGS = -march=native -Ofast -gcc-name=gcc-6 -Wall -c -fPIC
CC = gcc
CFLAGS = -Ofast -Wall -c -fPIC

.PHONY : all
.DEFAULT : all
.SECONDARY : 

all: 1DSchrodinger.so 1DThermal.so 1DMaxwell.so

1DSchrodinger_MP.so : 1DSchrodinger_MP.o
	$(CC) -shared -fPIC -fopenmp $< -o $@ -lm

%.so : %.o
	$(CC) -shared -fPIC $< -o $@ -lm

1DSchrodinger_MP.o : 1DSchrodinger.c science.h
	# $(CC) $(CFLAGS) -fopenmp -D __MP $< -o $@
	$(CC) $(CFLAGS) -D __DEBUG -fopenmp -D __MP $< -o $@

%.o : %.c science.h
	# $(CC) $(CFLAGS) $< -o $@
	$(CC) $(CFLAGS) -D __DEBUG $< -o $@

.PHONY : clean
clean :
	@$(RM) {1DSchrodinger,1DThermal,1DMaxwell}.{so,o}
