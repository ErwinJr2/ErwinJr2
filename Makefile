# CC = icc
# CFLAGS = -march=native -Ofast -gcc-name=gcc-6 -Wall -c -fPIC
CC = gcc
CFLAGS = -Ofast -Wall -c -fPIC

.PHONY : all
.DEFAULT : all
.SECONDARY : 

all: 1DSchrodinger.so 1DThermal.so 1DSelfConsistant.so

%.so : %.o
	$(CC) -shared -fPIC $< -o $@ -lm

%.o : %.c science.h
	# $(CC) $(CFLAGS) $< -o $@
	$(CC) $(CFLAGS) -D __DEBUG $< -o $@

.PHONY : clean
clean :
	@$(RM) {1DSchrodinger,1DThermal}.{so,o}
