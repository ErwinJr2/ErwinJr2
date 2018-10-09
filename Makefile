# CC = icc
# CFLAGS = -march=native -Ofast -gcc-name=gcc-6 -Wall -c -fPIC
CC = gcc
CFLAGS = -Ofast -Wall -c -fPIC

.PHONY : all
.DEFAULT : all
all: 1DSchrodinger.so

1DSchrodinger.so : 1DSchrodinger.o
	$(CC) -shared -fPIC $< -o $@ -lm

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY : clean
clean :
	$(RM) 1DSchrodinger.so 1DSchrodinger.o
