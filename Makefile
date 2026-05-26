CC = cc
CFLAGS = -O2 -Wall -std=c11
OBJ = matrices.o metrics.o solvers.o timing.o experiments.o main.o

lab: $(OBJ)
	$(CC) $(CFLAGS) -o lab $(OBJ) -lm

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f $(OBJ) lab

.PHONY: clean
