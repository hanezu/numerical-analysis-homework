
.PHONY: all clean
.PRECIOUS: %.c %.o

all: poisson-neumann-region.out

%.out: %.o
	$(CC) -Wall -O3 -lm -o $@ $<

%.o: %.c
	$(CC) -Wall -O3 -c -o $@ $<

clean:
	$(RM) *.o *.out

