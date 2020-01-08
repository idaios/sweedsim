CC=gcc
CFLAGS=-O2
DEPS=freqints.h
LIBS=-lm

%.o: %.c  $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: trajdemognpops stepftn2

trajdemognpops: popsize.o binomial.o rand1.o trajdemognpops.o
	gcc -o trajdemognpops popsize.o binomial.o rand1.o trajdemognpops.o $(LIBS)
	cp trajdemognpops ~/bin

stepftn2: stepftn2.o
	gcc -o stepftn2 stepftn2.o $(LIBS)
	cp stepftn2 ~/bin

clean:
	rm -f *.o trajdemog stepftn2





