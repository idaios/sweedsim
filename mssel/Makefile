CC=gcc
CFLAGS=-O3 -g
DEPS=mssel.h
LIBS=-lgsl -lgslcblas -lm

%.o: %.c  $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: mssel

mssel: rand1.o streecsel.o mssel.o kde.o
	gcc -o mssel rand1.o streecsel.o mssel.o kde.o $(LIBS)
	cp mssel ~/bin
clean:
	rm -f *.o mssel





