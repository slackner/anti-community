CFLAGS = -O
CC     = gcc
LIBS   = -lm -lgsl -lgslcblas

.PHONY: all
all: estimate communities

estimate: estimate.o readgml.o
	$(CC) $(CFLAGS) -o estimate estimate.o readgml.o $(LIBS)

estimate.o: estimate.c network.h readgml.h Makefile
	$(CC) $(CFLAGS) -c estimate.c

communities: communities.o readgml.o
	$(CC) $(CFLAGS) -o communities communities.o readgml.o $(LIBS)

communities.o: communities.c network.h readgml.h Makefile
	$(CC) $(CFLAGS) -c communities.c

readgml.o: readgml.c readgml.h network.h
	$(CC) $(CFLAGS) -c readgml.c

.PHONY: clean
clean:
	rm -f *.o estimate communities
