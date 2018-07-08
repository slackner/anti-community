LDFLAGS  = -Wl,--no-undefined -Wl,-rpath,'$$ORIGIN' -L. -I. -lgraph

.PHONY: all
all: libgraph.so

libgraph.so: $(wildcard *.c *.h)
	gcc -O2 -g -Werror -Wall -fstack-protector-strong -fPIC \
		-shared -Wl,-soname,libgraph.so -Wl,--no-undefined \
		-o libgraph.so $^

tests: tests.c graph.h internal.h libgraph.so
	gcc -O2 -g -Werror -Wall -fstack-protector-strong -o $@ $< $(LDFLAGS)

.PHONY: test
test: tests
	./tests

.PHONY: clean
clean:
	rm -f libgraph.so tests
