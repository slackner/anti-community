LDFLAGS  = -Wl,--no-undefined -Wl,-rpath,'$$ORIGIN/lib' -L./lib -I./lib -lgraph
SOURCES  = $(wildcard *.c)
PROGRAMS = $(patsubst %.c,%,$(SOURCES)) minmodularity antimodularity

.PHONY: all
all: $(PROGRAMS)

lib/libgraph.so: FORCE
	$(MAKE) -C lib libgraph.so

lib/test: FORCE
	$(MAKE) -C lib test

%: %.c lib/graph.h lib/libgraph.so
	gcc -O2 -g -Werror -Wall -fstack-protector-strong -o $@ $< $(LDFLAGS)

minmodularity: maxmodularity.c lib/graph.h lib/libgraph.so
	gcc -O2 -g -Werror -Wall -fstack-protector-strong -o $@ -DNEWMAN_FLAGS=NEWMAN_FLAGS_MINIMIZE  $< $(LDFLAGS)

antimodularity: maxmodularity.c lib/graph.h lib/libgraph.so
	gcc -O2 -g -Werror -Wall -fstack-protector-strong -o $@ -DNEWMAN_FLAGS=NEWMAN_FLAGS_ANTIMODULARITY  $< $(LDFLAGS)

FORCE:

.PHONY: test
test: lib/test

.PHONY: clean
clean:
	$(MAKE) -C lib clean
	rm -f $(PROGRAMS)
