SUBDIRS = $(dir $(wildcard */Makefile))

.PHONY: all
all: $(SUBDIRS)

.PHONY: $(SUBDIRS)
$(SUBDIRS):
	$(MAKE) -C $@

.PHONY: clean
clean:
	for dir in $(SUBDIRS); do \
		rm -f "$$dir/"*.txt; \
		rm -f "$$dir/"*.gml; \
	done
