SUBDIRS = $(dir $(wildcard */Makefile))

.PHONY: all
all: $(SUBDIRS)

datasets: src third-party

.PHONY: $(SUBDIRS)
$(SUBDIRS):
	$(MAKE) -C $@
