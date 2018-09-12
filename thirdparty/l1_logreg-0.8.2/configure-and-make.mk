# Make target for a recursive make.
.PHONY: doit

doit:
	./configure && $(MAKE)

clean:
	$(MAKE) clean

