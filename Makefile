# Compile third party libraries

# set targets
TARGETS = QPBO

all: $(TARGETS)
.PHONY: $(TARGETS)

# Recursive Make
QPBO:
	$(MAKE) -C thirdparty/QPBO-v1.32.src

