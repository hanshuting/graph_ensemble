# Compile third party libraries

# exclude blossom for now
#TARGETS = QPBO blossom csapp JustinsGraphicalModelsToolboxPublic l1graph l1_logreg
TARGETS = QPBO JustinsGraphicalModelsToolboxPublic l1graph l1_logreg fast_gibbs_sampler pmtk3

all: $(TARGETS)
.PHONY: $(TARGETS)

# Recursive Make
QPBO:
	$(MAKE) -C thirdparty/QPBO-v1.32.src

blossom:
	$(MAKE) -C thirdparty/blossom

csapp:
	$(MAKE) -C thirdparty/csa++

JustinsGraphicalModelsToolboxPublic:
	cd thirdparty/JustinsGraphicalModelsToolboxPublic && matlab -nojvm -r "compile; exit"

l1graph:
	$(MAKE) -C thirdparty/l1graph_clean

l1_logreg:
	$(MAKE) -C thirdparty/l1_logreg-0.8.2 -f configure-and-make.mk

#fast_gibbs_sampler:
#	$(MAKE) -C src/structure_learning/gibbs_sampler

pmtk3:
	cd thirdparty/pmtk3-master && matlab -nojvm -r "initPmtk3; exit"

