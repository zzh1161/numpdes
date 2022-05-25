CC = g++
RM = rm -rf
AUTO = -fconcepts-ts
LIB = -llapacke -llapack

default:
	$(MAKE) -C src

.PHONY: test
test:
	$(MAKE) -C test

.PHONY: run
run:
	./bin/main > ./output/convergence.txt
	./bin/plot > ./output/plot.txt

.PHONY: clean
clean:
	$(RM) ./bin/main ./bin/plot ./bin/test