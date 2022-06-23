RM = rm -rf

all: test plot

.PHONY: test plot story run clean
test:
	$(MAKE) -C test

plot:
	$(MAKE) -C src

story:
	$(MAKE) -C report

run:
	./bin/test
	./bin/plot

clean:
	$(RM) ./bin/test ./bin/plot ./output/*.txt