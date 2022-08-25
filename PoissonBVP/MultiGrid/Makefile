RM = rm -rf

default:
	$(MAKE) -C src

.PHONY:run
run:
	./bin/main1D > ./output/OneDimGrid.txt
	./bin/main2D > ./output/TwoDimGrid.txt
	./bin/display > ./output/display.txt

.PHONY:test
test:
	$(MAKE) -C test

.PHONY:clear
clear:
	$(RM) ./bin/main1D ./bin/main2D ./bin/display
	$(RM) ./output/OneDimGrid.txt ./output/TwoDimGrid.txt ./output/display.txt