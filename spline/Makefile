CC = g++
RM = rm -rf
COEF = -llapacke -llapack
INCL = -Iinclude

# default target
Exercise:
	$(CC) -o ./bin/exercise_B ./src/exercise_B.cpp $(INCL) $(COEF)
	$(CC) -o ./bin/exercise_CD ./src/exercise_CD.cpp $(INCL) $(COEF)
	$(CC) -o ./bin/exercise_E ./src/exercise_E.cpp $(INCL) $(COEF)
	$(CC) -o ./bin/exercise_F ./src/exercise_F.cpp $(INCL) $(COEF)

.PHONY: test
test:
	$(CC) -o ./bin/test ./test/test.cpp $(INCL) $(COEF)

.PHONY: latex
latex:
	$(MAKE) -C doc

.PHONY: clean 
clean:
	$(RM) ./bin/exercise_B ./bin/exercise_CD ./bin/exercise_E ./bin/exercise_F ./bin/test
	$(RM) ./output/*.txt