CXXFLAGS=-Werror -Wall -g3 -I../../PsimagLite/src
LDFLAGS=-L../../PsimagLite/lib -lpsimaglite -llapack -lblas

all: schwinger

schwinger: SchwingerExactDiag.o
	g++ -o schwinger SchwingerExactDiag.o $(LDFLAGS)

SchwingerExactDiag.o: SchwingerExactDiag.cpp
	g++ -c SchwingerExactDiag.cpp $(CXXFLAGS)

