CXXFLAGS=-Werror -Wall  -I../../PsimagLite/src -I../../PsimagLite -g3 \
		 -I ../../LanczosPlusPlus/src -I../../LanczosPlusPlus/src/Engine
LDFLAGS=../../LanczosPlusPlus/src/ProgramGlobals.o -L../../PsimagLite/lib -lpsimaglite -llapack -lblas 

all: schwinger

schwinger: SchwingerExactDiag.o
	g++ -o schwinger SchwingerExactDiag.o $(LDFLAGS)

SchwingerExactDiag.o: SchwingerExactDiag.cpp
	g++ -c SchwingerExactDiag.cpp $(CXXFLAGS)

clean:
	rm -f SchwingerExactDiag.o schwinger
