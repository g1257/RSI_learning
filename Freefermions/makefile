# Compiler and flags
CXX = g++
CXXFLAGS = -Wall -O2 -std=c++17

# Target executable
TARGET = Cdagger

# Source and object files
SRC = Cdagger.cpp
OBJ = Cdagger.o

# Libraries
LIBS = -llapack -lblas
STATIC_LIB = -lpsimaglite

# Include and lib paths (modify if needed)
INCLUDES = 
LIBDIRS = -L /home/georgebabus/codes/PsimagLite/lib
 
# Default target
all: $(TARGET)

# Linking step
$(TARGET): $(OBJ) 
	$(CXX) $(OBJ) $(STATIC_LIB) $(LIBDIRS) $(LIBS) -o $(TARGET)

# Compilation step
$(OBJ): $(SRC)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $(SRC) -o $(OBJ)


# Clean target
clean:
	rm -f $(OBJ) $(TARGET)

