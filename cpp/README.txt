Title: Difference between compilation and linking

See square.h and square.cpp included in this write-up.

Do the following compilation:

g++ -c square.cpp

This does compilation only because of the -c flag passed to the compiler named g++.

The compiled object is square.o
This file cannot be run because of various reasons.
First it doesn't have an "entry point", or starting point to start execution.
Second, it may not have required functions or libraries for basic things like printing.

If we want to execute something that uses the function, we need an entry point.
This is always int main().
So we write another file (main.cpp) where we call our function.

See main.cpp
Note that main.cpp includes only the declaration square.h but not its implementation.

Do the following compilation:

g++ -c main.cpp

which creates main.o which is an object but not executable.
main.o has an entry point (the function int main()) but does not
have required libraries or functions. It doesn't have square.o for example.

Do the following linking:

g++ -o main   main.o square.o

which creates the file main, which is executable.
This stage does linking not compilation. The compiler (in this case g++) understands
that it needs to link because of the .o given have already been compiled.
The compiler does some work but for the actual linking g++ calls another computer program
called (appropriately) the linker, ld. You may do man g++ and man ld to learn more about them.

--------------------------------------------------------------------------

More details about linking vs. compilation.

We could have done compilation and linking in one go. Do the following cleaning:
rm main.o square.o main

And then:
g++ main.cpp square.cpp -o main

which creates main. This command does first compilation and then linking.

We see then two approaches:
(1) compile in stages, and then link all stages
(2) compile and link in the same stage without intermediate compilations.

We use (2) if we want to divide compilation in parts to speed up compilation.
We can compile main.cpp and square.cpp in parallel at the same time, for example.

We use (1) to achieve better optimization. The compiler can optimize better
if it has all source files passed to it at the same time.
These days there's work being done for link-time-optimization or LTO.

We know what main.cpp, square.cpp and square.h contain, because we either wrote them ourselves,
or we can read the C++ code they contain.

What do main.o and square.o contain?
They are binary files that contain compiled code that the CPU can understand.

Let's look at square.o with the tool nm that lists the functions that it has:

Execute the following:

nm square.o

I get the following:

0000000000000000 T _Z6squarei

It tells me that it has only one function called _Z6squarei which is defined (letter T means defined in the text code section of the binary)
You may do
man nm
for more info on the nm tool

C++ allows function name overloading. But in machine code (assembly), functions are just memory locations
so they need a unique name. Because of this, the C++ function name (square) is mangled: the compiler uses extra info from the function,
like the fact that it takes an integer, to create squarei and mangles it even further to _Z6squarei
This is just a name for a memory address, which here appears as tons of zeros (first column above).

Do the same with the object main.o:
nm main.o

The first three lines I get are the following:

                 U atoi
0000000000000000 T main
                 U _Z6squarei

There is only one defined function: main, which is marked with the letter T.
Note that there are undefined functions marked with U: these functions are needed but not defined in main.o.
We see our _Z6squarei which is undefined in main.o, and that is one reason why main.o cannot execute: it needs linking with square.o among other things.
We also see atoi which is part of the standard library that gets linked in.

---------------------

Summary: Compilation and linking are different processes done by the compiler and the linker. Linking is always necessary to obtain an executable.


