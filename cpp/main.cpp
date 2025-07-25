#include <iostream>
#include <cstdlib>
#include "square.h"

// C++ entry point
int main(int argc, char* argv[])
{
	// Check that we received one argument
	if (argc != 2) {
		std::cout<<"USAGE: "<<argv[0]<<" number\n";
		return 1;
	}

	int n = atoi(argv[1]);
	int value = square(n);
	std::cout<<"The square of "<<n<<" is "<<value<<"\n";
}
