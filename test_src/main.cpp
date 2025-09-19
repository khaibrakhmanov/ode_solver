#include <iostream>
#include <fstream>
#include <cmath>
#include "test1.h"
#include "test2.h"
#include "test3.h"
#include "test4.h"
#include "test5.h"

int main(int argc, char* argv[])
{
	std::cout << "Press Enter to start a simulation" << std::endl;
	std::cin.get();

	std::cout << "Running test 1 ..." << std::endl;
	CTest1 test1;
	test1.Run();
	std::cout << "... Done!" << std::endl;

	std::cout << "Running test 2 ..." << std::endl;
	CTest2 test2;
	test2.Run();
	std::cout << "... Done!" << std::endl;

	std::cout << "Running test 3 ..." << std::endl;
	CTest3 test3;
	test3.Run();
	std::cout << "... Done!" << std::endl;

	std::cout << "Running test 4 ..." << std::endl;
	CTest4 test4;
	test4.Run();
	std::cout << "... Done!" << std::endl;

	std::cout << "Running test 5 ..." << std::endl;
	CTest5 test5;
	test5.Run();
	std::cout << "... Done!" << std::endl;

	return 0;
}
