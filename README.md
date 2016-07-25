# ParetoFrontEnumerationAlgorithm
This repository contains an implementation for an algorithm to enumerate all elements of a Pareto front for a multi-criterial optimization problem for which all optimization objectives have a finite range.

The algorithm is explained in the paper "Enumerating the Complete Pareto Front" that can be found [here](http://arxiv.org/abs/1512.05207).

Implementations for C++ and Python are available. Both of them are licensed under the relatively permissive MIT license, so that they can be used in other projects. The algorithm is not wrapped into a formal library as this would be overkill. Rather, it is suggested that the respective files are simply copied to projects that make use of the algorithm.

If you make use of the algorithm in your work, please cite the paper stated above. This will help to spread the word about it.


## C++ version
A C++ version of the algorithm can be found in the "c++" directory. The files "pareto_enumerator.hpp" and "pareto_enumerator.cpp" are the ones that need to be included in other projects in order to use the algorithm. A simple test program can be compiled (under Linux) by running

> g++ -g -std=c++14 -Wall -Wextra -O3 tests.cpp pareto_enumerator.cpp -o tests

in the "c++" directory. Alternative, it should be possible to compile the test program using an IDE, provided that C++14 support ist turned on. The algorithm itself (consisting of the files "pareto_enumerator.hpp" and "pareto_enumerator.cpp") should be usable under C++11.

The usage of the algorithm implementation is straight-forward. An easy-to-read example is given in function "doSimpleTest" of "tester.cpp"

A seed for the random number generator can be provided at the command line. In case you report a bug that can be found with the tester program, please include the seed value (that the tester program prints to the console) along with your bug report.

There is still a bit of room for optimizing the C++ implementation, as there is a tradeoff between readability and speed. In order to make porting the code to other languages easier, the focus for the implementation was more on the readability side.


## Python version
A Python version of the algorithm can be found in the "python" directory. It was tested both under Python 3.4.3 and Python 2.7.10. If the module is executed as main module, it performs some tests using the "unittest" module. The usage is pretty straight-forward and the tests in the module show how it can be used.
