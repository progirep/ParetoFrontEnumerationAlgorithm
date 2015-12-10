#include "pareto_enumerator.hpp"
#include <iostream>
#include <sstream>
#include <random>
#include <set>

/*
 * This is
 *   tests.cpp
 * that is part of the ParetoFrontEnumerationAlgorithm library, available from
 *   https://github.com/progirep/ParetoFrontEnumerationAlgorithm
 *
 * Is ia a library for enumerating all elements of a Pareto front for a
 * multi-criterial optimization problem for which all optimization objectives
 * have a finite range.
 *
 * The library and all of its files are distributed under the following license:
 *
 * -----------------------------------------------------------------------------
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2015 Ruediger Ehlers
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

//=================================================================================
// Helper functions for the tests
//=================================================================================
bool vectorOfIntIsSmaller(const std::vector<int>& a, const std::vector<int> &b) {
    bool diff = false;
    for (unsigned int i = 0; i<a.size();i++) {
        if (b.at(i)<a.at(i)) return false;
        diff |= b.at(i)!=a.at(i);
    }
    return diff;
}

bool vectorOfIntIsLeq(const std::vector<int>& a, const std::vector<int> &b) {
    for (unsigned int i = 0; i<a.size();i++) {
        if (b.at(i)<a.at(i)) return false;
    }
    return true;
}

std::list<std::vector<int> > cleanParetoFront(const std::list<std::vector<int> > &input) {
    std::list<std::vector<int> > cleanedElements;
    for (auto &it : input) {
        bool foundSmaller = false;
        for (auto &it2 : input) {
            foundSmaller |= (vectorOfIntIsSmaller(it2,it));
        }
        if (!foundSmaller) cleanedElements.push_back(it);
    }
    return cleanedElements;
}


//=================================================================================
// First test: Find a fixed Pareto front
//=================================================================================
bool simpleObjectiveFunction(const std::vector<int> &point) {
    if (point.size()!=3) throw "Function simpleObjectiveFunction was called with wrong number of elements.p";
    if (point.at(0)>5) return true;
    if (point.at(1)<3) return false;
    return (point.at(2)>7);
}

void doSimpleTest() {
    std::vector<std::pair<int,int> > limits;
    limits.push_back(std::pair<int,int>(0,10));
    limits.push_back(std::pair<int,int>(0,10));
    limits.push_back(std::pair<int,int>(0,10));
    std::list<std::vector<int> > front = paretoenumerator::enumerateParetoFront(simpleObjectiveFunction,limits);

    // Check Pareto optima
    if (front.size()>2) throw "Error: Found too many Pareto points in function doSimpleTest";
    bool foundA = false;
    bool foundB = false;
    for (auto a : front) {
        if ((a[0]==6) && (a[1]==0) && (a[2]==0)) {
            foundA = true;
        } else if ((a[0]==0) && (a[1]==3) && (a[2]==8)) {
            foundB = true;
        }
    }
    if ((!foundA) || (!foundB)) throw "Error: Not all Pareto points have been found in function doSimpleTest (plus possibly 1-2 too many points)";

}

//=================================================================================
// Second test: Find randomly generated Pareto fronts
//              -> Also check that no redundant calls to the feasibility
//                 function are made!
//=================================================================================
bool randomTestFeasibilityFunction(const std::vector<int> &point, const std::list<std::vector<int> > &paretoPoints, std::list<std::vector<int> > &positiveBuffer, std::list<std::vector<int> > &negativeBuffer ) {

    // Check if implied by the previous positive results
    for (auto &a : positiveBuffer) {
        if (vectorOfIntIsLeq(a,point)) throw "Error: Called the feasibility function on a point that is already known to map to TRUE.";
    }

    // Check if implied by the previous negative results
    for (auto &a : negativeBuffer) {
        if (vectorOfIntIsLeq(point,a)) throw "Error: Called the feasibility function on a point that is already known to map to FALSE.";
    }

    for (auto &a : paretoPoints) {
        if (vectorOfIntIsLeq(a,point)) {
            positiveBuffer.push_back(point);
            return true;
        }
    }
    negativeBuffer.push_back(point);
    return false;
}

void doRandomTest(unsigned int randomSeed) {
    std::mt19937 rand(randomSeed);
    std::vector<std::pair<int,int> > limits;
    std::list<std::vector<int> > paretoPoints;

    // Randomize the number of dimensions
    unsigned int nofDimensions = rand() % 7 + 5;

    // Randomize the number of points (prior to the removal of dominated solutions
    unsigned int nofPoints = rand() % 15 + 1;

    // Randomize the limits
    for (unsigned int i=0;i<nofDimensions;i++) {
        int min = rand() % 100 - 50;
        int max = rand() % 100 + min + 1;
        limits.push_back(std::pair<int,int>(min,max));
    }

    // Randomize the points
    for (unsigned int i=0;i<nofPoints;i++) {
        std::vector<int> newPoint;
        for (unsigned int j=0;j<nofDimensions;j++) {
            newPoint.push_back(rand() % (limits[j].second - limits[j].first) + limits[j].first);
        }
        paretoPoints.push_back(newPoint);
    }

    // Clean the points
    paretoPoints = cleanParetoFront(paretoPoints);

    // Enumerate points. The feasibility function
    // stores all results in the "positive" and "negative" Buffer
    // and checks for each of its calls whether the calling algorithm could have deduced
    // the results from results to earlier calls
    std::list<std::vector<int> > positiveBuffer;
    std::list<std::vector<int> > negativeBuffer;

    std::function<bool(const std::vector<int> &)> fun = [paretoPoints,&positiveBuffer,&negativeBuffer] (const std::vector<int> &point) { return randomTestFeasibilityFunction(point,paretoPoints,positiveBuffer, negativeBuffer); };
    std::list<std::vector<int> > front = paretoenumerator::enumerateParetoFront(fun,limits);
    std::set<std::vector<int> > frontSet(front.begin(),front.end());

    // Compare the correct Pareto front and the actual results
    for (auto a : paretoPoints) {
        if ((frontSet.count(a))==0) throw "Element not found in Pareto set";
    }
    if (frontSet.size()!=paretoPoints.size()) throw "Too many points found";
}

//=================================================================================
// Main function
//=================================================================================
int main(int nofArgs, const char **args) {
    try {
        unsigned int randomSeed = std::random_device{}();
        if (nofArgs>1) {
            std::istringstream is(args[1]);
            is >> randomSeed;
            if (is.fail()) {
                std::cerr << "Error: Expect at most one parameter: the random seed to random testing, which must be an integer." << std::endl;
                return 1;
            }
        }
        std::cout << "Random Seed: " << randomSeed << "\nTest progress: ";

        // Do simple test
        doSimpleTest();

        // Do random test
        for (unsigned int i=0;i<1000;i++) {
            if ((i % 20)==0) {
                std::cout << ".";
                std::cout.flush();
            }
            doRandomTest(randomSeed+i);
        }
        std::cout << "\nAll tests finished correctly.\n";
        return 0;
    } catch (const char *errorMsg) {
        std::cerr << errorMsg << std::endl;
        return 1;
    }
}
