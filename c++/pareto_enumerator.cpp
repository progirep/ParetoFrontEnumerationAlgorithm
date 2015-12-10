#include "pareto_enumerator.hpp"

/*
 * This is
 *   pareto_enumerator.cpp
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

namespace paretoenumerator {

    inline bool vectorOfIntIsSmaller(const std::vector<int>& a, const std::vector<int> &b) {
        const size_t size = a.size();
        for (size_t i = 0; i<size;i++) {
            if (b[i]<a[i]) return false;
            if (b[i]!=a[i]) {
                // We continue the outer loop now here as this is
                // faster than storing whether a smaller
                // element has been found in a flag.
                for (i++;i<size;i++) {
                    if (b[i]<a[i]) return false;
                }
                return true;
            }
        }
        return false;
    }

    bool vectorOfIntIsLeq(const std::vector<int>& a, const std::vector<int> &b) {
        const size_t size = a.size();
        for (size_t i = 0; i<size;i++) {
            if (b[i]<a[i]) return false;
        }
        return true;
    }

    /**
     * @brief Removes all dominating elements from a set of search space points
     * @param input The initial set of points
     * @return The cleaned set of points
     */
    std::list<std::vector<int> > cleanParetoFront(const std::list<std::vector<int> > &input) {
        std::list<std::vector<int> > cleanedElements;
        for (auto &it : input) {
            bool foundSmaller = false;
            for (auto &it2 : input) {
                if (vectorOfIntIsSmaller(it,it2)) {
                    foundSmaller = true;
                    break;
                }
            }
            if (!foundSmaller) cleanedElements.push_back(it);
        }
        return cleanedElements;
    }


    /**
     * @brief A class that buffers negative results from the feasibility function so that no
     * redundant calls are made to it.
     *
     * Dominated points are removed from the buffer
     */
    class NegativeResultBuffer {
        std::list<std::vector<int> > oldValueBuffer;
    public:
        bool isContained(const std::vector<int> &data) {
            for (auto &a : oldValueBuffer) {
                if (vectorOfIntIsLeq(data,a)) return true;
            }
            return false;
        }

        void addPoint(const std::vector<int> &data) {
            for (auto it = oldValueBuffer.begin();it!=oldValueBuffer.end();) {
                if (vectorOfIntIsLeq(*it,data)) {
                    it = oldValueBuffer.erase(it);
                } else {
                    it++;
                }
            }
            oldValueBuffer.push_back(data);
        }
    };


    /**
     * @brief Main function of the pareto front element enumeration algorithm
     * @param fn the feasibility function
     * @param limits the upper and lower bounds of the objective values. In every pair, the minimal value comes first.
     * @return the list of Pareto points.
     */
    std::list<std::vector<int> > enumerateParetoFront(std::function<bool(const std::vector<int> &)> fn, std::vector<std::pair<int,int> > limits) {

        // Buffer the number of dimensions of the search space
        unsigned int nofDimensions = limits.size();

        // Reserve the sets "P" and "S" from the paper
        std::list<std::vector<int> > paretoFront;
        std::list<std::vector<int> > coParetoElements;

        // Negative result buffer
        NegativeResultBuffer negativeResultBuffer;

        // Add the maximal element to the coParetoElements
        std::vector<int> maximalElement;
        for (auto i : limits) {
            maximalElement.push_back(i.second);
        }
        coParetoElements.push_back(maximalElement);

        // Main loop
        while (coParetoElements.size()>0) {
            std::vector<int> &testPoint = coParetoElements.front();
            if (!(negativeResultBuffer.isContained(testPoint))) {
                if (fn(testPoint)) {
                    // A Pareto point is missing. Let us find where exactly it is.
                    // We need to work on a copy of the point in order not to spoil
                    // the point form the coParetoElements
                    std::vector<int> x = testPoint;
                    for (unsigned int i=0;i<nofDimensions;i++) {
                        int max = x[i]+1;
                        int min = limits[i].first;
                        while ((max - min)>1) {
                            int mid = min + ((max-min-1)/2);
                            x[i] = mid;
                            if (negativeResultBuffer.isContained(x)) {
                                min = mid+1;
                            } else {
                                if (fn(x)) {
                                    max = mid+1;
                                } else {
                                    min = mid+1;
                                    negativeResultBuffer.addPoint(x);
                                }
                            }
                        }
                        x[i] = min;
                    }
                    paretoFront.push_back(x);

                    // Now update all points in the coParetoFront
                    std::list<std::vector<int> > coParetoElementsMod;
                    for (auto &y : coParetoElements) {
                        if (!vectorOfIntIsLeq(x,y)) {
                            coParetoElementsMod.push_back(y);
                        } else {
                            for (unsigned int i=0;i<nofDimensions;i++) {
                                if (x[i]>limits[i].first) {
                                    coParetoElementsMod.push_back(y);
                                    std::vector<int> &mod = coParetoElementsMod.back();
                                    mod[i] = x[i]-1;
                                }
                            }
                        }
                    }
                    coParetoElements = cleanParetoFront(coParetoElementsMod);

                } else {
                    // Get rid of this point in the co-Pareto front and add to the negative results buffer
                    negativeResultBuffer.addPoint(testPoint);
                    coParetoElements.pop_front();
                }
            } else {
                // Get rid of this point in the co-Pareto front
                coParetoElements.pop_front();
            }
        }
        return paretoFront;
    }

} // End of namespace
