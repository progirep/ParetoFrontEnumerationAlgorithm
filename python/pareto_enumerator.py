#!/usr/bin/env python3
"""
This is
   pareto_enumerator.py
that is part of the ParetoFrontEnumerationAlgorithm library, available from
   https://github.com/progirep/ParetoFrontEnumerationAlgorithm

Is ia a library for enumerating all elements of a Pareto front for a
multi-criterial optimization problem for which all optimization objectives
have a finite range.

The library and all of its files are distributed under the following license:

-----------------------------------------------------------------------------

The MIT License (MIT)

Copyright (c) 2015 Ruediger Ehlers

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import sys,curses, random, unittest, itertools

#=========================================================
# Tools
#=========================================================
def pointwiseLeq(a,b):
    assert len(a)==len(b)
    for i in range(0,len(a)):
        if a[i] > b[i]:
            return False
    return True
    
def cleanParetoSet(pointSet):
    ret = []
    for a in pointSet:
        found = False
        for b in pointSet:
            if pointwiseLeq(b,a) and (a!=b):
                found = True
        if not found:
            ret.append(a)
    return frozenset(ret)
    
def isParetoSet(pointSet):
    for a in pointSet:
        for b in pointSet:
            if pointwiseLeq(b,a) and (a!=b):
                return False
    return True

def cleanCoParetoSet(pointSet):
    ret = []
    for a in pointSet:
        found = False
        for b in pointSet:
            if pointwiseLeq(a,b) and (a!=b):
                found = True
        if not found:
            ret.append(a)
    return frozenset(ret)

def complement(dimensions,pointSet):
    pointSet = cleanParetoSet(pointSet)
    complement = [tuple([MAX for i in range(0,dimensions)])]
    for a in pointSet:
        oldComplement = complement
        complement = []
        for e in oldComplement:        
            foundSmaller = False
            for i in range(0,dimensions):
                if e[i]<a[i]:
                    foundSmaller = True
            if foundSmaller:
                complement.append(e)
            else:
                for i in range(0,dimensions):
                    if a[i]>0:
                        complement.append(e[0:i]+(a[i]-1,)+e[i+1:])
        complement = cleanCoParetoSet(complement)
    return complement                    
                
                
#=========================================================
# Buffer class for negative answers
#=========================================================
class NegativeAnswerBuffer:
    
    def __init__(self):
        self.buffer = set([])
        
    def isContained(self,testPoint):
        for a in self.buffer:
            if pointwiseLeq(testPoint,a):
                return True
        return False
        
    def add(self,testPoint):
        self.buffer.add(testPoint)
        
#=========================================================
# Pareto enumeration function
#=========================================================
def computeParetoFront(fn,bounds):
    paretoSet = []
    restSet = [tuple([b for (a,b) in bounds])]
    negativePoints = NegativeAnswerBuffer()
    while len(restSet)>0:
        thisOne = restSet[0]
        if negativePoints.isContained(thisOne):
            restSet = restSet[1:]
        else:
            if (not negativePoints.isContained(thisOne)) and fn(thisOne):
                # Binary search
                for i in range(0,len(bounds)):
                    maxi = thisOne[i]+1
                    mini = 0
                    while (maxi-mini)>1:
                        mid = mini+(maxi-mini-1)//2
                        testPoint = thisOne[0:i]+(mid,)+thisOne[i+1:]
                        found = False
                        if negativePoints.isContained(testPoint):
                            mini = mid+1
                        else:
                            if fn(testPoint):
                                maxi = mid+1
                            else:
                                negativePoints.add(testPoint)
                                mini = mid+1
                    thisOne = thisOne[0:i]+(mini,)+thisOne[i+1:]
                paretoSet.append(thisOne)

                # Change restSet
                oldRestSet = restSet
                restSet = []
                for e in oldRestSet:
                    foundSmaller = False
                    for i in range(0,len(bounds)):
                        if e[i]<thisOne[i]:
                            foundSmaller = True
                    if foundSmaller:
                        restSet.append(e)
                    else:
                        for i in range(0,len(bounds)):
                            if thisOne[i]>0:
                                restSet.append(e[0:i]+(thisOne[i]-1,)+e[i+1:])
                restSet = list(cleanCoParetoSet(restSet))
            else:
                negativePoints.add(thisOne)
                restSet = restSet[1:]
    return paretoSet


#======================================
# Testing
#======================================
class TestParetoRecurse(unittest.TestCase):

    #======================================
    # Example oracles
    #
    # All oracles need to be MONOTONE 
    # in all axes
    #======================================
    @staticmethod 
    def oracle2b(positions,values):
        assert len(values)==2
        contained = values[0]>=positions[0] and values[1]>=positions[1] or values[0]>=positions[2] and values[1]>=positions[3]
        sys.stdout.flush()
        return contained
       
    def testSimpleBinarySearch(self):
        # Test1 - Binary search on 0--20
        def oracle1(limit,values):
            assert len(values)==1
            return values[0]>=limit
            
        for i in range(0,20):
            paretoFront = computeParetoFront(lambda x,limit=i : oracle1(limit,x),[(0,20)])
            self.assertTrue(len(paretoFront)==1)
            self.assertTrue( (i,) in paretoFront)
    
    def testStairs(self):    
        # Test2 - Stairs
        def oracle2(limit,values):
            sys.stdout.flush()
            assert len(values)==2
            return values[0]+values[1] >= limit
            
        for i in range(0,40):
            paretoFront = computeParetoFront(lambda x,limit=i : oracle2(limit,x),[(0,20),(0,20)])
            self.assertTrue(len(paretoFront) == len(set(paretoFront)))
            nofElements = 0
            for j in range(0,21):
                data = i-j
                if data>=0 and data<=20:
                    assert (j,data) in paretoFront
                    nofElements += 1
            self.assertTrue(nofElements==len(paretoFront))

    def testTwoDPareto(self):
        # Test3 - Random Pareto elements
        for p1 in range(0,21):
            for p2 in range(0,21):
                for p3 in range(0,p1):    
                    for p4 in range(p2+1,21):    
                        paretoFront = computeParetoFront(lambda x,limit=(p1,p2,p3,p4) : self.oracle2b(limit,x),[(0,20),(0,20)])
                        self.assertTrue(len(paretoFront)==2)

    def testTinyWorkspaces(self):
    
        # Test 4: Tiny workspaces - X direction
        for p1 in range(0,21):
            for p2 in range(0,p1+2):
                paretoFront = computeParetoFront(lambda x,limit=(p2,0,p2,0) : self.oracle2b(limit,x),[(0,p1),(0,0)])
                if p2==p1+1:
                    self.assertTrue(len(paretoFront)==0)
                else:
                    self.assertTrue(len(paretoFront)==1)
                    self.assertTrue((p2,0) in paretoFront)

        # Test 5: Tiny workspaces - Y direction
        for p1 in range(0,21):
            for p2 in range(0,p1+2):
                paretoFront = computeParetoFront(lambda x,limit=(0,p2,0,p2) : self.oracle2b(limit,x),[(0,0),(0,p1)])
                if p2==p1+1:
                    self.assertTrue(len(paretoFront)==0)
                else:
                    self.assertTrue(len(paretoFront)==1)
                    self.assertTrue((0,p2) in paretoFront)


    def testThreeDimensions(self):
        def oracle3(limit,values):
            assert len(values)==3
            return values[0]+values[1]+values[2] >= limit
              
        # Test 6: Three dimensions        
        for i in range(0,60):
            paretoFront = computeParetoFront(lambda x,limit=i : oracle3(limit,x),[(0,20),(0,20),(0,20)])
            assert len(paretoFront) == len(set(paretoFront))
            nofElements = 0
            for x in range(0,21):
                for y in range(0,21):
                    data = i-x-y
                    if data>=0 and data<=20:
                        self.assertTrue((x,y,data) in paretoFront)
                        nofElements += 1
            self.assertTrue(nofElements==len(paretoFront))
            
    def testThreeDRandomParetos(self):
        # This test also tests that no dominated oracle calls are made.
        def oracle3b(limit,data,values):
            assert len(values)==3
            for (a1,a2,a3) in limit:
                if a1<=values[0] and a2<=values[1] and a3<=values[2]:
                    # Check that there is no known data value that implies this one
                    for b in data[0]:
                        if (b[0]<=values[0]) and (b[1]<=values[1]) and (b[2]<=values[2]):
                            assert False
                    data[0].append(tuple(values))
                    return True
            # Check that there is no known data value that implies this one
            for b in data[1]:
                if (b[0]>=values[0]) and (b[1]>=values[1]) and (b[2]>=values[2]):
                    assert False
            data[1].append(tuple(values))
            return False
        
        # Test 7: Three-D triple optima
        for i in range(0,100):
            # Coordinates for the test
            a = [(10,20,30),(20,10,30),(1,30,20)]
            # Randomize...
            a = [(a1+random.randrange(0,5),a2+random.randrange(0,5),a3+random.randrange(0,5)) for (a1,a2,a3) in a]
            # a = [(14, 23, 30), (21, 12, 34), (5, 32, 20)]
            # Prepare an array that stores all requests to the oracle
            oracleDataStorage = ([],[])        
            paretoFront = computeParetoFront(lambda x,limit=a,data=oracleDataStorage : oracle3b(limit,data,x),[(0,40),(0,40),(0,40)])
            self.assertTrue(len(paretoFront)==3)


# If called as main function, run Some tests
if __name__ == "__main__":
    unittest.main(verbosity=2,failfast=True,buffer=True)
    sys.exit(0)

