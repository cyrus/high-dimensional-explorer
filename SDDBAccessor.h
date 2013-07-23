/* 

Copyright (C) 2004,2005,2006,2007,2008,2009,2010,2011,2012,2013  Cyrus Shaoul and Geoff Hollis 

This file is part of HiDEx.

    HiDEx is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HiDEx is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HiDEx in the COPYING.txt file.
    If not, see <http://www.gnu.org/licenses/>.


*/

//
// SEMANTIC DISTANCE DATABASE ACCESSOR
// 
// This class allows you to create a new chunk of the
// semantic distance database. It creates a file that
// will hold N Semantic Distance vectors, each with a
// set length, and independant-length behind-and-ahead
// window sizes.
//
#ifndef SDDBACCESSOR
#define SDDBACCESSOR

#include <cassert>
#include <string>
#include <utility>
#include <map>
#include <vector>
#include "Matrix.h"

#define MAX_FILENAME_LENGTH          (200)
#define RECORD_EXTENTION         ".record"

//Frequency sorter is for sorting words by frequency

typedef pair<size_t,int> FreqEntry;
//typedef pair<int,size_t> SizeEntry;
typedef vector<FreqEntry> FreqSorter;
//typedef vector<SizeEntry> SizeData;
//typedef map<int, int> FrequencyMap;

// Constructor
//    AccessorData(int numVectors, int vectorLen,  int windowLenBehind, int windowLenAhead);

class FreqSort3
{
public:
    bool operator() (const FreqEntry & a, const FreqEntry & b) const
    {
        return a.first < b.first;
    }
};

class SDDBAccessor {
  //
  // A class that stores all of the variables/private functions
  // we might possibly use. 
  //

 public:
  //
  // Create a new SDDBAccessor, with the given parameters
  //
  SDDBAccessor(const char *accessorName, const int numVectors,
	       const int vectorLen, 
	       const int windowLenBehind, 
	       const int windowLenAhead, const size_t maxMemory);

  //
  // Close/delete the SDDBAccessor
  // 
  ~SDDBAccessor();


  //
  // Increment the co-occurance value between vectorNum and vectorPos
  // at distance, by amount
  //
  void add(const int vectorNum, const int vectorPos,
	   const int distance, const int amount);


  //
  // add everything about the passed matrix to our current info
  // about our matrix at vectorNum
  //
  void add(const int vectorNum, const Matrix<int> *matrix);


  //
  // Set the co-occurance value between vectorNum and vectorPos
  // at distance, to amount
  //
  void set(const int vectorNum, const int vectorPos,
	   const int distance, int amount);


  //
  // Get the co-occurance value between vectorNum and vectorPos at distance,
  // by amount, but don't do anything to it.
  //
  int get(const int vectorNum, const int vectorPos,
		     const int distance);


  void show(const int vectorNum) const;

  //
  // Returns a pointer to a two-dimensional array for the vector of vectorNum,
  // in the specified window range.
  //
  Matrix<int> *getMatrix(const int vectorNum, 
                         const int windowLenBehind,
                         const int windowLenAhead);
  

  //
  // Load the entry into RAM for quick read/writing
  // return true if load was successful
  // return false if vector is already in RAM or load
  // was unsuccessful
  //
  void loadToRAM(const int vectorNum);
  
  //
  // returns true of the vector is in RAM, and false
  // otherwise
  //
  bool vectorInRAM(const int vectorNum);

  // Flush the Cache if it is too big.
  void FlushCache();

  //
  // These functions return a copy of the corresponding private variables
  //
  int myWindowLenBehind();
  int myWindowLenAhead();
  int myVectorLen();
  int myNumVectors();

  //
  // Flush anything we might have waiting
  //
  void flush();


  //
  // Flush what info we have about the specified vector to disk.
  // remove the vector from RAM afterwards
  //
  void flush(const int vectorNum);


  //
  // Flush the entry and remove it from RAM
  // return true if the removal was successful
  // return false if the vector was not in RAM, or the
  // removal was unsuccessful
  //
  bool removeFromRAM(const int vectorNum);

  int getMatrixElement(const int vectorNum) const
      {
          if (_vectors[vectorNum] == NULL) 
              return 0;
          else { 
              for (int i = 0; i < 10; i++) {
                  for (int j = 0; i < 10; i++) {
                      cerr << _vectors[vectorNum]->get(i,j) << " " ;
                  }
              }
              return 1;
          }
      };

  //
  // Close the filestream down
  // 
  void close();

 private:

  string _name;
        
  int _numVectors;       // referenced as an array: 0, 1, ..., _numVectors-1
        
  int _vectorLen;        // referenced as an array: 0, 1, ..., _vectorLen-1

  int _windowLenBehind;  // referenced as neg distance: -1, -2, ..., _windowLenBehind
        
  int _windowLenAhead;   // referenced as actual distance: 1, 2, ..., _windowLenAhead
  
  size_t _MemorySize;   // Amount of memory used by all matricies in object

  size_t _maxMemory;   // maximum amount of memory to use.

  vector<size_t> _freq; //  Frequency of Matrix usage.  

  FreqSorter _MatrixCache; // Sorted Freq of matrix usage.

  vector<size_t> _MatrixSizes; // Sorted Freq of matrix usage.

  Matrix<int> **_vectors;    // the entries we have put into memory

  int distToPos(int distance) const;

  bool entryInRAM(const int entry);

  bool NeedToDumpData();

  void DumpLowestFreqData();

  
};

#endif // SDDBACCESSOR
