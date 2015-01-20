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

#define _LARGE_FILES
// ****************************************************************************
//
// These are the functions neccessary for a rather simple
// implementation of the SDDBAccessor with a hash-like
// style
//
//
// How does this class work?
//
// Each record in the accessor has its own file. However, there are many parts
// to the record - all of the same length - and not all neccessarily have to be
// saved (we know exactly what the "default" or "uninitialized" piece looks
// like, so there's no point in saving it when its easier and more 
// space-efficient to simply recontstruct it every time we need it). 
//
// So, instead, what we do is store each piece in a file, in order of which 
// piece was created first. Then, to figure out where a certain piece of a 
// record is located in the file, we see when it was created, with respect to 
// the other records.
//
// ****************************************************************************


#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "utilities.h"
#include "SDDBAccessor.h"
#include "SDDB.h"

using namespace std;

// ****************************************************************************
// Defines, constants, typedefs, and macros
// ****************************************************************************


// #define _numVectors              (_data->numVectors)
// #define _vectorLen               (_data->vectorLen)
// #define _windowLenBehind         (_data->windowLenBehind)
// #define _windowLenAhead          (_data->windowLenAhead)
// #define _vectors                 (_data->vectors)



// ****************************************************************************
// Structures and private data----Moving to Class SDDBData
// ****************************************************************************



// struct sddb_accessor_data
// {
//     string name;

//   int numVectors;       // referenced as an array: 0, 1, ..., _numVectors-1
//   int vectorLen;        // referenced as an array: 0, 1, ..., _vectorLen-1
//   int windowLenBehind;  // referenced as neg distance: -1, -2, ..., _windowLenBehind
//   int windowLenAhead;   // referenced as actual distance: 1, 2, ..., _windowLenAhead
//   Matrix<int> **vectors;    // the entries we have put into memory
// };

//****************************************************************************
// AccessorData functions
//****************************************************************************




// ****************************************************************************
// Local functions
// ****************************************************************************



// ****************************************************************************
// Public functions (SDDBAccessor API)
//
// Documentation in:
// SDDBAccessor.h
//
// ****************************************************************************
SDDBAccessor::SDDBAccessor(const char *accessorName, const int numVectors,
			   const int vectorLen, 
			   const int windowLenBehind, 
			   const int windowLenAhead) {

  //cerr << "SDDBAccessor::SDDBAccessor()" << endl;
    
    
    _numVectors      = numVectors;
    _vectorLen       = vectorLen;
    _windowLenAhead  = windowLenAhead;
    _windowLenBehind = windowLenBehind;
    _name = accessorName;
    //    _maxMemory = maxMemory;

    _vectors = new Matrix<int> *[_numVectors];

    _MatrixSizes.reserve(_numVectors);
    _freq.reserve(_numVectors);

    for(int i = 0; i < _numVectors; i++)
    {
        _vectors[i] = NULL;
        _freq[i] = 0;
        _MatrixSizes[i] = 0;
    }
}


SDDBAccessor::~SDDBAccessor() { 
  //  cerr << "SDDBAccessor::~SDDBAccessor()" << endl;
  for(int i = 0; i < _numVectors; i++) {

    if(!entryInRAM(i))
      continue;

    delete _vectors[i];
    _vectors[i] = NULL;
  }

  delete[] _vectors;
  //  delete _data;
}


void SDDBAccessor::add(const int vectorNum, const int vectorPos,
		       const int distance, const int amount) {
  //cerr << "SDDBAccessor::add()" << endl;
  loadToRAM(vectorNum);
  //    cerr << "SDDBAccessor::add() " << vectorNum << "\t" 
  //         << vectorPos << endl;
  _vectors[vectorNum]->add(vectorPos, distToPos(distance), amount);
}


void SDDBAccessor::add(const int vectorNum, const Matrix<int> *matrix)
{
    loadToRAM(vectorNum);
    _vectors[vectorNum]->add(matrix);
}

int SDDBAccessor::distToPos(int distance) const
{
    if (distance < 0)
        return _windowLenBehind + distance;
    else
        return _windowLenBehind + distance - 1;
}

void SDDBAccessor::set(const int vectorNum, const int vectorPos,
		       const int distance, int amount)
{
    //cerr << "SDDBAccessor::set()" << endl;
    loadToRAM(vectorNum);
    _vectors[vectorNum]->set(vectorPos, distToPos(distance), amount);
}


//int SDDBAccessor::get(const int vectorNum, const int vectorPos,const int distance) {
//  cerr << "SDDBAccessor::get() has not been completed." << endl;
//  return 0;
//}

void SDDBAccessor::show(const int vectorNum) const 
{
    if (_vectors[vectorNum] != NULL) {
        _vectors[vectorNum]->show(cout);
        cout << endl;
    }
}

bool SDDBAccessor::entryInRAM(const int entry) {
    if (_vectors[entry] != NULL) 
        return true;
    else
        return false;
}

// void SDDBAccessor::FlushCache() {
//     if (NeedToDumpData()) {
//         DumpLowestFreqData();
//     }
// }

// Manage Matrix Cache.
void SDDBAccessor::loadToRAM(const int vectorNum) {
    if(entryInRAM(vectorNum)) {
        //        cerr << "SDDBAccessor::loadToRAM() Vector " << vectorNum << " alread loaded to RAM" << endl;
        _freq[vectorNum]++;
    } else {
        //if we need to clear some data from cache, do it.
        // Check if data is on disk right now.
        ostringstream filename;
        filename << _name << vectorNum << RECORD_EXTENTION;
        //    cout << "Filname for db to open = "<< filename.str() << endl;
        Matrix<int>* M = Matrix<int>::fromFile(filename.str().c_str(), _vectorLen,
                                               _windowLenBehind + _windowLenAhead);
        // if not there, create an empty
        if(M == NULL) {
            _vectors[vectorNum] = new Matrix<int>(_vectorLen, _windowLenBehind+_windowLenAhead);
            _freq[vectorNum]++;
            //cerr << "SDDBAccessor::loadToRAM(), Created empty Matrix for " << vectorNum << endl;
        }
        else {
            _vectors[vectorNum] = M;
            _freq[vectorNum]++;
            //            cerr << "SDDBAccessor::loadToRAM(), loaded Matrix from disk " << vectorNum << endl;
        }
    }
}

// bool SDDBAccessor::NeedToDumpData() {
//     size_t matrixCount = 0;

//     _MemorySize = 0;
//     //cerr << "\n NeedToDump:: Getting Sizes:" << timestamp() << "\n";
//     for(int i = 0; i < _numVectors; i++) {
//         if (_vectors[i] != NULL) {
//             matrixCount++;
//             size_t MatrixSize =  _vectors[i]->size();
//             _MemorySize += MatrixSize;
//             _MatrixSizes[i] = MatrixSize;
//             //            cerr << ".";
//         }
//     }
//     //    cerr << "NeedToDump:: Got Memory Usage:" << timestamp() << "\n";
//     //    cerr << "NeedToDump:: No. of matricies in RAM =" << matrixCount << ", memory used= " << _MemorySize << "Max = " << _maxMemory <<endl;
//     if (_MemorySize > _maxMemory) {
//         //    if (static_cast<Float>(matrixCount) > (static_cast<Float>(_numVectors)*1)) {
//         cerr << "NeedToDump:: No. of matricies in RAM =" << matrixCount << ", memory used= " << _MemorySize << ", Max = " << _maxMemory <<endl;
//         return true;
//         }
//     else {
//         return false;
//     }
// }

// void SDDBAccessor::DumpLowestFreqData(){
    
//     cerr << "DumpLowest::Dumping excess data from Cache" << endl;
//     //    cerr << "DumpLowest::Sorting cache by usage frequency " << timestamp() <<  endl;
    
//     //Collect Matrix Sizes

//     // if memory is too big, drop low freq matricies
//     _MatrixCache.clear();
//     for(int i = 0; i < _numVectors; i++) {
//         if (_freq[i] != 0) {
// 	  _MatrixCache.push_back(FreqEntry(_freq[i],i));
//             //            sum += static_cast<Float>(i->second);
//         }
//     }
    
//     // sort matricies by frequency of usage.
//     sort(_MatrixCache.begin(), _MatrixCache.end(), FreqSort3());
//     // << "DumpLowest::Finished sort: " << timestamp() << endl;
    
//     size_t counter = 0;
//     size_t dumped = 0;
//     for (FreqSorter::iterator i = _MatrixCache.begin(); i !=  _MatrixCache.end(); i++) {
//         //        if (static_cast<Float>(i->first) < threshold) {
//         // subject matrix size from total memory usage
//         _MemorySize -= _MatrixSizes[i->second];
//         dumped +=  _MatrixSizes[i->second];
//         if (_MemorySize > _maxMemory) {
//             counter++;
//             flush(i->second);
//             //                cerr << "DumpLowest::Dropping Matrix " << i->second << " from cache with Freq = " << i->first  << " current MemorySize =" << MemorySize << " Max is: " << MAX_MEMORY<<"\n";
//         } else 
//             break;
//     }
//     Float percentage = static_cast<Float>(counter) / static_cast<Float>(_numVectors);
//     cerr << "DumpLowest:: Dropped " << counter << " vectors = " << dumped  << " bytes.\n DumpLowest:: I dumped ";
//     cerr.precision(10); 
//     cerr << percentage * 100.0;
//     cerr << " percent of the vectors " << endl;
// }



bool SDDBAccessor::vectorInRAM(const int vectorNum) {
  return entryInRAM(vectorNum);
}


Matrix<int> *SDDBAccessor::getMatrix(const int vectorNum, 
                                     int windowLenBehind,
                                     int windowLenAhead)
{
    //cerr << "SDDBAccessor::getMatrix()" << endl;
    if(entryInRAM(vectorNum)) {
        return matrixCopy<int, int>(_vectors[vectorNum], 
                                    0, _vectorLen - 1,
                                    _windowLenBehind - windowLenBehind,
                                    _windowLenBehind + windowLenAhead - 1);
        //    cout << "Found this matrix in RAM : "<< vectorNum << endl;
    }
    else
    {
        ostringstream filename;
        filename << _name << vectorNum << RECORD_EXTENTION;
        //    cout << "Filname for db to open = "<< filename.str() << endl;
        Matrix<int>* M
            = Matrix<int>::fromFile(filename.str().c_str(), _vectorLen,
                                    _windowLenBehind + _windowLenAhead);
        if(M == NULL)
            return new Matrix<int>(_vectorLen, windowLenBehind+windowLenAhead);
        
        if(windowLenBehind == _windowLenBehind &&
           windowLenAhead  == _windowLenAhead)
            return M;
        
        Matrix<int> *M2 = matrixCopy<int, int>(M, 
                                               0, _vectorLen - 1,
                                               _windowLenBehind - windowLenBehind,
                                               _windowLenBehind + windowLenAhead - 1);
        delete M;
        return M2;
    }
}


int SDDBAccessor::myWindowLenBehind() { return _windowLenBehind; }
int SDDBAccessor::myWindowLenAhead()  { return _windowLenAhead;  }
int SDDBAccessor::myVectorLen()       { return _vectorLen;       }
int SDDBAccessor::myNumVectors()      { return _numVectors;      }


void SDDBAccessor::flush() {
  //cerr << "SDDBAccessor::flush()" << endl;
    size_t flushCount = 0;
#pragma omp parallel for    
    for(int i = 0; i < _numVectors; i++) {
        if(_vectors[i] != NULL) {
            //            cerr << "Flushing Vector " << i << endl;
            flush(i);
#pragma omp critical
            flushCount++;
        }
    }
    //    cerr << "Flushed: " << flushCount << " Matricies" << endl;
}


void SDDBAccessor::flush(const int vectorNum) {
    //  cerr << "SDDBAccessor::writeMatrix()" << endl;
  if(_vectors[vectorNum] == NULL)
    return;

  ostringstream filename;
  filename << _name << vectorNum << RECORD_EXTENTION;

  // load up what we've flushed out
  //  cerr << "loading previously flushed data" <<endl;
  //  Matrix<int> *m = Matrix<int>::fromFile(filename.str().c_str(), _vectorLen, _windowLenBehind + _windowLenAhead);
  // add it to what we've currently got
//  if(m != NULL)
//    add(vectorNum, m);
  // flush everything again

  //  cerr << "flush is writing the last version to disk" << endl;
  _vectors[vectorNum]->toFile(filename.str().c_str());
  //  cerr << "Matrix Object Size : " << _vectors[vectorNum]->size() << " for vector " << vectorNum << "\n"; 

  delete _vectors[vectorNum];
  _vectors[vectorNum] = NULL;

  //  if(m != NULL)
  //    delete m;
}


bool SDDBAccessor::removeFromRAM(const int vectorNum) {
  if(!entryInRAM(vectorNum))
    return false;

  flush(vectorNum);
  return true;
}


void SDDBAccessor::close() {
  //  cerr << "SDDBAccessor::close() " << _name << endl;
  flush();
}

