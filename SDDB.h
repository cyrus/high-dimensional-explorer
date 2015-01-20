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
// This class handles the creation and maintenance
// of Semantic Distance databases.
//

#ifndef SDDB_H
#define SDDB_H

#include <algorithm>
#include <map>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <string>
#include <list>
#include <set>
#include <sstream>
#include <vector>
#include <algorithm>
#include <utility>
#include "SDDBAccessor.h"
// #include "dict.h"



// ****************************************************************************
// Defines, typedefs, and macros
// ****************************************************************************
#define DBDIR_TAG                      ".data"
#define DICT_TAG                       ".dict"
#define VAR_TAG                        ".var"
#define ACCESSOR_TAG                   ".db"
#define DBINFO_TAG                     ".dbinfo"
#define GCM_TAG                        ".gcm"
#define CONTEXT_TAG                        ".context"
#define MIN_WORD_VAL                       0
#define UNRECOGNIZED_WORD                 -1
#define END_OF_DOCUMENT                   -2
#define NO_WORD                           -3
#define SCALEFACTOR                       10.0 
#define MAXNEIGHBOURS                     80000
#define BLOCKSIZE                         8
#define DOC_BATCH_SIZE                    100000
#define MIN_WORDS_PER_DOC                 2
//for preprocessor
#define MAX_WORDLEN                        40
#define NONWORD_MARKER                    -1


struct Settings
{
  int contextSize;
  string corpusFilename;
  string dbname;
  string dictFilename;
  string eod;
  int maxWindowAhead;
  int maxWindowBehind;
  int neighbourhoodSize;
  double percenttosample;
  bool separate;
  int stepsize;
  bool usezscore;
  int weightingScheme;
  string metric;
  bool saveGCM;
  string normalization;
  int windowLenAhead;
  int windowLenBehind;
  string wordlistfilename;
  int wordlistsize;
  string outputpath;
  string dbpath;
  string multipleFiles;
  size_t maxMemory;
  bool normCase;
  bool englishContractions;
  bool useVariance;
};


//typedef float Float;
typedef double Float;

//structs
struct resultdata {
  string word;
  Float LDRT;
  Float ANS;
  Float InverseNcount;
};

struct pairdata {
  string word1;
  string word2;
  Float distance;
  Float ANS1;
  Float ANS2;
  Float InverseNcount1;
  Float InverseNcount2;
};

struct wordpair
{
  string main;    // main word
  string possessive;     // posessive
};

/* class Numpair */
/* { */
/*  public: */
/*     int _x;  */
/*     int _y;  */

/*     Numpair(int x, int y) */
/*         : _x(x), */
/*           _y(y) */
/*     {         */
/*     } */

/*     bool operator<(const Numpair& num) const */
/*     { */
/*         return (_x < num._x || (!(num._x < _x) && _y < num._y)); */
/*     } */
/* }; */

typedef pair<unsigned int,unsigned int> Numpair;
typedef set<Numpair> numberpairset;
typedef map<Numpair,bool> numberpairmap;
typedef vector<Numpair> numberpairvector;

// Dictionary has keys of words and values of ID

typedef map<string, int> Dictionary;

// Maps word to their frequencies
typedef map<int, size_t> FrequencyMap;
typedef map<int, Float> VarianceMap;

//Context type for sorting words by frequency or variance
typedef pair<Float, size_t> ContextEntry;
typedef vector<ContextEntry> ContextSorter;

//idMap has keys of IDs and values of words (For reverse lookup)
typedef map<int,string> idMap;

// Vector of distances and words
typedef pair<Float, int> NeighborhoodEntry;
typedef vector<NeighborhoodEntry> NeighborsVector;

// Vector of distances and words
typedef pair<Float, int> NeighborhoodEntry;
typedef vector<NeighborhoodEntry> NeighborsVector;


//typedef vector<const DictPair *> DictVector;

//***************************************************

class RevSort
{
public:
    bool operator() (const NeighborhoodEntry & a, const NeighborhoodEntry & b) const
    {
        return a.first > b.first;
    }
};

class FreqSort
{
public:
    bool operator() (const NeighborhoodEntry & a, const NeighborhoodEntry & b) const
    {
        return a.first < b.first;
    }
};

class FreqSort2
{
public:
    bool operator() (const ContextEntry & a, const ContextEntry & b) const
    {
        return a.first > b.first;
    }
};

class SDDB
{
public:
  static const int FLAT = 0;
  static const int RAMPED_LINEAR = 1;
  static const int RAMPED_QUADRATIC = 2;
  static const int FORWARD_RAMP = 3;
  static const int BACKWARD_RAMP = 4;
  static const int INVERSE_RAMP = 5;
  static const int INVERSE_QUADRATIC = 6;
  static const int SECOND_WORD = 7;
  static const int THIRD_WORD = 8;
  static const int FOURTH_WORD = 9;


  //
  // Instantiate the object
  //
  SDDB(const string dbname, const string dbpath);

  //
  // close everything down and delete everything from RAM
  //
  ~SDDB();

  //
  // Create an SDDB with these settings
  //
  //
  // Load the SDDB defined by this filename into memory
  //
  void load(const string eod, const size_t maxMemory);

  void initialize(const string& dictfile, const int windowLenAhead, const int windowLenBehind, const string& eod);
  //
  // Process the co-occurances in this file.
  // If this is a concatination of more than one file, we will stop processing
  // it as if it were one file as soon as eof is encountered
  //
  void update(istream& in, const int testmode);

  void process2(istream& in, string EOD);
  //
  // Return a Semantic distance matrix for the given word.
  //
  Matrix<int> *getMatrix(const char *word,
                         const int windowLenBehind,
                         const int windowLenAhead);


  //
  // How many entries are in this database?
  //
  int rows();


  //
  // How long is each vector?
  //
  int columns();


  //
  // What is the size of the behind sliding window?
  //
  int windowLenBehind();


  //
  // What is the size of the ahead sliding window?
  // 
  int windowLenAhead();

  //
  // Return a copy of every dictionary entry we have.
  // sets num to the number of entries returned
  //
  //  vector<string> getEntries(int &num);

  //
  // Flush all pending actions
  //
  void flushDB();

  //
  // Close the SDDB
  // 
  void close();

  //
  // Set the minimum vector num we are collecting co-occruance
  // counts for
  // 
  void setCurrentStep(const int step);


  //
  // Set the maximum vector num we are collecting co-occurance
  // counts for, where max = currentStep + stepsize - 1
  //
  void setOptions(const Settings settings);

  // sets the end of document marker
  void setEOD(const string EOD);


  //
  // Increment to the next step of vectors. Returns true if there
  // is another step to work on, and returns false otherwise
  //
  bool stepUp();
  
  //
  // Ahhh, at long last! This is what we've been waiting for!
  // Print all of the SDs for all of our words, plus neighbours,
  // to our output directory
  //
  void printPairs(istream &in,
                  const int context_size, int weightingScheme, 
                  const int windowLenBehind, const int windowLenAhead,
                  const int separate, 
                  const string outputpath,
		  const string metric, 
		  const string normalization,
		  const int saveGCM
		  );

  int printSDs(istream &in,
	       const int context_size, 
	       int weightingScheme, 
	       const string metric, 
	       const string normalization,
	       const int windowLenBehind, 
	       const int windowLenAhead,
	       const size_t neighbourhood_size,
	       const int usezscore,
	       const int separate,
	       const double percenttosample, 
	       const int wordlistsize, 
	       const string outputpath,
	       const int saveGCM,
	       const string configdata
	       ); 

  int printVects(istream &in,
                 const int context_size, 
                 int weightingScheme, 
                 const int windowLenBehind, const int windowLenAhead,
                 const int wordlistsize, 
                 const int separate, const string outputpath,
		 const string normalization,		     
		 const int saveGCM
		 );
  
  Float GenerateStandardDev(const Float percenttosample, 
                            const vector<Float*> &vectors, Float &average, Float &stddev,
			    const string metric);
  
  vector<int> GenerateContext(const size_t context_size, const bool separate);
  
  void AggregateVectors(vector<Float*> &vectors, const bool separate, vector<int>& context, const int behind, const int ahead, const vector<int> weightScheme, const string normalization);
  
  // Reads in a new document from the corpus file.
  bool ConvertADocument(istream& in, vector<int>& wordsInDocument, const size_t behind, const size_t ahead, const int testmode, string lang);
  
  //get a set of N documents from the corpus
  void GetDocuments(istream &in, const int number, vector<string> &documents);

  // tests a document to see if it is empty
  bool documentIsNotEmpty(vector<int>& wordsInDocument, size_t ahead) ;
  
  //  make the window when there is no window
  void makeWindow(vector<int>& window, vector<int>& wordsInDocument) ;
  
  // slide the window over one word.
  void slideWindow(vector<int>& window, vector<int>& wordsInDocument) ;
  
  // Adds Cooccurences to raw cooccrence database
  void addCooccurrences(vector<int>& window, size_t target) ;

  // creates Directories needed for output of results.
  pair<string,string> createDirectories (const string outputpath, const bool wordsout);
  
  // nornamlizes vectors
  Float* normalizeRawVector(int TargetWord, Float cooccurenceVector[], const string normalization, 
			  vector<int> &context, const bool separate);

  // calculates the similarity between two vectors
  Float CalcSimilarity(const vector<Float*> &vectors, int w1, int w2, const string algorithm);

  // Saves a copy of the GCM to disk.
  void SaveMatrix(const vector<Float*> &vectors);

  // Loads a copy of the GCM to disk.
  void LoadMatrix(vector<Float*> &vectors);

 private:
  string _dbname;
  string _dbpath;
  
  /** The accessor for all of the words */
  SDDBAccessor* _accessor;
  
  /** The dictionary we use to hold word - id key-value mappings */
  Dictionary _dict;
  
  /** Frequency and variance counts for all of the words */
  FrequencyMap _frequency;

  // Variables for calculating variance
  VarianceMap _variance;
    
  // map from ID back to word.
  idMap _idMap;
  
  /** How many words have we processed so far? */
  long _wordNum;
  
  /** What step are we working on? */
  int _currStep;

  /** How big are our steps? */
  int _stepsize;
  
  /** How big is our corpus? */
  long _corpussize;
  
  /** How many words are there in our lexicon? */
  size_t _numwords;

  /** How many context dimensions will we use? */
  size_t _numdimensions;

  /** What is the DB Forw window available */

  int _realAhead;

  /** What is the DB Back window available */
  int _realBehind;
  
  // End of Document string
  string _eod;

  //possessive
  string _possessive;

  //norm case?
  bool _normCase;

  //possessive
  bool _englishContractions;

  //possessive
  bool _useVariance;


};

void removeAllFiles(const string& dbname, const string& dbpath);


template <class T>
Float * collapseMatrix(Matrix<T> *M, 
                        const int *weight_scheme,
                        const int realBehind, const int behind,
                        const int realAhead,  const int ahead,
                        const int *context, const int context_size, const int num_dimensions, 
                        const Float wordfrequency, const int separate);


/* 
void addtoresults(vector<resultdata> &results, string word, double ANS, int Ncount);
std::string getcorrelation(vector<resultdata> &results);
*/

#endif // SDDB_H

























