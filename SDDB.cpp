/* 
Copyright (C) 2004,2005,2006,2007,2008,2009  Cyrus Shaoul and Geoff Hollis 

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

//_LARGE_FILES is required for using large files on AIX.... Ha!
#define _LARGE_FILES

#ifdef _OPENMP
/* using conditional compilation to let sequential compilers ignore the omp.h header*/
#include <omp.h>
#endif
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
#include <stdlib.h>
#include <fenv.h>
#include <time.h>
#include "SDDB.h"
#include "MatrixUtils.h"
#include "utilities.h"
#include "Exception.h"
#include "sys/stat.h"

using namespace std;

// ****************************************************************************

// destructor: Remove accessor
SDDB::~SDDB()
{ 
    if(_accessor != NULL)
    {
        _accessor->close();
        delete _accessor;
    }
}

SDDB::SDDB(const string dbname, const string dbpath) {
    _dbname = dbname;
    _dbpath = dbpath;
    _accessor = NULL;
}

// constructor: open a database to use for update or printsds
void SDDB::load(const string eod, const size_t maxMemory) { 

    _eod = eod;
    
    string dbBase = _dbpath + _dbname + DBINFO_TAG;

    int numVectors;
    int vectorLen;
    int windowLenBehind;
    int windowLenAhead;
    long corpussize;

    readMetaData(dbBase, numVectors, vectorLen, windowLenBehind, windowLenAhead, corpussize);

    cerr << "Current Corpus Size = " << corpussize << " words."<< endl;
    _corpussize = corpussize;

    cerr << "Checking for DB......"<< endl;
    if (!dbReady(_dbname)) {
        ostringstream buffer;
        buffer << "An intact database was not found for " << _dbname << ". Cannot continue without database. Exiting." ;
        throw Exception(buffer.str());
    }

    string accname;

    //  sprintf(accname, "%s/%s%s", accname, _dbname, ACCESSOR_TAG);
    accname = _dbpath + _dbname + DBDIR_TAG + "/" + _dbname + ACCESSOR_TAG;

    _accessor = new SDDBAccessor(accname.c_str(), numVectors, vectorLen, 
                                 windowLenBehind, windowLenAhead, maxMemory);

    _realAhead = windowLenAhead;
    _realBehind = windowLenBehind; 

    ostringstream dictname;
    dictname << _dbpath << _dbname << DICT_TAG;
    cerr << "Lexicon to be used = " << dictname.str() << endl;

    build_dict_and_freqs(_dict, dictname.str(), _frequency, eod);
    
    build_idMap(_dict, _idMap);

    _numwords = _dict.size();
    cerr << "Number of words in Dictionary = " << _numwords << endl;
    _wordNum = 0;
    _currStep = 0;
    _stepsize = numVectors;
}

// called when creating a brand new database
void SDDB::initialize (const string& dictfile, const int windowLenBehind, const int windowLenAhead, const string& eod) { 
    ofstream out;
    _eod = eod;
    _corpussize = 0;

    string dbBase = _dbpath + _dbname + DBINFO_TAG;

    string dictname = _dbpath + dictfile;
    
    if (dbExists(dbBase)) {
        ostringstream buffer;
        buffer << "A database called " << _dbname << " cannot be initializeed. A file called : "<< dbBase <<" alread exists. Please remove the old data before creating a new database. Exiting." ;
        throw Exception(buffer.str());
    }
    
    build_starting_dict(_dict, dictname, _frequency);
    
    size_t entries = _dict.size(); 
    cerr << "Built a lexicon with " << entries << " words\n";
    
    ostringstream accname;
    accname << _dbpath << _dbname << DBDIR_TAG;
    makeDir(accname.str());
    
    ostringstream LexiconFileName;
    
    LexiconFileName << _dbpath << _dbname << DICT_TAG;
    
    write_dict_and_freqs(_dict, LexiconFileName.str() , _frequency, _eod);
    
    _wordNum = 0;
    
    int numVectors = entries;
    int vectorLen  = entries;
    long corpussize = _corpussize;
    
    writeMetaData(dbBase, numVectors, vectorLen, windowLenBehind, windowLenAhead, corpussize);
    
    _currStep = 0;
    _stepsize = numVectors;
    
    _accessor = 0;
}

void SDDB::GetDocuments(istream &in, const int number, vector<string> &documents) {
    string word;
    ostringstream buffer;
    documents.reserve(number);
    for (int i = 0; i < number; i++) {
        in >> word;
        while ((word != _eod) && (in.good())) {
            buffer << word << " ";
            in >> word;
        }
        documents.push_back(buffer.str());
        //        cout << buffer.str() << endl << ":::::::::"<< endl;
        buffer.str("");
    }
} 

// Reads in a new document from the corpus file. Returns false if document is too small.
bool SDDB::ConvertADocument(istream& in, vector<int>& wordsInDocument, const size_t behind,
                            const  size_t ahead, const int testmode) 
{
    vector<int> processedwords;
    vector<string> compoundword;
    string word;
    Numpair ids;
    // clear last document
    wordsInDocument.clear();
    // add empty window to begining of words in document
    for (size_t i = 0; i < (behind+1); i++ ) {
            wordsInDocument.push_back(NO_WORD);
    }
    // take a word out of the document
    in >> word;
    while ((word != _eod) && (in.good())) {    
	// Clean up the word, including splitting possessives.
	ids = CleanWord(word,_dict);
	if (ids.first) {
	    wordsInDocument.push_back(ids.first);
	    if (ids.second) {
		wordsInDocument.push_back(ids.second);
	    }
        } else {
            //            cerr << word << ";" << _dict[word]<< " ";
	    wordsInDocument.push_back(UNRECOGNIZED_WORD);
	}
        _wordNum++;
        in >> word;
    }
    // add empty window to end of words in document
    for (size_t i = 0; i < (ahead+2); i++ ) {
            wordsInDocument.push_back(END_OF_DOCUMENT);
    }
    //    cerr << " reached the end of document...It had " << wordsInDocument.size()<< " words.\n";

    // for debugging
    //    for (vector<int>::iterator i = wordsInDocument.begin(); i != wordsInDocument.end(); ++i ) {
    //        cerr << *(i) << " ";
    //    }

    size_t minwords;
    if (testmode)
        minwords = 4;
    else
        minwords = MIN_WORDS_PER_DOC;

    if (wordsInDocument.size() > minwords)
    {
        //        cerr << "Doc of " << wordsInDocument.size() << " words." << endl;
        return true;
    }
    else
        return false;
}




// tests a document to see if it is empty
bool SDDB::documentIsNotEmpty(vector<int>& window, size_t behind) {
    if (window.at(behind) == END_OF_DOCUMENT) 
        return false; 
    else
        return true;
}

//  make the window when there is no window
void SDDB::makeWindow(vector<int>& window, vector<int>& wordsInDocument) {
  //    size_t i;
  //  cerr << "Words in Window = " << wordsInDocument.size() << endl;
  assert(window.size() < wordsInDocument.size());
    for (vector<int>::iterator i = window.begin(); i != window.end(); ++i) {
        *(i)=wordsInDocument[0];
        //pop element off front of vector
        wordsInDocument.erase(wordsInDocument.begin(), wordsInDocument.begin()+1);
    }
    //for debugging
    //    for (vector<int>::iterator i = window.begin(); i != window.end(); ++i) {
    //            cerr << *(i) << " ";
    //        }
    //    cerr << endl;
}

// slide the window over one word.
void SDDB::slideWindow(vector<int>& window, vector<int>& wordsInDocument) {
    //remove first element from beginning of window
    window.erase(window.begin(), window.begin()+1);
    //add new element to end of window.
    window.push_back(wordsInDocument[0]);
    //remove first element from Words in Document
    wordsInDocument.erase(wordsInDocument.begin(), wordsInDocument.begin()+1);
    //for debugging
    //    for (vector<int>::iterator i = window.begin(); i != window.end(); ++i) {
    //      cerr << *(i) << " ";
    //    }
    //    cerr << endl;
}

// Adds Cooccurences to raw cooccrence database
void SDDB::addCooccurrences(vector<int>& window, size_t target) {
    
    if(window[target] != UNRECOGNIZED_WORD) {
        // exclude words that are not in current step.
        if ((window[target] >= _currStep) && (window[target] < (_currStep + _stepsize - 1))) {
            //load matrix for target word into RAM if it is not already in RAM
            //            cerr << window[target]  << " must be between" <<  _currStep<< " and " << _currStep + _stepsize - 1 << endl;
            //_accessor->loadToRAM(window[target]);
            
            for(size_t i = 0; (i < window.size()) && (window[i] != END_OF_DOCUMENT); i++){
                
                // we still haven't moved far enough into the 
                // file to have a behind window
                if(window[i] == NO_WORD)
                    continue;
                // if we're at a word we don't recognize, skip it
                else {
                    if(window[i] == UNRECOGNIZED_WORD)
                        continue;
                // we're at the current word ... up its frequency and skip it
                    else {
                        if (i == target){
                            //                            cerr << "Freq for " << _idMap[window[i]] << " increased by 1" << endl;
			  _frequency[window[i]]++;
                        }
                // otherwise let's update the co-occurance database
                        else {
                            // cerr << "adding data for word " << window[target] << " with " << window[i] <<  endl; 
                            _accessor->add(window[target], window[i], i-target, 1);
                        }
                    }
                }
            }
        }
    }
}

void SDDB::update(istream& in, const int testmode) {

    size_t ahead = _accessor->myWindowLenAhead();
    size_t behind = _accessor->myWindowLenBehind(); 
    size_t coOccurSpan = (ahead + behind + 1); // +1 for current word
    vector<int> window(coOccurSpan,-6);
    vector<int> wordsInDocument;
    
    size_t documentCount = 0;
    size_t badDocumentCount = 0;

    bool validDocument;
    vector<string> documents;
    int x = 0;
    _possessive = "'S";
    size_t loopCount = 0;

    while(in.good()) {
	loopCount++;
        documents.clear();
	// Please ignore this stuff.... It didn't work so well......
	//
        // get a large number of documents.
        //        GetDocuments(in,DOC_BATCH_SIZE,documents);
        //        loopCount = loopCount + DOC_BATCH_SIZE;
        // Progress Meter
        //        size_t docBatchSize = documents.size();
        //        cerr << "Read " << docBatchSize << " documents from disk " << timestamp() << endl;
        //#pragma omp parallel for
        //        while (documents[x] != "") {
        //        for (int x=0; x < static_cast<int>(docBatchSize); x++) {
        //            string temp = documents[x];
        //            if (temp.size() < 1) {
        //                cerr << "found empty doc " << endl;
        //                continue;
        //            }

        // Progress Meter
        if(((x+1) % 1000) == 0) {
	  cerr << "Reading Doc number: " << loopCount << " " << timestamp() << endl;
	  //	  cerr << ".";
	  cerr.flush();
        }
        validDocument = ConvertADocument(in, wordsInDocument, behind, ahead, testmode);
        if (validDocument) {
            documentCount++;
            makeWindow(window, wordsInDocument);
            while (documentIsNotEmpty(window, behind)) {
                slideWindow(window, wordsInDocument);
                //#pragma omp critical
                addCooccurrences(window, behind);
            }
        } else {
            badDocumentCount++;
        }
	// increment progress meter
        x++;
    }
    //    cerr << endl << " Finished Processing Doc number: " << x << endl;
    
    cerr << "Finished this processing step.\n " ;
    _corpussize = _wordNum;
    
    //for debugging.. show Matrix
    //for (size_t i = 0; i < _dict.size(); i++) {
    //        cerr << "Word number " << i << ":" << _idMap[i] << " \n";
    //        _accessor->show(i);
    //    }
}


Matrix<int> *SDDB::getMatrix(const char *word,
                             const int windowLenBehind,
                             const int windowLenAhead)
{
  int num;
  if( _dict.find(word) != _dict.end()){
    num = _dict[word];
    return _accessor->getMatrix(num, windowLenBehind, windowLenAhead);
  }
  else
    return NULL;
}


int SDDB::rows() { return _accessor->myNumVectors(); }
int SDDB::columns() { return _accessor->myVectorLen(); }
int SDDB::windowLenBehind() { return _accessor->myWindowLenBehind(); }
int SDDB::windowLenAhead() { return _accessor->myWindowLenAhead(); }


void SDDB::flush() {
    assert(_accessor);
    _accessor->flush();
    _wordNum = 0;
}


//called as the close the db and save state to disk
void SDDB::close() {

  ostringstream dictname;
  dictname << _dbpath << _dbname << DICT_TAG;

  cerr << "Writing out current lexicon and frequencies." << endl;

  write_dict_and_freqs(_dict, dictname.str(), _frequency, _eod);
  

  int numVectors = _accessor->myNumVectors();
  int vectorLen = _accessor->myVectorLen();
  int windowLenBehind = _accessor->myWindowLenBehind();
  int windowLenAhead = _accessor->myWindowLenAhead();
  int corpussize = _corpussize;

  string dbBase = _dbpath + _dbname + DBINFO_TAG;;

  cerr << "Writing out meta-data." << endl;
  writeMetaData(dbBase, numVectors, vectorLen, windowLenBehind, windowLenAhead, corpussize);
 
  _accessor->close();
  delete _accessor;

  _accessor = NULL;
}

bool SDDB::stepUp() {
  _currStep += _stepsize;
  if(_currStep >= rows())
    return false;
  return true;
}


void SDDB::setCurrentStep(int step) {
    _currStep = step;
}


void SDDB::setStepSize(int stepsize) {
    _stepsize = stepsize;
}

pair<string,string> SDDB::createDirectories (const string outputpath, const bool wordsout) {
    // creates Directories needed for output of results.

    string currenttime;
    currenttime = timestamp();
    
    string outputdir;
    outputdir = outputpath + "output." + currenttime;

    int resultcode = 0;
    resultcode = makeDir(outputdir);

    // Code added to prevent time-based file name collision problems
    // This can happen when multiple HiDEx processes begin in the same working 
    // directory at the same time!
    if (resultcode < 0) {
        wait((rand()%10)+20);
        currenttime = timestamp();
        outputdir = outputpath + "output." + currenttime;
        resultcode = makeDir(outputdir);
        if (resultcode < 0) {
            ostringstream buffer;
            buffer << "Could not create directory: " << outputdir << " ...Exiting.";
            throw Exception(buffer.str());
        }
    }
    
    string subdir = outputdir + "/wordneighborhoods";
    if (wordsout) {
      // Subdirectory for word neighborhoods
      resultcode = makeDir(subdir);
      if (resultcode < 0) {
        ostringstream buffer;
        buffer << "Could not create directory: " << subdir << " ...Exiting.";
        throw Exception(buffer.str());
      }
    }


    pair<string,string> outputloc;
    outputloc.first = outputdir;
    outputloc.second = subdir;
    return outputloc;
}


void SDDB::printPairs(istream &in,
		      const int context_size, int weightingScheme, 
		      const int windowLenBehind, const int windowLenAhead,
		      const int separate, 
		      const string outputpath,
		      const string metric,
		      const string normalization,
		      const int saveGCM)
{

    size_t pairs_to_print = 0;
    //    size_t num_dimensions = 0;
    int behind = min(windowLenBehind, _realBehind);
    int ahead  = min(windowLenAhead, _realAhead);
    int windowLen = behind + ahead;

    // get all words from lexicon, put them in array of strings.
    idMap words;
    words = _idMap;

    // create the directories
    pair<string,string> outputloc = createDirectories(outputpath,0);

    // collect all of the word pairs we'll be printing distances for
    cerr << "Time is: " <<  timestamp() << endl;
    cerr << "Reading in word pairs of interest" << endl;

    // Load word pairs
    vector<pairdata> results;
    LoadPairs(in, results);
    cerr << "Loaded " << results.size() << " word pairs." << endl;


    vector<pairdata> finalresults;

    // Find words that are not in our lexicon from the word pairs list.
    vector<pairdata>::iterator results_iter;
    for (results_iter=results.begin(); results_iter != results.end(); results_iter++) {
        string key1 = (*results_iter).word1;
        string key2 = (*results_iter).word2;
        if ( _dict.find(key1) == _dict.end() ) {
            cerr << "Skipping Words " << key1 << " & " << key2 << ". Word " << key1 << " does not exist in our dictionary." << endl;
        } else {
            if ( _dict.find(key2) == _dict.end()) {
                cerr << "Skipping Words " << key1 << " & " << key2 << ". Word " << key2 << " does not exist in our dictionary." << endl;
            } else { 
                if ((_frequency[_dict[key1]] > 0) && (_frequency[_dict[key2]] > 0)) {
                    pairs_to_print++;
                    finalresults.push_back(*results_iter);
                }
                else {
                    cerr << "Skipping Words " << key1 << " & " << key2 <<". A Word did not appear in the corpus." << endl;
                }
            }
        }
    }
    results = finalresults;
    cerr << pairs_to_print << " target word pairs found" << endl;

    vector<Float*> vectors(_numwords);
    if (GCMexists(_dbname)) {
      LoadMatrix(vectors);
    } else {
      //create weighting scheme
      vector<int> weightScheme = createWeightScheme(windowLen, behind, weightingScheme);
      //create context vector
      vector<int> context = GenerateContext(context_size, separate);
      //Aggregate Vectors
      AggregateVectors(vectors, separate, context, behind, ahead, weightScheme, normalization);
      if (saveGCM) {
	SaveMatrix(vectors);
      }
    }

    //start calculating distances
    for (results_iter=results.begin(); results_iter != results.end(); results_iter++) {
        string word1 = (*results_iter).word1;
        string word2 = (*results_iter).word2;
        
        int id1 = _dict[word1]; 
        int id2 = _dict[word2]; 
        
        if ((vectors[id1] == NULL) || (vectors[id2] == NULL)) {
            cerr << "Tried to print SDs for, " << word1 << " and " << word2 << ", but we did not have a vector!" << endl;
            continue;
        }
	// get inter-word distance 
        (*results_iter).distance  = CalcSimilarity(vectors, id1, id2, metric);
    }

    string currenttime;
    wait((rand()%10)+5);
    currenttime = timestamp();
    string filename = outputloc.first + "/" + "pair.dists." + currenttime + ".txt";
    ofstream dists;
    dists.open(filename.c_str());
    if(!dists.good()) {
        ostringstream buffer;
        buffer << "Could not create output file: " << filename << " ...Exiting.";
        throw Exception(buffer.str());
    }
    dists << "WORD1\tWORD2\tDIST" << endl;

    for (results_iter=results.begin(); results_iter != results.end(); results_iter++) {
        string word1 = (*results_iter).word1;
        string word2 = (*results_iter).word2;
        Float distance = (*results_iter).distance;
        dists << word1 << "\t" << word2 << "\t" << distance << endl; 
    }
    dists.close();

    // Open global output file.
    filename = outputloc.first + "/global.txt";
    ofstream global;
    global.open(filename.c_str());
    if(!global.good()) {
        ostringstream buffer;
        buffer << "Could not create output file " << filename <<" Exiting.";
        throw Exception(buffer.str());
    }

    global << outputloc.first << "\t" << context_size << "\t" << weightingScheme << "\t" << windowLenBehind << "\t" <<  windowLenAhead << "\t";
    global << separate << "\t" << endl;
    global.close();

    cerr << "Done writing results.\nBeginning to delete matrix from memory." << endl;
    for(size_t i = 0; i < _numwords; i++)
        if(vectors[i] != NULL)
            delete [] vectors[i];
    cerr << "Success. Finished processing pairs at " << timestamp() << endl;
    return;
}


int SDDB::printSDs(istream &in,
		   const int context_size, 
		   int weightingScheme, 
 		   const string metric, const string normalization,
		   const int windowLenBehind, const int windowLenAhead,
		   const size_t neighbourhood_size, const int usezscore, 
		   //		   const int useldrt, 
		   const int separate, 
		   const double percenttosample, const int wordlistsize,
		   const string outputpath, 	       
		   const int saveGCM, 
		   const string configdata
		   )
{
    int words_to_print = 0;
    int behind = min(windowLenBehind, _realBehind);
    int ahead  = min(windowLenAhead, _realAhead);
    int windowLen = behind+ahead;
    Float threshold = 0.0;
    Float average = 0.0;
    Float stddev = 0.0;
    // corpus size per million
    Float PerMillionDivisor = static_cast<Float>(_corpussize)/1000000.0;

    // get all words from lexicon, put them in array of strings.
    idMap words;
    words = _idMap;

    // create the directories
    pair<string,string> outputloc = createDirectories(outputpath,1);

    cerr << "Time is: " <<  timestamp() << endl;
    cerr << "Reading in words of interest" << endl;

    // collect all of the words we'll be printing distances for
    vector<resultdata> results;
    //load words
    LoadWords(in, wordlistsize, results);
    cerr << "Loaded " << results.size() << " words." << endl;

    vector<resultdata> finalresults;

    //Dictionary *words_of_interest = dictNew(50000);
    Dictionary words_of_interest; 

    // Find words that are not in our lexicon from word list.
    vector<resultdata>::iterator results_iter;
    for (results_iter=results.begin(); results_iter != results.end(); results_iter++) {
        string key = (*results_iter).word;
        if ( _dict.find(key) != _dict.end()) {
            if (_frequency[_dict[key]] > 0) {
                words_to_print++;
                //      cerr << key << " " << id << endl;
                words_of_interest[key] = _dict[key];
                finalresults.push_back(*results_iter);
            }
            else {
                cerr << "Skipping Word, " << key << ". Word did not appear in corpus." << endl;
            }
        }
        else {
            cerr << "Skipping Word, " << key << ". Word does not exist in our dictionary." << endl;
        }
    }
    results = finalresults;

    cerr << words_to_print << " target words found" << endl;

    // Storage structure for GCM.
    vector<Float*> vectors(_numwords);

    // Look for pre-calculated GCM, and if it exists, load it.
    if (GCMexists(_dbname)) {
      LoadMatrix(vectors);
    } else {
      // If there is no pre-calculated GCM, begin creating one.
      vector<int> weightScheme = createWeightScheme(windowLen, behind, weightingScheme);
      if (weightScheme[0] == 0) {
        return -1;
      }
      // Create context vector
      vector<int> context = GenerateContext(context_size, separate);
      // Aggregate all the raw vectors in the database, and create a GCM.
      AggregateVectors(vectors, separate, context, behind, ahead, weightScheme, normalization);
      // If the user desires, save the CGM for later re-use.
      if (saveGCM) {
	SaveMatrix(vectors);
      }
    }

    //debugging....
    //    for (int x = 0; x < _numwords; x++) 
    //        for (int y = 0; y < num_dimensions; y++) 
    //            if (vectors[x][y] != 0.0) {
    //                cerr << "At: " << x << " " << y << " value: " << vectors[x][y] << endl; 
    //            }

    if (usezscore) {
        cerr << "Generating random distances. Time is: " << timestamp() << endl;
        threshold = GenerateStandardDev(percenttosample,vectors,average,stddev,metric);
        cerr << " Got Threshold, it is: " << threshold  <<  " . Time is: " << timestamp() << endl;
    }
    else    {
        threshold = 0.0;
    }
  
    // ****************************************
    // START CALCULATING VECTOR SIMILARITIES.
    // ****************************************
    // Calculate the similarity of each vector
    // with all other vectors and print out the results

    // Generate our list of words we want similarities

    string filename;

    // Open parameter output file.
    filename = outputloc.first + "/parameters.txt";
    ofstream parms;
    parms.open(filename.c_str());
    if(!parms.good()) {
        ostringstream buffer;
        buffer << "Could not create output file." << filename << "  Exiting.";
        throw Exception(buffer.str()); 
    }
    parms << configdata << endl;
    parms << "Threshold: " << threshold << endl;
    parms.close();

    // Open arc output file.
    filename = outputloc.first + "/arcs.txt";
    ofstream arcs;
    arcs.open(filename.c_str());
    if(!arcs.good()) {
        ostringstream buffer;
        buffer << "Could not create output file." << filename << "  Exiting.";
        throw Exception(buffer.str()); 
    }
    arcs << "WORD\tOFREQ\tARC\tNCount\tInverseNCount" << endl;
    // set precision
    arcs.precision(15);

    // Open global output file.
    filename = outputloc.first + "/global.txt";
    ofstream global;
    global.open(filename.c_str());
    if(!global.good()) {
        ostringstream buffer;
        buffer << "Could not create output file " << filename <<" Exiting.";
        throw Exception(buffer.str());
    }


    // Loops through the list of words to have SDs calculated.
    //    size_t counter = 1;
    cerr << "Completed Processing these words: ";




    // Output loop
    for (Dictionary::iterator i = words_of_interest.begin(); i != words_of_interest.end(); i++) {
      string word = i->first;
      size_t id   = static_cast<size_t>(i->second);

      NeighborsVector neighbours;

        if(vectors[id] == NULL) {
	  cerr << endl << "Tried to print SDs for, " << word
	       << ", but we did not have a vector for it!" << endl;
	  continue;
        }

        // Main distance calculation loop. 
        int limit = static_cast<int>(_numwords);
#pragma omp parallel for 
        for (int j = 0; j < limit; j++) {
            // don't get the distance from ourself
            if(id == static_cast<size_t>(j))
                continue;
            // don't get the distance from words with no vectors
            if(vectors[j] == NULL)
                continue;
	    // get distance
            Float dist =  CalcSimilarity(vectors, id, j, metric);
	    // atomic operation when writing.
#pragma omp critical
            neighbours.push_back(NeighborhoodEntry(dist,j));
        }
        sort(neighbours.begin(), neighbours.end(), RevSort());
        //    remove all neighbors we don't need (gt than MAXNEIGHBOURS)
        if (neighbours.size() > MAXNEIGHBOURS) {
            size_t extras = neighbours.size() - MAXNEIGHBOURS;
            for (size_t m = 0; m < extras; m++) {
                neighbours.pop_back();
            }
        }
        // print out our neighbours.
        filename = outputloc.second + "/" + word + ".nbr.txt" ;
        ofstream nbrs;
        nbrs.open(filename.c_str());
        nbrs.precision(10);
        Float neighbourhood_distance = 0;
        Float neighbourhood_frequency = 0;
        Float freqpermillion = 0;
	//        Float word_magnitude = 0;
        int clustercount = neighbourhood_size;
        Float ARC = 0;
        if (!usezscore) {
	  nbrs << "WORD\tNEIGHBOR\tOFREQ\tSIMILARITY" << endl;
	  // Get standard neighborhood
	  for(size_t j = 0; j < neighbourhood_size; j++) {
	    if(neighbours[j].second != -1) {
	      neighbourhood_distance += neighbours[j].first;
	      freqpermillion = static_cast<Float>(_frequency[neighbours[j].second]) / PerMillionDivisor;
 	      nbrs << word << "\t" 
 		   << words[neighbours[j].second] << "\t"
 		   << freqpermillion << "\t"
 		   << neighbours[j].first << endl;
	    }
	  }
        }
        else {
            // get z-cluster neighborhood
            nbrs << "WORD\tNEIGHBOR\tOFREQ\tSIMILARITY" << endl;
            clustercount = 0;
            for (size_t j=0; j < MAXNEIGHBOURS; j++) {
                if(neighbours[j].second != -1) {
                    if (neighbours[j].first >= threshold) {
                        clustercount += 1;
                        neighbourhood_distance += neighbours[j].first;
                        freqpermillion = static_cast<Float>(_frequency[neighbours[j].second]) / PerMillionDivisor;
                        neighbourhood_frequency += freqpermillion;
                        nbrs << word << "\t" 
			     << words[neighbours[j].second] << "\t"
                             << freqpermillion << "\t"
                             << neighbours[j].first << endl;
                    } 
                    else   {
                        break;
                    }
                }
            }
        }
        // Calculate ARC
        if (clustercount == 0) {
            ARC = neighbours[0].first;
            // Get ARC Zscore
            //      ARC = (ARC - average)/stddev;
        } 
        else {
            if (!usezscore) {
                ARC = ( neighbourhood_distance / static_cast<Float>(neighbourhood_size) ) ;
            } else {
                ARC = ( neighbourhood_distance / static_cast<Float>(clustercount) ) ;
            }
        }
        nbrs.close();

        // print log entry
        Float inverseclustercount = 1.0 / (1.0 + static_cast<Float>(clustercount));
	//#pragma omp critical
        arcs << word << "\t" << (static_cast<Float>(_frequency[id]) / PerMillionDivisor) << "\t" << ARC << "\t" << clustercount << "\t" << inverseclustercount << endl;
	cerr << word << " , ";
    }
    cerr << endl;

    // Write out parameters.
    global << "OutputDir\tContextSize\tWeightingSheme\tWinBehind\tWinAhead\tSeparate\tPercentToSample\tMetric\tNormalization" << endl;
    global << outputloc.first << "\t" << context_size << "\t" << weightingScheme << "\t" << windowLenBehind << "\t" <<  windowLenAhead << "\t";
    global << separate << "\t" << percenttosample << "\t";
    global << metric << "\t" << normalization << endl;
    //    cerr << "Time is: " << timestamp() << endl;
    arcs.close();
    global.close();

    //    delete [] neighbours, vectors;
    //    cerr << "Trying to delete GCM in memory....." << endl;
    for(size_t i = 0; i < _numwords; i++)
        if(vectors[i] != NULL)
	  delete [] vectors[i];
    cerr << "Success. Finished processing all words' neighborhoods at " << timestamp() << endl;
    return 0;
}

int SDDB::printVects(istream &in,
                     const int context_size, 
                     int weightingScheme, 
                     const int windowLenBehind, const int windowLenAhead,
                     const int wordlistsize, 
                     const int separate, const string outputpath, const string normalization, 	       
		     const int saveGCM)
{
    int words_to_print = 0;
    int behind = min(windowLenBehind, _realBehind);
    int ahead  = min(windowLenAhead, _realAhead);
    int windowLen = behind+ahead;
   
    // corpus size per million
    Float PerMillionDivisor = static_cast<Float>(_corpussize)/1000000.0;
    cerr << "Corpus size / 1000000 = " <<  PerMillionDivisor << endl;
    //size_t num_dimensions = 0;
    // get all words from lexicon, put them in array of strings.
    idMap words;
    words = _idMap;
    // create the directories
    pair<string,string> outputloc = createDirectories(outputpath,0);
    // collect all of the words we'll be printing distances for
    cerr << "Time is: " <<  timestamp() << endl;
    cerr << "Reading in words of interest" << endl;
    vector<resultdata> results;
    //load words
    LoadWords(in, wordlistsize, results);
    cerr << "Read " << results.size() << " lines." << endl;
    vector<resultdata> finalresults;
    //Dictionary *words_of_interest = dictNew(50000);
    Dictionary words_of_interest; 
    // Find words that are not in our lexicon from word list.
    vector<resultdata>::iterator results_iter;
    for (results_iter=results.begin(); results_iter != results.end(); results_iter++) {
        string key = (*results_iter).word;
        if ( _dict.find(key) != _dict.end()) {
            if (_frequency[_dict[key]] > 0) {
                words_to_print++;
                //      cerr << key << " " << id << endl;
                words_of_interest[key] = _dict[key];
                finalresults.push_back(*results_iter);
            }
            else {
                cerr << "Skipping Word, " << key << ". Word did not appear in corpus." << endl;
            }
        }
        else {
	  if (key != "") {
            cerr << "Skipping Word, " << key << ". Word does not exist in our dictionary." << endl;
	  } 
        }
    }
    results = finalresults;

    cerr << words_to_print << " target words found" << endl;

    vector<Float*> vectors(_numwords);

    if (GCMexists(_dbname)) {
      LoadMatrix(vectors);
    } else {
      vector<int> weightScheme
        = createWeightScheme(windowLen, behind, weightingScheme);
      if (weightScheme[0] == 0) {
        return -1;
      }
      //    return -7;
      //create context vector
      vector<int> context = GenerateContext(context_size, separate);
      //    cerr << "Actual Context size = " << context.size() << endl;
      AggregateVectors(vectors, separate, context, behind, ahead, weightScheme, normalization);
      cerr << "Writing vectors to disk." << endl;
      if (saveGCM) {
	SaveMatrix(vectors);
      }
    }

    string filename;
    // Open Word output file.
    filename = outputloc.first + "/words.txt";
    ofstream wordsOut;
    wordsOut.open(filename.c_str());
    if(!wordsOut.good()) {
      ostringstream buffer;
      buffer << "Could not create output file." << filename << "  Exiting.";
      throw Exception(buffer.str()); 
    }

    // Open vector output file.
    filename = outputloc.first + "/vectors.mat";
    ofstream vectOut;
    vectOut.open(filename.c_str());
    if(!vectOut.good()) {
        ostringstream buffer;
        buffer << "Could not create output file " << filename <<" Exiting.";
        throw Exception(buffer.str());
    }

    // write header line with size of matrix

    vectOut << words_of_interest.size() << " " << _numdimensions  << endl;
    vectOut.precision(8);
    cerr << "Starting Output. Time is: " << timestamp() << endl;

    // Loops through the list of words to have vectors output

    for(Dictionary::iterator i = words_of_interest.begin(); i != words_of_interest.end(); i++) {

        string word  = i->first;
        size_t  id   = static_cast<size_t>(i->second);

        // print word to word list file.
        wordsOut << word << endl;

        if (vectors[id] == NULL) {
            cerr << "Tried to print vector for, " << word
                 << ", but we did not have a vector for it!" << endl;
            continue;
        }

        vectOut << fixed;
       //        vectOut << word << " ";
        // print values to vector file
        for(size_t k = 0; k < _numdimensions; k++) {
            vectOut << vectors[id][k] << " ";
        }
        vectOut << endl;
    }
    wordsOut.close();
    vectOut.close();
    for(size_t i = 0; i < _numwords; i++)
        if(vectors[i] != NULL)
            delete [] vectors[i];
    //  delete[] vectors;
    // no need to delete words_of_interest;
    cerr << "Success. Finished processing at " << timestamp() << endl;
    return 0;
}


Float SDDB::GenerateStandardDev(const Float percenttosample, const vector<Float*> &vectors, Float &average, Float &stddev, const string metric)
{

  Float total = 0.0;
  average = 0.0;

  numberpairvector listpairs;
  // find all pairs that co-occur
  cerr << "Looking for words that actually co-occurred ..."<< endl;
  listpairs.reserve(100000);

  int numwords = static_cast<int>(_numwords);

  #pragma omp for 
  for (int x = 0; x < numwords; x = x + 3) {
    // Skipping every two to reduce memory footprint.
      Numpair newpair;
      if ((x % 5000) == 0) {
          cerr << "Finished looking for cooccuring words in column: " << x << "." << endl; 
      }
      if (vectors[x] == NULL)
          continue;
      //      if (listpairs.size() > 10000000)
      //          break;
      for (size_t y = 0; y < _numdimensions; y++){
          //for (int y = 0; y < static_cast<int>(num_dimensions); y++){
          // get all pairs that co-occure, and are not the same word!
          if ((vectors[x][y] != 0.0) && (x != static_cast<int>(y))) {
              newpair.first = x;
              newpair.second = y;
              //              cerr << "Found Pair:" << newpair.first << " " << newpair.second << " ";
              #pragma omp critical
              {
                  listpairs.push_back(newpair);
              }
          }
      }
  }

  size_t NumPairs = listpairs.size();
  cerr << "Number of Pairs found: "<< NumPairs <<endl;
  cerr << "Time is: " << timestamp() << endl;
  
  Float numpossible;
  numpossible = static_cast<Float>(_numwords) * static_cast<Float>(_numdimensions);
  cerr << "Number of possible : "<< numpossible << ", percent of space with data: " << static_cast<Float>(NumPairs)/numpossible <<endl;
  
  
  // Calculate the amount to calculate based on the percent to sample.
  size_t NumDesired = static_cast<size_t>(NumPairs * percenttosample);
  cerr << "Percent to Sample: "<< percenttosample << endl; 
  cerr << "Number Desired to find similarity for: "<< NumDesired << endl; 
  
  // make random subsampling of set of possible pairs.
  //seed random number generator
  int seed = time(0);
  srand(seed);
  
  unsigned int randompairnumber = 0;
  //  numberpairset pairset;
  //  numberpairmap pairmap;
  numberpairvector pairvector;

  //  cerr << "Maximum size of set = " << pairset.max_size() << endl;
  // Insert unique sample of word pairs
  size_t i = 0;
  size_t setsize = 0;
  //  numberpairset::iterator listpair_iter;  
  while ((setsize < NumDesired) && (i < NumPairs))  {
      i++;
      setsize = pairvector.size();
      if(!i || (i+1) % 10000000 == 0) {
          cerr << "Choosing Random Pair number " << (i+1) << ", Setsize = " << setsize << endl;
      }
      randompairnumber = static_cast<size_t>(rand() % NumPairs);
      Numpair newpair;
      newpair.first = listpairs[randompairnumber].first;
      newpair.second = listpairs[randompairnumber].second;
      pairvector.push_back(newpair);
  }

  cerr << "Found " << setsize << " pairs" << endl;
  cerr << "Time is: " << timestamp() << endl;

  //  delete listpairs;
  listpairs.clear();
  // allocate memory for scores
  vector<Float> scores;
  scores.reserve(10000);
  

  //  ofstream out;
  //  out.open("distdata.txt");
  
  //  int loopcount = 0;
  int pairvectorsize = static_cast<int>(pairvector.size());

#pragma omp parallel for reduction(+:total)
  for (int pindex = 0; pindex < pairvectorsize; ++pindex) {

      //#pragma omp atomic
      //      loopcount++;
      if (!pindex || ((pindex+1) % 100000 == 0)) {
          cerr << "Calculated cooccurring pair similarity " << pindex + 1 << endl; 
      }
      int num1 = pairvector[pindex].first;
      int num2 = pairvector[pindex].second;

      // Get similarity.
      Float similarity =  CalcSimilarity(vectors, num1, num2, metric); 
#pragma omp critical
      scores.push_back(similarity);
      total += similarity;
  }

  //  out.close();
  //  delete pairset;

  cerr << "Time is: " << timestamp() << endl;

  size_t numSimilarities = static_cast<size_t>(scores.size());
  cerr << "Number of randomly selected word pairs = " << numSimilarities << " which is ";
  cerr.precision(6);
  cerr << ((static_cast<Float>(numSimilarities) / static_cast<Float>(NumPairs))*100.0) << "% of those possible" << endl; 
  
  average = total/(static_cast<Float>(numSimilarities));
  cerr << "Average similarity between randomly co-occuring words is : " << average << endl;  

  //calculating the StdDev:  std dev = sqrt (ss/N)
  Float SumSquares = 0.0;
  stddev = 0.0;
  
  for(size_t i = 0; i < numSimilarities; i++) {
    SumSquares += ((scores[i] - average) * (scores[i] - average));
  }
  
  stddev = sqrt(SumSquares / static_cast<Float>(numSimilarities));
      
  cerr << "STDEV is : " << stddev << endl;
  Float threshold = 0.0;
  Float multiplier = 2.0;

  //similarity should be one half stddev below the average
  threshold = (average + (stddev * multiplier));
  while (threshold >= 1.0) {
     multiplier = multiplier - 0.25;  
      cerr << "WARNING: Threshold was too large. Using smaller threshold,  " << multiplier << " STD above the mean." << endl;
      threshold = (average + (stddev * multiplier));
    }
  

  //
  //find threshold in the top 1%
  //  Float rank = 0.000001;
  //  unsigned int indx = static_cast<unsigned int>(numSimilarities * rank );
  //  threshold = scores[indx];

  cerr << "Threshold  was " << threshold << ", which is " << multiplier << " standard deviations above the mean. " << endl;
  //  cerr << "Threshold  is set to " << threshold << ", which is in the top " << rank * 100  << "% of the distances." << endl;
  cerr << "Finished Threshold Calculations"<< endl;

  // Send back the threshold,
  return threshold;
}

void removeAllFiles(const string& dbname, const string& dbpath) {
    removeDBFiles(dbname, dbpath);
}


vector<int> SDDB::GenerateContext(const size_t context_size, const bool separate)
{
  // get a list of all our words
  //  keep only the N most frequent
  //  
  cerr << "Number of words in lexicon = " << _numwords << endl;

  if (_numwords < context_size) {
    ostringstream buffer;
    buffer << "Context Size ("<< context_size << "is too large. Must be less than number of words in lexicon (" << _numwords << ". Please change your settings.";
    throw Exception(buffer.str());
  }

  // Sort our words by frequency
  cerr << "Sorting words in lexicon by frequency" << endl;
  //  cerr << "Time is: " << timestamp() << endl;

  FrequencySorter SortedByFreq;

  for(FrequencyMap::iterator i = _frequency.begin(); i != _frequency.end() ; ++i) {
    if (i->second > 0 ) {
      SortedByFreq.push_back(FrequencyEntry(i->second,i->first));
    }
  }
  cerr << "Size of Non-zero frequency word list = " << SortedByFreq.size() << endl;
  // checking for sufficient data for context size
  
  if (SortedByFreq.size() < context_size) {
    ostringstream buffer;
    buffer << "The size of your context, " << context_size << " is larger than the size of of the number of words used in the corpus, "  <<  SortedByFreq.size() <<"." ;
    throw Exception(buffer.str());
  }

  sort(SortedByFreq.begin(), SortedByFreq.end(), FreqSort2());

  // Double the size of vectors for original HAL (seperate forward and back vectors.)
  cerr << "Building vector context, " ;
  if (((2*context_size) > SortedByFreq.size()) & separate) {
    ostringstream buffer;
    buffer << "You need a smaller context size to use separate forward and backward contexts. Please edit the config file and reduce your context size to something below " <<  SortedByFreq.size()/2 <<" words." ;
    throw Exception(buffer.str());
  }
  if (separate) {
    cerr << "Using separate forward and backward contexts.\n";
    _numdimensions = 2 * context_size;
  }   else {
    _numdimensions = context_size;
  }
    
  cerr << "Context size is " << context_size << " words " << endl;
  cerr << "Size of Global Co-occurrence Matrix is " << _numwords << " by " << _numdimensions << ", containing " << _numwords * _numdimensions << " elements."  << endl;

  // build context array. Put in all id's for more frequent words
  //    vector<int> context;
  vector<int> context(context_size,0);

  size_t z = 0;
  for(FrequencySorter::iterator i = SortedByFreq.begin(); i != SortedByFreq.end(); i++) {
    if (z < context_size) {
      context[z] = (i->second);
    } else {
      break;
    }
    z++;
  }
  //  cerr << "The context vector is size: " << context.size() << endl;
  //  cerr << "The context is : [";
  // FOR DEBUGGING... print context.
  //  for (size_t i=0; i < context.size(); i++) {
  //    cerr << i << " : " << _idMap[context[i]] << "\n";
  //  } 
  //  cerr  << "]" << endl;
  return context;
}


void SDDB::AggregateVectors(vector<Float*> &vectors, const bool separate, vector<int> &context, const int behind, const int ahead, const vector<int> weightScheme, const string normalization)
  // Collect raw counts and create word vectors from them!
{
    // ***STEPS***
    // 1. get a matrix in array form
    // 2. apply weights to raw data in matrix
    // 3. copy data to our collapsed vector
    // 4. normalize vector
    cerr << "Processing word matrices." << endl;
    cerr << "Loading, Aggregating and Normalizing vectors..." << endl;
    cerr << "Time is: " << timestamp() << endl;
    cerr << "Printing one dot per 1000 vectors aggregated: ";

    // Much speeded by parallel execution!
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(_numwords); i++) {
        if(!i || (i+1) % 1000 == 0)
            cerr << "." ;
        Matrix<int> *M = _accessor->getMatrix(i,_realBehind,_realAhead);
        if(!M) 
        {
            vectors[i] = NULL;
            continue;
        }
        if (separate) {
            vectors[i] = computeWeightedCooccurenceVectorSeparate(M, 
                                                                  weightScheme,
                                                                  _realBehind, behind,
                                                                  _realAhead, ahead,
                                                                  context);
        } else 
        {
            vectors[i] = computeWeightedCooccurenceVectorUnified(M, 
                                                                 weightScheme,
                                                                 _realBehind, behind,
                                                                 _realAhead, ahead,
                                                                 context);
        }        
	// Normalize the vector, using the correct normalization algorithm
	normalizeRawVector(i, vectors[i], normalization, context, separate);
        //  Free up some memory and delete the raw matrix.
        delete M;
    }

    cerr << endl;
    //    cerr << "Count of Words with Matricies "<< tempcounter << endl;
    for (size_t j=0; j < _numwords; j++) {
        if (vectors[j] == NULL) {
            cerr << "Found null vector number = " << j << endl;
        }
    }
    cerr << "Finished processing all vectors. The GCM is now ready.\nTime is: " << timestamp() << endl;
}

bool SDDB::acceptable_word_char(char& ch, bool alphaOnly)
//
// tells whether a character is an acceptable
// component of a word or not
//
{
  // spaces separate words
  if(isspace(ch))
    return false;
  else if(!alphaOnly)
    return true;
  else
    return isalpha(ch);
}

string SDDB::strip_non_alpha(string& word, bool sidesOnly)
//
// creates a duplicate of the word, with all or trailing/leading
// non-alphanumeric characters removed.
//
{
    string stripped = "";
    size_t wordlength = word.length();
    //cout << "Stripping word: " << word << endl; 
    if (sidesOnly) {
        size_t i, j;
        // find the starting position...
        for(i = 0; i < wordlength - 1; i++) {
            if (isalpha(word[i])) {
                break; 
            }
        }
        // find the ending position...
        for (j = wordlength - 1; j > i; j--) {
            if (isalpha(word[j])) {
                break; 
            }
        }
        // copy!
        stripped = word.substr(i,(j-i+1));
        //cout << "Stripped word: " << stripped << endl; 
        return stripped;
    }
    else {
        stripped = word;
        for (size_t i = 0; i < stripped.length(); i++) {
            if ( !isalpha(stripped[i])) {     
                stripped.erase(i,1); // remove char
                i = i - 1;
            }
        }
        //cout << "Stripped word: " << word << endl; 
        return stripped;
    }
}

vector<string> SDDB::iscompound(string& localword){
  // Decides of a word is a compound.
  // Currently unused
    vector<string> output;
    string::size_type pos;

    pos = localword.find ("-",1);
    if ((pos != string::npos)) {
        output.push_back(localword.substr(0, pos));
        output.push_back(localword.substr(pos, localword.length()));
        return output;
    } else {
        pos = localword.find ("/",1);
        if ((pos != string::npos)) {
            output.push_back(localword.substr(0, pos));
            output.push_back(localword.substr(pos, localword.length()));
            return output;
        } 
    }
    output.push_back(localword);
    return output;
}



wordpair SDDB::extractWord(string& localword)
// 
// Extracts the next word from the input stream. This function has the same
// functionality as "cin >> word", excepts it additionally makes sure that
// the word buffer is not over-run! And it separates possessive endings.
//
{
    wordpair output;
    output.main = "";
    output.possessive = "";     // this is mainly used for holding "'S" since we can't putback

    size_t localwordlength = localword.length();
    //uppercase
    
    for (size_t j=0; j < localwordlength; ++j) {
      localword[j]=toupper(localword[j]);   
    }
    
    // cout << "localword: " << localword << endl;
    //check if the special posessive words are there.
    if (!((localword == "HE'S") || (localword =="SHE'S") || (localword == "IT'S") || (localword == "HERE'S") || (localword == "THERE'S") || (localword == "WHAT'S"))) {
        //cout << "Looking for possessive edning for: " << localword << endl;
        // we have to check if our word ends in "'S" and it is longer than
        // "'S" alone. If that's the case, we throw back the "'S"
        string::size_type pos;
        pos = localword.find (_possessive,0);
        if ((localwordlength > 2) && (pos != string::npos)) {
            output.possessive = _possessive;
            localword.erase(pos,2);
            //cout << "Found possessive ending for: " << localword << endl;
        }
    }
    //  cout << "adding word" << localword << endl; 
    output.main = localword;
    return output;
}



vector<int> SDDB::preprocess(string& word)
// PREPROCESS:
// Take a string, and search through it for space-separated
// strings of letters. Compare them to entries in our dictionary.
// If the exist, then add them to the output. 
// add a nonword marker to show that we encountered a nonword.
//
{
    wordpair output;
    vector<int> words;

    // cout << "working on word: " << word << endl;
    if ((word.length() < MAX_WORDLEN) || (word.length() > 0)) { 
        output = extractWord(word);
    } else {
        return words;
    }
    
    // If word in Dict
    if (_dict.find(output.main) != _dict.end()) {
        //cout << "Found word in dict: " << output.main << endl;
        words.push_back(_dict[output.main]);
    } else {
        // If word is dirty Try stripping non-letters off the sides of the word (like 500Club -> Club)
        //	cout << word << endl;
        string cleanWord = strip_non_alpha(output.main, true);
        // cout << "Tried cleaning into: " << cleanWord << endl;
        if (cleanWord.empty()) {
            // we encountered no letters ... skip this
            // cout << "word was empty after cleaning" << cleanWord << endl;
            return words; 
        } else {
            // try adding clean word to output
            if (_dict.find(cleanWord) != _dict.end()) {
                words.push_back(_dict[cleanWord]);
                // cout << "Found cleaned word in dict: " << cleanWord << endl;
            } else {
                // Strip all non letters out.
                string newWord = strip_non_alpha(output.main, false);
                // cout << "Tried deep cleaning into: " << newWord << endl;
                if (newWord.empty()) {
                    // cout << "word was empty after more cleaning" << newWord << endl;
                    return words;
                } else {
                    if (_dict.find(newWord) != _dict.end()) {
                        // cout << "Found deepcleaned word in dict: " << newWord << endl;
                        words.push_back(_dict[newWord]);
                    } else {
                        // cout << "Couldn't Find word in dict: " << newWord << endl;
                        words.push_back(NONWORD_MARKER);
                    }
                }
            }
        }
    }
    if (!output.possessive.empty()) {
        words.push_back(_dict[output.possessive]);
    }
    //    for (vector<int>::iterator i = words.begin(); i != words.end(); ++i) {
    //        cout << *(i) << endl;
    //    }
    return words;
}


Float* SDDB::normalizeRawVector(int TargetWord, Float cooccurenceVector[], const string normalization, vector<int> &context, const bool separate )
{
  // Depending on the setting in the configfile, a normalization procedure is chosen. Currently only PPMI and Frequency Ratio are supported.

  // Inverse corpus size for efficient division by N.
  Float InvN = 1.0/(static_cast<Float>(_corpussize));

  // Set size of matrix, depending on separate forw/back vectors or not.
  size_t contextSize = context.size();
  size_t size;

  if (separate) {
    size =  contextSize* 2;
  } else {
    size = contextSize;
  }
  //  cerr << "About to test with " << size << " dimensions" << endl; 
  Float InvTargetWordFrequency = 1.0/(static_cast<Float>(_frequency[TargetWord]) + 1.0);
  Float pContextGivenTarget = 0.0;
  Float pContext = 0.0;
  Float PPMI = 0.0;

  // Do Positive Pointwise Mutual Information 
  if (normalization == "PPMI") {
    for (size_t j=0; j < size; j++) {
      if (cooccurenceVector[j]) {
	pContextGivenTarget = cooccurenceVector[j]*InvTargetWordFrequency;
	//	cerr << "pContext|Target: " << pContextGivenTarget << endl;
	if (j > contextSize) {
	  pContext = (static_cast<Float>(_frequency[context[j-contextSize]]) + 1.0) * InvN;
	} else {
	  pContext = (static_cast<Float>(_frequency[context[j]]) + 1.0) * InvN;
	}
	//	cerr << "pContext: " << pContext << endl;
	PPMI = log(pContextGivenTarget/pContext);
	if (PPMI < 0.0)
	  PPMI = 0.0;
	  // cerr << "Negative PMI :" << PMI << " for item " << TargetWord << endl;
	cooccurenceVector[j] = PPMI;
      }
    }
  }
  else {
    // Default: normalize vectors by dividing by ofreq.
    // First: make reciprocal of 1 over non-zero frequncy (Optimization to avoid many division ops)
    for(size_t i = 0; i < size; i++) {
      if (cooccurenceVector[i]) {
	//      cerr << "Changed unnormalized " << cooccurenceVector[i] << " to ";
	cooccurenceVector[i] *= InvTargetWordFrequency;
	//      cerr << cooccurenceVector[i] << endl;
      }
    }
  } 
  return cooccurenceVector;
}



// Calculate the similarity of two vectors using choice of algorithms.
Float SDDB::CalcSimilarity(const vector<Float*> &vectors, int w1, int w2, const string algorithm) 
{
  Float similarity = 0.0;

  if (algorithm == "Cosine") {
    // Take the dot product of the two vectors, and divide it by the product of
    // the two magnitudes.
    Float n1 = 0.0;
    Float len1 = 0.0;
    Float len2 = 0.0;
    for(size_t k = 0; k < _numdimensions; k++) {
      n1 += (vectors[w1][k] * vectors[w2][k]);
      len1 += (vectors[w1][k] * vectors[w1][k]);
      len2 += (vectors[w2][k] * vectors[w2][k]);
    }
    if ((len1 > 0.0) && (len2 > 0.0)) {
      similarity =  (n1/(sqrt(len1) * sqrt(len2)));
    } else {
      similarity = 0.0;
    }
  } else if (algorithm == "Correlation" ) {
    // Take the correllation between the two vectors.
    Float prod = 0.0;
    Float sum1 = 0.0;
    Float sum2 = 0.0;
    Float sumsq1 = 0.0;
    Float sumsq2 = 0.0;
    Float N = static_cast<Float>(_numdimensions);
    for(size_t k = 0; k < _numdimensions; k++) {
      prod += vectors[w1][k] * vectors[w2][k];
      sum1 += vectors[w1][k];
      sum2 += vectors[w2][k];
      sumsq1 += vectors[w1][k] * vectors[w1][k];
      sumsq2 += vectors[w2][k] * vectors[w2][k];
    }
    if ((sum1 != 0) && (sum2 != 0)) {
      // Final calculation of correlation
      similarity = abs(((N * prod) - (sum1 * sum2))/(sqrt(((N * sumsq1) - (sum1 * sum1)) * ( (N * sumsq2) - (sum2 * sum2))))); 
    } 
  } else if (algorithm == "CityBlock") {
    // Inverse City Block
    for(size_t k = 0; k < _numdimensions; k++) {
      similarity += abs(vectors[w1][k] - vectors[w2][k]);
    }
    similarity = 1/((similarity * similarity) + 1.0);
  } else {
    // Inverse Euclidean Distance: inverse of the sqrt of the sum of the squared differences
    for(size_t k = 0; k < _numdimensions; k++) {
      similarity += (vectors[w1][k] - vectors[w2][k]) * (vectors[w1][k] - vectors[w2][k]);
    }
    similarity = 1/(sqrt(similarity) + 1.0);
  }
  return similarity;
}

// Save a copy of the GCM to disk for later use.
void SDDB::SaveMatrix(const vector<Float*> &vectors) {
  //open file handlers
  string GCMpath = _dbpath + _dbname + GCM_TAG;
  ofstream GCMOut;
  GCMOut.open(GCMpath.c_str());

  if(!GCMOut.good()) {
    ostringstream buffer;
    buffer << "Could not create output file " << GCMpath <<" Exiting.";
    throw Exception(buffer.str());
    }

  cerr << "Size of GCM = " << _numwords << " by " << _numdimensions << endl;
  // Save Matrix size.
  GCMOut << _numwords << " " << _numdimensions << endl;
  GCMOut.precision(10);
  cerr << "Saving GCM requested. Starting GCM output. This may take a while. Time is: " << timestamp() << endl;
  for (size_t i = 0; i < _numwords; i++) {
    if(((i+1) % 1000) == 0) {
      cerr << ".";
      cerr.flush();
    }
    if (vectors[i]) {
      for(size_t k = 0; k < _numdimensions; k++) {
	GCMOut << vectors[i][k] << " ";
      }
    }
    else {
      for(size_t k = 0; k < _numdimensions; k++) {
	GCMOut << "0 ";
      }
    }
    GCMOut << endl;
  }
  GCMOut.close();
  cerr << "GCM output completed. Time is: " << timestamp() << endl;
}

void SDDB::LoadMatrix(vector<Float*> &vectors) {

  size_t GCMwords = 0;
  size_t GCMcontext = 0;
  Float tempValue;
  string GCMpath = _dbpath + _dbname + GCM_TAG;
  ifstream GCMIn;
  

  // Open GCM file and check for failure to open.
  GCMIn.open(GCMpath.c_str());
  if(!GCMIn.good()) {
    ostringstream buffer;
    buffer << "Could not open GCM file " << GCMpath <<" Exiting.";
    throw Exception(buffer.str());
    }
  
  cerr << "Starting to load GCM file. Please be patient as this could take a while. A dot is printed for every 1000 vectors loaded.\n Time is: " << timestamp() << endl;
  // read matrix size information from file header.
  GCMIn >> GCMwords >> GCMcontext;

  if(GCMwords != _numwords) {
    ostringstream buffer;
    buffer << "GCM size mismatch. Size of matrix in file " << GCMpath << " is " << GCMwords << " which does not match lexicon size, " << _numwords << ". Please check GCM file.";
    throw Exception(buffer.str());
    }
  else {
    //    cerr << "NumWords = " << GCMwords << " Context = " << GCMcontext << endl;
    _numdimensions = GCMcontext;
    // Zero out memory.
    // Read through file and add values to matrix in memory.
    for (size_t i = 0; i < _numwords; i++) {
        if(((i+1) % 1000) == 0) {
	  cerr << ".";
	  cerr.flush();
        }
      Float *tempvect = new Float[_numdimensions];
//       for (size_t h = 0; h < GCMcontextsize; h++) {
// 	tempvect[h] = 0.0;
//       }
      vectors[i] = NULL;
      for (size_t j = 0; j < _numdimensions; j++) {
	GCMIn >> tempValue;
	if ((tempValue >= 0.0) && (!GCMIn.eof())) {
	  tempvect[j] = tempValue;
	} else {
	  throw Exception("GCM file read error! Check your GCM file for corruption. Exiting.");
	}
      }
      vectors[i] = tempvect;
      //      delete [] tempvect;
    }
    cerr << "Finished loading GCM file. Time is: " << timestamp() << endl;
  }
  GCMIn.close();
}
