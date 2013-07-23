/* 
Copyright (C) 2004,2005,2006,2007,2008  Cyrus Shaoul and Geoff Hollis 

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

#define _LARGE_FILES 1

// Utility functions for HIDEX
// ****************************************************************************
// 
// ****************************************************************************

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
#include <errno.h>
#include "SDDB.h"
#include "Exception.h"
#include "Filesystem.h"
#include "utilities.h"
#include "string.h"

#include <unistd.h>
#include <unistr.h>
#include <uninorm.h>
#include <unitypes.h>
#include <unictype.h>
#include <unistdio.h>
#include <unicase.h>

bool
// 
// Extracts the next word from the input stream. This function has the same
// functionality as "cin >> word", excepts it additionally makes sure that
// the word buffer is not over-run!
//
getWord(istream& in, string& word)
{
    int  i = 0;
    char ch;
    word = "";

    while(in.get(ch))
    {
        if(!isspace(ch)) {
            if(i > MAX_WORDLEN) {
                in.putback(ch);
                break;
            } else {
                word += ch;
                i++;
            }
        }
        else {
            break;
        }
    }
        
    if (i > 0) 
        return true;
    else
        return false;
}
    
void
//
// Copy the file at source, to dest... Not currently used.
//
copy_file(const char *source, const char *dest) {

  ifstream in(source);
  ofstream out(dest);

  char ch;
  while(in.get(ch))
    out << ch;

  in.close();
  out.close();
}


void 
// Loads pairs of words to find distances for.
LoadPairs(istream& in, vector<pairdata> &results, const bool normCase) {
  pairdata temp;
  cerr << "Reading in pairs from wordlist" << endl;
  
  string lang = getLang();
  while (!in.eof()){
    if (in.fail()) {
      throw Exception("Error Reading Input File.. Make sure it is in the correct format. Exiting\n");
    }
    string tempstring1 = "";
    string tempstring2 = "";
    in >> tempstring1 >> tempstring2;
    if (tempstring1 == "---END.OF.DOCUMENT---")
      break;
    if ((tempstring1.length() > MAX_WORDLEN) || (tempstring2.length() > MAX_WORDLEN))   {
      throw Exception("A Word was too long. Exiting\n");
    }
    if (normCase) {
      temp.word1 = downstring(tempstring1,lang);
      temp.word2 = downstring(tempstring2,lang);
    } else {
      temp.word1 = tempstring1;
      temp.word2 = tempstring2;
    }
    temp.distance = 0.0;
    results.push_back(temp);
  }
}

void
//
// Loadwords: loads a file with target words for similarity measurements.
// Will generate a random subsample if desired.
//
//LoadWords(istream& in, const int useldrt, const int wordlistsize, vector<resultdata> &results) {
LoadWords(istream& in, const int wordlistsize, vector<resultdata> &results, const bool normCase) {
  resultdata temp;
  vector<resultdata> wordlist;
  //  bool done_reading_file = false;
  //    cerr << "Reading in words from wordlist" << endl;
  string lang = getLang();
  while (!in.eof()){
    if (in.fail()) {
      throw Exception("Error Reading Input File.. Make sure it is in the correct format. Exiting\n");
    }
    string tempstring = "";
    Float tempdouble = 0.0;
    in >> tempstring;
    if (tempstring == "---END.OF.DOCUMENT---")
      break;
    if (tempstring.length() > MAX_WORDLEN)  {
      throw Exception("Maximum Word Length was exceeded in your word list input file. Exiting");
    }
    temp.ARC = 10000000;
    temp.InverseNcount = 0.0;
    if (normCase) {
      temp.word = downstring(tempstring,lang);
    } else {
      temp.word = tempstring;
    }
    temp.LDRT = tempdouble;
    wordlist.push_back(temp);
  }
  unsigned int numwords = wordlist.size();
  cerr << "Got " << numwords << " from word list. " << endl ;
  if ( numwords< (unsigned int)wordlistsize) {
    results = wordlist;
  } else {
    cerr << " Only " << wordlistsize << " were desired. " << endl;
    cerr << "Generating random subsample of words." << endl;
    // Create subsample
    int seed = time(0);
    unsigned int randomnumber;
    srand(seed);
    set<int> wordset;
    vector<resultdata> subsamplewordlist;
    int i;
    // Generate random set of unique indicies
    for (i = 0; wordset.size() < (unsigned int)wordlistsize; i++) {
      randomnumber = (unsigned int)(rand() % numwords);
      wordset.insert(randomnumber);
      //      cerr << randomnumber << " " << i <<endl;
    }
      // copy random words into new vector
    set<int>::iterator iter;
    cerr << endl << "Picked: ";
    for (iter = wordset.begin(); iter != wordset.end(); iter++) {
      subsamplewordlist.push_back(wordlist[*iter]);
      //      cerr << *iter << endl;
    }
    results = subsamplewordlist;
  }
}
  
void addtoresults(vector<resultdata> &results, string word, Float ARC, Float InverseNcount)
{
  // Add distance data to results vector for later correlation
  vector<resultdata>::iterator results_iter;
  for (results_iter=results.begin(); results_iter != results.end(); results_iter++) {
    if ((*results_iter).word == word) {
      (*results_iter).ARC = ARC;
      (*results_iter).InverseNcount = InverseNcount;
      break;
    }
  }
}

std::string getcorrelation(vector<resultdata> &results) {
  char buffer[128];
  vector<resultdata>::iterator results_iter;
  Float size = static_cast<Float>(results.size());
  // sum of scores, sum of squared scores sss
  Float sumARC = 0.0;
  Float sssARC = 0.0;
  Float sumNcount = 0.0;
  Float sssNcount = 0.0;
  Float sumLDRT = 0.0;
  Float sssLDRT = 0.0;
  Float sumprodAL = 0.0;
  Float sumprodNL = 0.0;
  Float sign = 1.0;
  // calculate all the sums.
  for (results_iter=results.begin(); results_iter != results.end(); results_iter++) {
    sumARC += (*results_iter).ARC ;
    sssARC += (*results_iter).ARC * (*results_iter).ARC;
    sumNcount += (*results_iter).InverseNcount;
    sssNcount += (*results_iter).InverseNcount * (*results_iter).InverseNcount;
    sumLDRT += (*results_iter).LDRT ;
    sssLDRT += (*results_iter).LDRT * (*results_iter).LDRT ;
    sumprodAL += (*results_iter).ARC  * (*results_iter).LDRT ;
    sumprodNL += (*results_iter).InverseNcount * (*results_iter).LDRT ;
  }

  // calculate corr btween ARC and LDRT
  Float numeratorAL = 0.0;
  Float denominatorAL = 0.0;
  numeratorAL = (size * sumprodAL ) - (sumARC * sumLDRT);
  if (numeratorAL < 0)
    sign = -1.0;
  numeratorAL *= numeratorAL;
  denominatorAL = ((size * sssARC) - (sumARC * sumARC)) * ((size * sssLDRT) - (sumLDRT * sumLDRT));
  Float r2AL = numeratorAL / denominatorAL;
  Float rAL  = sqrt(r2AL) * sign;

  //reset sign
  sign = 1.0;
  
  // calculate corr btween Ncount and LDRT
  Float numeratorNL = 0.0;
  Float denominatorNL = 0.0;
  numeratorNL = (size * sumprodNL ) - (sumNcount * sumLDRT);
  if (numeratorNL < 0)
    sign = -1.0;
  numeratorNL *= numeratorNL;
  denominatorNL = ((size * sssNcount) - (sumNcount * sumNcount)) * ((size * sssLDRT) - (sumLDRT * sumLDRT));
  Float r2NL = numeratorNL / denominatorNL;
  Float rNL  = sqrt(r2NL) * sign;

  string reply = "";
  //  string reply = "The correlation between ARC and LDRT: r2 = ";
  sprintf(buffer,"%.4f\t%.4f\t%.4f\t%.4f", r2AL, rAL, r2NL, rNL);
  reply += buffer;
  
  //  reply += "The correlation between Ncount and LDRT: r2 = ";
  //  sprintf(buffer,"%.8f r = %.8f\n", r2NL, rNL);
  //  reply += buffer;

  return reply;
  
}


//
// Sort the elements of pairs[] by their value (ascending order). 
// If two  words have the same value, resolve ties by which
// words' letters come first in the alphabet
//
// Now uses qsort.

// int compare (const void * a, const void * b)
// {
//   return ( ((const DictPair*)b)->value - ((const DictPair*)a)->value );
// }

// void value_sort2(const DictPair *pairs[], const int size) {
//     qsort (pairs, size, sizeof(DictPair*), compare);
// }

// void value_sort(const DictPair *pairs[], const int size)
// {
//   for(int i = 0; i < size; i++)
//   {
//     int lowestValPos = i;

//     for(int j = size-1; j > i; j--)
//     {
//       if(pairs[lowestValPos]->value > pairs[j]->value)
// 	lowestValPos = j;
//       else if(pairs[lowestValPos]->value == pairs[j]->value && 
// 	      strcmp(pairs[lowestValPos]->key, pairs[j]->key) > 0)
// 	lowestValPos = j;
//     }

//     if(lowestValPos != i)
//     {
//       const DictPair *temp = pairs[i];
//       pairs[i] = pairs[lowestValPos];
//       pairs[lowestValPos] = temp;
//     }
//   }
// }

void
//
// Reads the contents of the file to the specified dictionary,
// and give each a unique numeric value to use for reference in
// the SDDBAccessors
//
build_starting_dict(Dictionary& D, string filename, FrequencyMap& frequencies, const bool normCase)
{
  ifstream inStream(filename.c_str());
  string word = "";
  string lang = getLang();
  int i = MIN_WORD_VAL;
  cerr << "normCase = " << normCase << endl;
  if (!inStream.good()) {
    throw Exception("Dictionary file could not be opened");
  } else {
    do {
      inStream >> word;
      if (normCase) { 
	word = downstring(word,lang);
      } 
      if(D.find(word) != D.end()) {
	cerr << "Word, " << word << ", exists in dictionary already... Skipping it." << endl;
	continue;
      }
      D[word] = i;
      frequencies[i] = 0; 
      i++;
    } while (inStream);
  }
}


void
//
// Builds a dictionary and gives each word an id. Also fills frequencies
// up with the counts of the words
//
build_dict_and_freqs(Dictionary& D, string filename, FrequencyMap& frequencies, const string eod) 
{
    ifstream inStream(filename.c_str());
    string word;
    int id;
    int freq;
    while (inStream.good()){
        id = 0;
        freq = 0;
        inStream >> word;
        if (word == eod) {
            break;
        }
        inStream >> id;
        inStream >> freq;
        D[word] = id;
        frequencies[id] = freq;
        //        cerr << "Got Word:" << word << " id:" << id << " Freq:" << freq << endl;
    }
}

void build_idMap(Dictionary& D, idMap& words) {
  for(Dictionary::iterator i = D.begin(); i != D.end() ; ++i) {
    words[i->second] = i->first;
  }
}

void
//
// Write a dictionary and the word frequencies to file
//
write_dict_and_freqs(Dictionary& D, string filename, FrequencyMap& frequencies, const string eod)
{
    //    cerr << eod << " is the EOD" << endl;
    ofstream out(filename.c_str());
    if (out.good()) 
    {
        for (Dictionary::iterator i = D.begin(); i != D.end(); i++) {
            out << i->first << "\t" << i->second << "\t" << frequencies[i->second] << endl;
        }
        out << eod << endl;
        out.close();
    }
    else {
        throw Exception("Could not write to the .dict file. Exiting");
    }
}



//
// Sorts an array of type T into ascending order
//
template<class T>
void binary_sort(T array[], int len) {

  if(len <= 1)
    return;
  else if(len == 2) {
    if(array[1] > array[0])
      return;
    else {
      T tmp = array[0];
      array[0] = array[1];
      array[1] = tmp;
      return;
    }
  }
  else {

    T arr1[len/2];
    T arr2[len - (len/2)];
    int i;
    for(i = 0; i < len/2; i++)
      arr1[i] = array[i];
    for(; i < len; i++)
      arr2[i - (len/2)] = array[i];

    binary_sort(arr1, len/2);
    binary_sort(arr2, (len - len/2));

    int j;
    for(i = j = 0; i < len/2 && j < (len - len/2);) {
      if(arr1[i] < arr2[j]) {
	array[i+j] = arr1[i];
	i++;
      }
      else {
	array[i+j] = arr2[j];
	j++;
      }
    }
    for(; j < (len - len/2); j++)
      array[i+j] = arr2[j];
    for(; i < len/2; i++)
      array[i+j] = arr1[i];
  }
}


//
// returns the place in the array the element exists
// at. returns len if the element does not exist
//
template<class T>
int array_place(T *array, T elem, int len) {

   int i;
  for(i = 0; i < len; i++)
    if(array[i] == elem)
      return true;

  return i;
}


void wait ( size_t seconds ) {
    clock_t endwait;
    endwait = clock () + seconds * CLOCKS_PER_SEC ;
    while (clock() < endwait) {}
}

string
// TimeStamper
timestamp () {
  time_t rawtime;
  int month;
  int year;
  struct tm * tinfo;
  stringstream tstamp;

  time( &rawtime );
  tinfo = localtime (&rawtime);
  month = tinfo->tm_mon + 1;
  year = tinfo->tm_year + 1900;
  tstamp << setfill('0');
  tstamp << year << "." <<  setw(2) <<  month << "." << setw(2) << tinfo->tm_mday << "." 
	 << setw(2) <<tinfo->tm_hour << "." << setw(2) << tinfo->tm_min << "." <<  setw(2) << tinfo->tm_sec;
  return tstamp.str(); 
}


/**
 * C++ version std::string style "itoa":
*/

std::string itoa(int value, unsigned int base) {
  const char digitMap[] = "0123456789abcdef";
  std::string buf;
  // Guard:
  if (base == 0 || base > 16) {
    // Error: may add more trace/log output here
    return buf;
  }

  // Take care negative int:
  std::string sign;
  int _value = value;
  if (value < 0) {
    _value = -value;
    sign = "-";
  }
  
  // Translating number to string with base:
  for (int i = 30; _value && i ; --i) {
    buf = digitMap[ _value % base ] + buf;
    _value /= base;
  }
  return sign.append(buf);
}





/////


void add_smallest(word_distance_pair *entries, int len, int id, Float distance)
{
  // add the smallest distance values to the neighbor list.
  // does mini-sort
    for(int i = 0; i < len; i++) {
        if(entries[i].id == -1 || distance < entries[i].distance) {
            for(int j = len-1; j > i; j--) {
                entries[j].id = entries[j-1].id;
                entries[j].distance = entries[j-1].distance;
            }
            entries[i].id = id;
            entries[i].distance = distance;
            return;
        }
    }
}

Float
// Calculate the magnitude of a word vector
get_magnitude(Float *vect, int context_size) {
  Float magnitude = 0;
  for (int i = 0; i < context_size; i++) {
    magnitude += (vect[i] * vect[i]);
  }
  magnitude = sqrt(magnitude);
  return magnitude;
}


void removeDBFiles (const string& dbname, const string& dbpath) {
    string dbBase = dbpath + dbname + DBINFO_TAG;
    string dictfilename = dbpath + dbname + DICT_TAG;
    bool errorState = false;
    string datadirname = dbpath + dbname + DBDIR_TAG;
    string gcmfilename = dbpath + dbname + GCM_TAG;
    
    cerr << "Removing all relevant files now....Please be patient."<< endl;

    if (unlink(dbBase) == 0) {
        cerr << "Removed " << dbBase << endl;
    } else {
        cerr << "Could not remove " << dbBase << ". File does not exist" << endl;
        errorState = true;
    }

     if (unlink(dictfilename) == 0) {
         cerr << "Removed " << dictfilename << endl;
     } else {
         cerr << "Could not remove " << dictfilename << ". File does not exist" << endl;
         errorState = true;
     }

    
    if (rmdir(datadirname) == 0) {
        cerr << "Removed " << datadirname << endl;
    } else {
        cerr << "Could not remove " << datadirname << ". Directory does not exist" << endl;
        errorState = true;
    }
    
    if (unlink(gcmfilename) == 0) {
        cerr << "Removed " << gcmfilename << endl;
    } else {
        cerr << "Could not remove " << gcmfilename << ". File does not exist" << endl;
    }
    
    if (errorState) {
        throw Exception("There were errors detected. Database files were not able to be successfully removed.");
    }
    
}

bool dbExists (const string& dbname) {
    string dictfilename = dbname + DICT_TAG;
    string datadirname = dbname + DBDIR_TAG;

    if ( file_exists(dbname) || file_exists(dictfilename) || dir_exists(datadirname)) {
        cerr << "Found a previous database." << endl;
        return true;
    } else {
        return false;
    }
}

bool GCMexists(const string& dbname) {
    string GCMfilename = dbname + GCM_TAG;
    //    string ContextFilename = dbname + CONTEXT_TAG;
    //    if ( file_exists(GCMfilename) || file_exists(ContextFilename)) {
    if ( file_exists(GCMfilename)) {
        cerr << "Found a GCM file." << endl;
        return true;
    } else {
        cerr << "Did not find a GCM file." << endl;
        return false;
    }
}


bool dbReady (const string& dbname) {
    string dictfilename = dbname + DICT_TAG;
    string datadirname = dbname + DBDIR_TAG;
    string dbinfofilename = dbname + DBINFO_TAG;

    if ( file_exists(dbinfofilename) && file_exists(dictfilename) && dir_exists(datadirname)) {
        return true;
    } else {
      cerr << "Database is not ready. These file and directories are missing: "  << dictfilename << " & " << datadirname << endl;
        return false;
    }
}

int makeDir (const string& dirName) {
    if (mkdir(dirName, 0755) >= 0) {
        cerr << "Created directory " << dirName << endl;
        return 0;
    } 
    else {
        ostringstream message;
        message << "Could not create directory " << dirName << " : " << strerror(errno) << " .. Aborting! " << endl;
        throw Exception(message.str());
    }
}

// void removeDBFiles (const string& dbname, const string& dbpath) {
//     string dbBase = dbpath + dbname;
//     string dictfilename = dbpath + dbname + ".dict";
//     bool errorState = false;
//     unsigned long fileCount = 0;
//     string datadirname = dbpath + dbname + ".data";
    
//     cerr << "Removing database now....Please be patient."<< endl;
//     fs::path dbFile( dbBase);
//     fs::path dictLoc( dictfilename);
//     fs::path dataLoc( datadirname);

//     if (fs::remove(dbFile)) {
//         cerr << "Removed " << dbFile.string() << endl;
//     } else {
//         cerr << "Could not remove " << dbFile.string() << ". File does not exist" << endl;
//         errorState = true;
//     }

//     if (fs::remove(dictLoc)) {
//         cerr << "Removed " << dictLoc.string() << endl;
//     } else {
//         cerr << "Could not remove " << dictLoc.string() << ". File does not exist" << endl;
//         errorState = true;
//     }

//     fileCount = fs::remove_all(dataLoc);
//     if (fileCount > 0) {
//         cerr << "Removed " << dataLoc.string() << ", " << fileCount << " files." << endl;
//     } else {
//         cerr << "Could not remove " << dataLoc.string() << ". Directory does not exist" << endl;
//         errorState = true;
//     }
    
//     if (errorState) {
//         throw Exception("There were errors detected. Database was not able to be removed.");
//     }
    
// }

// bool dbExists (const string& dbname) {
//     bool test1, test2, test3;
//     string dictfilename = dbname + ".dict";
//     string datadirname = dbname + ".data";

// //     fs::path dbFile( fs::initial_path() / dbname);
// //     fs::path dictFile( fs::initial_path() / dictfilename);
// //     fs::path dataLoc( fs::initial_path() / "data");

//     test1=false;
//     test2=false;
//     test3=false;

//     try {
//         fs::path dbFile(dbname); 
//         test1 = fs::exists( dbFile ); 
//     }
//     catch(std::exception const& e) {}
//     catch(...) { 
//         std::cout << "whoops!" << std::endl; 
//     }
//     try {
//         fs::path dataLoc(datadirname);
//         test2 = fs::exists( dataLoc );
//     }
//     catch(std::exception const& e) {}
//     catch(...) { 
//         std::cout << "whoops!" << std::endl; 
//     }
//     try {
//         fs::path dictFile(dictfilename);
//         test3 = fs::exists(dictFile);
//     }
//     catch(std::exception const& e) {}
//     catch(...) { 
//         std::cout << "whoops!" << std::endl; 
//     }
//     if (test1 || test2 || test3) {
//         return true;
//     } else {
//          return false;
//     }
// }

// void makeDir (const string& dirName) {
 
//   //    fs::path dirLoc( fs::initial_path() / dirName);
//   fs::path dirLoc( dirName);
//     if (fs::create_directory (dirLoc)) {
//         cerr << "Created directory " << dirLoc.string() << endl;
//     } else {
//         cerr << "Did not create directory " << dirLoc.string() << " as it already existed." << endl;
//     }
// }


void split_words(string text, vector<string> &words, char ws)
{
    if (text.length() > 0 && text[0] == ws) {
        text.erase(0, text.find_first_not_of(ws));
    }
    while (text.length() > 0) {
        words.push_back(text.substr(0, text.find_first_of(ws)));
        text.erase(0, text.find_first_of(ws));
        text.erase(0, text.find_first_not_of(ws));
    }
}


void ps (string& s){ cout << s << endl; }

void pd (Dictionary& D) {
    for (Dictionary::iterator i = D.begin(); i != D.end(); i++) {
        cout << i->first << ":" << i->second << " , ";
    }
    cout << endl;
}


vector<int> createWeightScheme(int windowLen, int behind, int scheme)
{
    int ahead = windowLen - behind;
    assert(windowLen >= 0);
    vector<int> weights(windowLen);
    //    vector<int> error(1);
    //    error[0]=0;
    if(scheme == SDDB::RAMPED_LINEAR)
    {
        cerr << "Creating ramped, linear weighting scheme\n";
        for(int i = 0; i < behind; i++)
            weights[i] = i + 1;
        for(int i = behind; i < windowLen; i++)
            weights[i] = windowLen - i;
    }
    else if(scheme == SDDB::RAMPED_QUADRATIC)
    {
        cerr << "Creating ramped, quadratic weighting scheme\n";
        for(int i = 0; i < behind; i++)
            weights[i] = (i + 1)*(i + 1);
        for(int i = behind; i < windowLen; i++)
            weights[i] = (windowLen - i)*(windowLen - i);
    }
    else if (scheme == SDDB::FORWARD_RAMP)
    {
        cerr << "Creating forward ramp, flat backwards\n";
        for(int i = 0; i < behind; i++)
            weights[i] = 1;
        for(int i = behind; i < windowLen; i++)
            weights[i] = windowLen - i;
    }
    else if (scheme == SDDB::BACKWARD_RAMP)
    {
        cerr << "Creating backwards ramp, flat forwards\n";
        for(int i = 0; i < behind; i++)
            weights[i] = i + 1;
        for(int i = behind; i < windowLen; i++)
            weights[i] = 1;
    }
    else if (scheme == SDDB::INVERSE_RAMP)
    {
        cerr << "Creating inverse ramped linear weighting scheme\n";
        for(int i = 0; i < behind; i++)
            weights[i] = behind - i;
        for(int i = behind; i < windowLen; i++)
            weights[i] = ((i - behind) + 1);
    }
    else if (scheme == SDDB::INVERSE_QUADRATIC)
    {
        cerr << "Creating inverse quadratic weighting scheme\n";
        for(int i = 0; i < behind; i++)
            weights[i] = (behind - i) * (behind - i) ;
        for(int i = behind; i < windowLen; i++)
            weights[i] = ((i - behind) + 1) * ((i - behind) + 1);
    }
    else if (scheme == SDDB::SECOND_WORD)
    {
        cerr << "Creating second-word heavy weighting scheme\n";
        for(int i = 0; i < windowLen; i++)
            weights[i] = 1;
        if (ahead > 2)
            weights[windowLen - ahead + 1] = 10;
        if (behind > 2)
            weights[windowLen - ahead - 2] = 10;
    }
    else if (scheme == SDDB::THIRD_WORD)
    { 
        cerr << "Creating third-word heavy weighting scheme\n";
        for(int i = 0; i < windowLen; i++)
            weights[i] = 1;
        if (ahead > 3)
            weights[windowLen - ahead + 2] = 10;
        if (behind > 3)
            weights[windowLen - ahead - 3] = 10;
    }
    else if (scheme == SDDB::FOURTH_WORD)
    { 
        cerr << "Creating third-word heavy weighting scheme\n";
        for(int i = 0; i < windowLen; i++)
            weights[i] = 1;
        if (ahead > 4)
            weights[windowLen - ahead + 3] = 10;
        if (behind > 4)
            weights[windowLen - ahead - 4] = 10;
    }
    else {
        cerr << "Creating flat weighting scheme\n";
        for(int i = 0; i < windowLen; i++)
            weights[i] = 1;
    } 
    cerr << "List of weights is: ";
    for(int i = 0; i < windowLen; i++)
        cerr << weights[i] << " ";
    cerr << '\n';
    return weights;
}


// for PREPROCCESSING 

bool
//
// tells whether a character is an acceptable
// component of a word or not
//
acceptable_word_char(char ch, bool alphaOnly)
{
  // spaces separate words
  if(isspace(ch))
    return false;
  else if(!alphaOnly)
    return true;
  else
    return isalpha(ch);
}

string
//
// creates a duplicate of the word, with all or trailing/leading
// non-alphanumeric characters removed.
//
strip_non_alpha(string word, bool sidesOnly)
{
  string stripped = "";
  //cout << "Stripping word: " << word << endl; 
  if(sidesOnly) {
    unsigned int i, j;
    // find the starting position...
    for(i = 0; i < word.length() - 1; i++) {
      if(isalpha(word[i])) {
	break; }
    }
    // find the ending position...
    for(j = word.length() - 1; j > i; j--) {
      if(isalpha(word[j])) {
	break; }
    }
    // copy!
    stripped = word.substr(i,(j-i+1));
    //cout << "Stripped word: " << stripped << endl; 
    return stripped;
  }
  else {
    for(unsigned int i = 0; i < word.length(); i++) {
      if(!isalpha(word[i])) {     
	word.erase(i,1); // remove char
	i = i - 1;
      }
    }
    //cout << "Stripped word: " << word << endl; 
    return word;
  }
}

wordpair
// 
// Extracts the next word from the input stream. This function has the same
// functionality as "cin >> word", excepts it additionally makes sure that
// the word buffer is not over-run! And it separates possessive endings.
//
ExtractWord(string localword, const bool normCase, const bool englishContractions, string lang)
{
  wordpair output;
  output.main = "";
  output.possessive = ""; 
  // this is mainly used for holding "'S" since we can't put strings back

  // Normalize the case.
  if (normCase) {
    localword = downstring(localword,lang);
  }
  
  //check if the special posessive words are there.

  //NOTE: english Contractions currently only works when normCase is TRUE.
  if (englishContractions) {
    if (!((localword == "he's") || (localword =="she's") || (localword == "it's") || (localword == "here's") || (localword == "there's") || (localword == "what's") || (localword == "let's") || (localword == "that's") || (localword == "where's") || (localword == "who's"))) 
      {
	//cout << "Looking for possessive edning for: " << localword << endl;
	// we have to check if our word ends in "'S" and it is longer than
	// "'S" alone. If that's the case, we throw back the "'S"
	string::size_type pos;
	pos = localword.find(POSSESSIVESUFFIX,0);
	if ((localword.length() > 2) && (pos != string::npos)) {
	  output.possessive = POSSESSIVESUFFIX;
	  localword.erase(pos,2);
	  //cout << "Found possessive ending for: " << localword << endl;
	}
      }
    //  cout << "adding word" << localword << endl; 
  }
  output.main = localword;
  return output;
}

Numpair
CleanWord(string word, Dictionary& dict,const bool normCase, const bool englishContractions, string lang)
{
  wordpair extractedWords;
  Numpair cleanedPair;
  cleanedPair.first = 0;
  cleanedPair.second = 0;

  if (word.length() < MAX_WORDLEN) { 
    extractedWords = ExtractWord(word,normCase,englishContractions,lang);
    if (extractedWords.possessive != "") {
      cleanedPair.second= dict[extractedWords.possessive];
    }
    // If word is not in lexicon
    if (dict.find(extractedWords.main) == dict.end()) {
      // If word is dirty Try stripping non-letters off the sides of the word (like 500Club -> Club)
      //	cout << word << endl;
      string cleanWord = strip_non_alpha(extractedWords.main, true);
      //cout << "Tried cleaning into: " << cleanWord << endl;
      if(cleanWord.empty()) {
	// we encountered no letters ... skip this one
	return(cleanedPair);
      }
      else  {
	// if clean word is in lexicon, add it.
	if (dict.find(cleanWord) != dict.end()) {
	  cleanedPair.first = dict[cleanWord];
	  return(cleanedPair);
	}
	else {
	  // Strip all non letters out.
	  string newWord = strip_non_alpha(extractedWords.main, false);
	  if (dict.find(newWord) != dict.end()) {
	    cleanedPair.first = dict[newWord];
	    return(cleanedPair);
	  }
	  else {
	    return(cleanedPair);
	  }
	}
      }
    } else {
      cleanedPair.first = dict[extractedWords.main];
      return(cleanedPair);
    }
  }
  else {
    return(cleanedPair);
  }
}


void readMetaData(string dbBase, int& numVectors, int& vectorLen, int& windowLenBehind, int& windowLenAhead, long& corpussize) 
{
    ifstream in;

    cerr << "Using database : " << dbBase << endl;
    in.open(dbBase.c_str());
    if (!in.good()) {
        ostringstream message;
        message << "Could not open database MetaData file " << dbBase << ". Aborting." << endl;
        throw Exception(message.str());
    }

    try {
        in >> numVectors;
        in >> vectorLen;
        in >> windowLenBehind;
        in >> windowLenAhead;
        in >> corpussize;
    }
    
    catch(...) { 
        cerr << "Meta data was missformated. Aborting." << endl;
    }

    in.close();
}

void writeMetaData(string dbBase, int numVectors, int vectorLen, int windowLenBehind, int windowLenAhead, long corpussize) 
{
    ofstream out;

    cerr << "Using database : " << dbBase << endl;
    out.open(dbBase.c_str());
    if (!out.good()) {
        ostringstream message;
        message << "Could not write to database MetaData file " << dbBase << ". Aborting." << endl;
        throw Exception(message.str());
    }

    try {
        out << numVectors << " ";
        out << vectorLen<< " ";
        out << windowLenBehind << " ";
        out << windowLenAhead << " ";
        out << corpussize << " ";
    }
    
    catch(...) { 
        cerr << "MetaData could not be written . Aborting." << endl;
    }

    out.close();
}

string getLang() {
  char* langVar;
  langVar = getenv("LANG");
  if (langVar!=NULL) {
    string output(langVar,2);
    cerr << "Language is set to: " << output << endl; 
    return(output);
  } else {
    return("en");
  }
}

string downstring(string localword, string lang) {
  // old Way to do it, not unicode aware.....
  //
  //  for (unsigned int j=0; j < localword.length(); ++j)    {
  //    localword[j]=toupper(localword[j]);   
  //  }
  //  const uint8_t * word = static_cast<const uint8_t*>(localword.c_str());


  // New way to do it using libunicode
  //
  //Get string length
  size_t length = localword.size();
  // create correct type for c-style unicode string
  const uint8_t * word = (const uint8_t*)localword.c_str();
  // create output buffer
  uint8_t output[200];
  uint8_t * errCode;
  uint8_t val;
  errCode = &val;
  // create output length location
  size_t outLength = 200;
  // make lowercase, normalize and put output in the output buffer, length in the outLength variable
  if (u8_check(word, length)) {
    cerr << "Problem processing the word: "<< word << endl;
    throw Exception("This is an invalid UTF8 in string. Please make sure that you are using UTF8 encoding in all input files. Exiting.");
  }
  if (!u8_tolower(word, length, lang.c_str(), UNINORM_NFKD, output, &outLength))  {
    cerr << word << endl;
    throw Exception("Error during case conversion (in downstring) ");
    //    return(" ");
  }
  // return a c++ string, using begining and end pointers to the c-style string!
  return(string((const char *)output,(const char *)output+outLength));
}
