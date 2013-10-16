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

#ifndef utilities_H
#define utilities_H

#include "SDDB.h"

#define POSSESSIVESUFFIX       "'s"
#define MAX_WORDLEN                      30

struct word_distance_pair
{
    int id;          // id code
    Float distance;    // distance
};

//struct wordpair
//{
//    string main;    // main word
//    string possessive;     // posessive
//};

vector<int> createWeightScheme(int windowLen, int behind, int scheme);

void LoadPairs(istream& in, vector<pairdata> &results, const bool normCase);

void LoadWords(istream& in, const int wordlistsize, vector<resultdata> &results,const bool normCase);

void addtoresults(vector<resultdata> &results, string word, Float ARC, Float InverseNcount);

std::string getcorrelation(vector<resultdata> &results);

bool getWord(istream& in, string& word);

template<class T> int array_place(T *array, T elem, int len);

template<class T> void binary_sort(T array[], int len);

void write_dict_and_freqs(Dictionary& D, string filename, FrequencyMap& frequencies, const string eod);

void write_vars(Dictionary& D, string filename, VarianceMap& frequencies, const string eod);

int makeDir (const string& dirName); 

void wait ( size_t seconds );

// ****************************************************************************
// Private functions
//
//
// ****************************************************************************


std::string timestamp ();

std::string itoa(int value, unsigned int base);

void build_starting_dict(Dictionary& D, string filename,  FrequencyMap& frequencies, const bool normCase);

void build_dict_and_freqs(Dictionary& D, string filename, FrequencyMap& frequencies, const string eod);

void build_idMap(Dictionary& D, idMap& words);
//void value_sort(const DictPair *pairs[], const int size);

Float get_magnitude(Float *vect, int context_size);

void add_smallest(word_distance_pair *entries, int len, int id, Float distance);

bool dbExists(const string& dbname);

bool dbReady (const string& dbname, const string& dbpath );

void removeDBFiles (const string& dbname, const string& dbpath);

bool acceptable_word_char(char ch, bool alphaOnly);

string strip_non_alpha(string word, bool sidesOnly);

wordpair ExtractWord(string localword, const bool normCase, const bool englishContractions, string lang);

Numpair CleanWord(string word, Dictionary& D, const bool normCase, const bool englishContractions, string lang);

void readMetaData(string dbBase, int& numVectors, int& vectorLen, int& windowLenBehind, int& windowLenAhead, long& corpussize);

void writeMetaData(string dbBase, int numVectors, int vectorLen, int windowLenBehind, int windowLenAhead, long corpussize);  
//Float CalcSimilarity(const vector<Float*> &vectors, size_t num_dimensions, int w1, int w2, const string algorithm); 

bool GCMexists(const string& dbname);

string downstring(string localword, string lang);

string getLang();

#endif // utilities_H
