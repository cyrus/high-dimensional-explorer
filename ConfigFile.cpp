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


/** @file ConfigFile.cpp
    Class ConfigFile.
    This Class is used to process configuration settings when the program begins
    execution. It can handle text, integer and floating point configuration settings.
*/



#define _LARGE_FILES
#include "ConfigFile.h"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include "Exception.h"

using namespace std;

namespace {

/** Check if non-whitespace characters are left in stream. */
bool checkNothingRemains(istream& in)
{
    string remains;
    in >> remains;
    return (in && remains != "");
}

/** Trim leading and trailing whitespaces from string. */
void trim(string& str)
{
    char const* whiteSpace = " \t\r";
    size_t pos = str.find_first_not_of(whiteSpace);
    str.erase(0, pos);
    pos = str.find_last_not_of(whiteSpace); 
    str.erase(pos + 1);
}

} // namespace

ConfigFile::ConfigFile(const string& filename,
                       const vector<string>& allowedKeys, string& multi)
    : _filename(filename),
      _allowedKeys(allowedKeys)
{
    if (multi.size() > 0) {
        readstring(multi);
    } else {
        fstream in(filename.c_str());
        if (! in.is_open())
            throw Exception("Could not open " + filename);
        read(in);
    }
}

ConfigFile::ConfigFile(istream& in, const vector<string>& allowedKeys)
    : _allowedKeys(allowedKeys)
{
    read(in);
}


bool ConfigFile::contains(const std::string& key) const
{
    assert(isAllowedKey(key));
    return (_map.find(key) != _map.end());
}

string ConfigFile::get(const string& key) const
{
    assert(isAllowedKey(key));
    map<string,string>::const_iterator pos = _map.find(key);
    if (pos == _map.end())
        throw Exception(_filename + ": Missing key: " + key );
    return pos->second;
}

string ConfigFile::get(const string& key, const string& defaultValue) const
{
    if (! contains(key))
        return defaultValue;
    return get(key);
}

bool ConfigFile::getBool(const string& key) const
{
    string value = get(key);
    istringstream in(value);
    int result;
    in >> result;
    if (! in || checkNothingRemains(in) || (result != 0 && result != 1))
        throw Exception("Expected boolean value (0/1) for " + key);
    return result;
}

bool ConfigFile::getBool(const string& key, bool defaultValue) const
{
    if (! contains(key))
        return defaultValue;
    return getBool(key);
}

double ConfigFile::getDouble(const string& key) const
{
    string value = get(key);
    istringstream in(value);
    double result;
    in >> result;
    if (! in || checkNothingRemains(in))
        throw Exception("Expected double value for " + key);
    return result;
}

double ConfigFile::getDouble(const string& key, double defaultValue) const
{
    if (! contains(key))
        return defaultValue;
    return getDouble(key);
}

int ConfigFile::getInt(const string& key) const
{
    string value = get(key);
    istringstream in(value);
    int result;
    in >> result;
    if (! in || checkNothingRemains(in))
        throw Exception("Expected integer value for " + key);
    return result;
}

int ConfigFile::getInt(const string& key, int defaultValue) const
{
    if (! contains(key))
        return defaultValue;
    return getInt(key);
}

int ConfigFile::getPositiveInt(const string& key, int defaultValue) const
{
    int value = getInt(key, defaultValue);
    if (value < 0)
        throw Exception("Value for " + key + " must be positive");
    return value;
}

size_t ConfigFile::getSizet(const string& key) const
{
    string value = get(key);
    istringstream in(value);
    size_t result;
    in >> result;
    if (! in || checkNothingRemains(in))
        throw Exception("Expected size_t value for " + key);
    return result;
}

size_t ConfigFile::getSizet(const string& key, size_t defaultValue) const
{
    if (! contains(key))
        return defaultValue;
    return getSizet(key);
}

unsigned long long int ConfigFile::getLong(const string& key) const
{
    string value = get(key);
    istringstream in(value);
    long long int result;
    in >> result;
    if (! in || checkNothingRemains(in))
        throw Exception("Expected long int value for " + key);
    return result;
}

unsigned long long int ConfigFile::getLong(const string& key, unsigned long long int defaultValue) const
{
    if (! contains(key))
        return defaultValue;
    return getLong(key);
}

void ConfigFile::handleLine(const string& line, int lineNumber)
{
    size_t pos = line.find('=');
    if (pos == string::npos)
        throwError("Missing =", lineNumber);
    string key = line.substr(0, pos);
    trim(key);
    if (! isAllowedKey(key))
        throwError("Unknown key: " + key, lineNumber);
    if (contains(key))
        throwError("Duplicate entry: " + key, lineNumber);
    string value = line.substr(pos + 1);
    trim(value);
    _map[key] = value;
}

bool ConfigFile::isAllowedKey(const string& key) const
{
    return (find(_allowedKeys.begin(), _allowedKeys.end(), key)
            != _allowedKeys.end());
}

void ConfigFile::read(istream& in)
{
    _map.clear();
    int lineNumber = 0;
    while (! in.eof())
    {
        string line;
        getline(in, line);
        trim(line);
        ++lineNumber;
        if (line == "" || line[0] == '#')
            continue;
        handleLine(line, lineNumber);        
    }
}

void ConfigFile::readstring(string& inputstring)
{
    _map.clear();
    int itemNumber = 0;
    while (! inputstring.empty())
    {
        string line;
        line = getItem(inputstring);
        trim(line);
        ++itemNumber;
        if (line == "" || line[0] == '#')
            continue;
        handleLine(line, itemNumber);        
    }
}

string ConfigFile::getItem(string& inputstring)
{
    size_t pos = inputstring.find(',');
    if (pos == string::npos) {
        string item = "";
        inputstring = item;
        return(item);
    }
    string item = inputstring.substr(0, pos);
    string tempstring = inputstring;
    inputstring.erase();
    inputstring = tempstring.substr(pos+1);
    trim(item);
    return(item);
}
void ConfigFile::throwError(const string& message, int lineNumber) const
{
    ostringstream buffer;
    if (_filename != "")
        buffer << _filename << ": ";
    buffer << "Line " << lineNumber << ": " << message;
    throw Exception(buffer.str());
}
