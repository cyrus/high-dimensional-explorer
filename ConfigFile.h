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

/** @file ConfigFile.h
    Class ConfigFile.
*/

#ifndef CONFIGFILE_H
#define CONFIGFILE_H

#include <iosfwd>
#include <map>
#include <string>
#include <vector>

/** Parser for config files.
    The format is similar to Windows or KDE INI-files.
    - Comment lines begin with \#
    - Empty lines are ignored
    - Keys are case-sensitive

    Key value pairs are defined with lines like:
    @verbatim key = value @endverbatim
    The value can contain spaces; leading and trailing spaces are stripped.
*/
class ConfigFile
{
public:
    /** Create a configuration from a file.
        @throws Exception If file contains syntax errors or not allowed keys.
    */
    ConfigFile(const std::string& filename,
               const std::vector<std::string>& allowedKeys, std::string& multi);

    /** Create a configuration from a stream.
        @throws Exception If file contains syntax errors or not allowed keys.
    */
    ConfigFile(std::istream& in, const std::vector<std::string>& allowedKeys);


    /** Check if file contains a key. */
    bool contains(const std::string& key) const;

    /** Get value.
        @throws Exception If file does not contain the key.
    */
    std::string get(const std::string& key) const;

    /** Get value or return default.
        @return The value from the file or the default value, if the file
        does not contain the key.
    */
    std::string get(const std::string& key,
                    const std::string& defaultValue) const;

    /** Get boolean value.
        @throws Exception If file does not contain the key or value is not a
        boolean (0 or 1).
    */
    bool getBool(const std::string& key) const;

    /** Get boolean value or return default.
        @return The value from the file or the default value, if the file
        does not contain the key.
        @throws Exception If value is not a boolean (0 or 1).
    */
    bool getBool(const std::string& key, bool defaultValue) const;

    /** Get double value.
        @throws Exception If file does not contain the key or value is not a
        double.
    */
    double getDouble(const std::string& key) const;

    /** Get double value or return default.
        @return The value from the file or the default value, if the file
        does not contain the key.
        @throws Exception If value is not a double.
    */
    double getDouble(const std::string& key, double defaultValue) const;

    /** Get integer value.
        @throws Exception If file does not contain the key or value is not an
        integer.
    */
    int getInt(const std::string& key) const;

    /** Get integer value or return default.
        @return The value from the file or the default value, if the file
        does not contain the key.
        @throws Exception If value is not an integer.
    */
    int getInt(const std::string& key, int defaultValue) const;

    /** Get positive integer value or return default.
        @return The value from the file or the default value, if the file
        does not contain the key.
        @throws Exception If value is not a positive integer or 0.
    */
    int getPositiveInt(const std::string& key, int defaultValue) const;

    size_t getSizet(const std::string& key) const;

    /** Get size_t value or return default.
        @return The value from the file or the default value, if the file
        does not contain the key.
        @throws Exception If value is not an integer.
    */
    size_t getSizet(const std::string& key, size_t defaultValue) const;

    /** Get size_t value or return default.
        @return The value from the file or the default value, if the file
        does not contain the key.
        @throws Exception If value is not a positive integer or 0.
    */

    unsigned long long int getLong(const std::string& key) const;

    /** Get size_t value or return default.
        @return The value from the file or the default value, if the file
        does not contain the key.
        @throws Exception If value is not an integer.
    */
   unsigned long long int getLong(const std::string& key, unsigned long long int defaultValue) const;

private:
    const std::string _filename;

    const std::vector<std::string> _allowedKeys;

    std::map<std::string,std::string> _map;

    void handleLine(const std::string& line, int lineNumber);

    std::string getItem(std::string& inputstring);

    bool isAllowedKey(const std::string& key) const;

    void read(std::istream& in);

    void readstring(std::string& inputstring);

    void throwError(const std::string& message, int lineNumber) const;
};

#endif // CONFIGFILE_H
