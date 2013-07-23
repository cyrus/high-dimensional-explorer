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

/*
This is a helper class to process textual messages in thrown exceptions that pass
up the chain to the top level.
*/


#define _LARGE_FILES

#include "Exception.h"

using namespace std;

Exception::Exception()
{
}

Exception::Exception(const string& message)
{
    _message = message;
}

Exception::~Exception() throw()
{
}

const char* Exception::what() const throw()
{
    return _message.c_str();
}
