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

#ifndef EXCEPTION_H
#define EXCEPTION_H

#include <exception>
#include <string>

/** Base class for exceptions. */
class Exception
    : public std::exception
{
public:
    /** Construct with no message. */
    Exception();

    /** Construct with message. */
    Exception(const std::string& message);
    
    /** Destructor. */
    virtual ~Exception() throw();

    /** Implementation of std::exception::what(). */
    const char* what() const throw();

private:
    std::string _message;
};

#endif // EXCEPTION_H
