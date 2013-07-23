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

#ifndef Filesystem_H
#define Filesystem_H

#include "SDDB.h"
#include <unistd.h>
#include <ctype.h>
#include <fcntl.h>

bool file_exists(const string &fname);

bool dir_exists(const string &fname);

int rmdir(const string &filename);

int mkdir(const string &filename, mode_t mode);

int stat(const string &filename, struct stat *buf);

void touch(const string &filename);

int unlink(const string &filename);

int system(const string &command);

#endif
