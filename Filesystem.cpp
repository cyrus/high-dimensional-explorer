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

#define _LARGE_FILES

/* 

This is a collection of functions that handle filesystem operations needed by
other parts of HiDEx. 

*/
#include "sys/stat.h"
#include "Filesystem.h"

using namespace std;

bool file_exists(const string& fname)
{
    struct stat sbuf;
    // exists && is a regular file
    return stat(fname, &sbuf) == 0 && S_ISREG(sbuf.st_mode);
}

bool dir_exists(const string& fname)
{
    struct stat sbuf;
    // exists && is a dir
    return stat(fname, &sbuf) == 0 && S_ISDIR(sbuf.st_mode);
}

int rmdir(const string& filename)
{
    // Check filename exists and is actually a directory
    struct stat sb;
    if (stat(filename, &sb) != 0 || !S_ISDIR(sb.st_mode)) {
        cerr << filename << " does not exist or is not a directory" << endl;
        return -1;
    }

    string safefile = filename;
    string::size_type p = 0;
    while (p < safefile.size()) {
        // Don't escape a few safe characters which are common in filenames
        if (!isalnum(safefile[p]) && strchr("/._-", safefile[p]) == NULL) {
            safefile.insert(p, "\\");
            ++p;
        }
        ++p;
    }
    system("rm -rf " + safefile);
    return 0;
}

int mkdir(const string& filename, mode_t mode) {
    return mkdir(filename.c_str(), mode);
}

int stat(const string& filename, struct stat *buf) {
    return stat(filename.c_str(), buf);
}

void touch(const string& filename) {
   int fd = open(filename.c_str(), O_CREAT|O_WRONLY, 0644);
   if (fd >= 0) close(fd);
}

int unlink(const string& filename) { return unlink(filename.c_str()); }

int system(const string& command) { return system(command.c_str()); }


