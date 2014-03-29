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

//****************************************************************************
//
// hidex.cpp
//
// This is where the master control program - getting options, reading in new corpora, building new
// databases, calculating measures based on word vectors.
// 
// Based on the software called "sddbmaker" written  by Geoff Hollis.
//
// Currently under continuous development by Cyrus Shaoul. 
// Renamed from "sddbmaker" to "hidex (High Dimensional EXplorer)"
// on Aug 9th, 2005.
//
// For instructions on how to use HiDEx, please see the documentation which can be found in the file
//                HowToUseHiDEx.txt
// in the same directory as this source code.
// 
//****************************************************************************

#define _LARGE_FILES
#include <new>
#include <cstdlib>
#include <iostream>
#include <exception>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>
#include "ConfigFile.h"
#include "Exception.h"
#include "Filesystem.h"
#include "SDDB.h"
#include "time.h"

using namespace std;

namespace {
  
Settings settings;

void readSettings(const string& configString, string& multi)
{
    vector<string> allowedKeys;
    allowedKeys.push_back("contextSize");
    allowedKeys.push_back("corpusFilename");
    allowedKeys.push_back("dbname");
    allowedKeys.push_back("dictFilename");
    allowedKeys.push_back("eod");
    allowedKeys.push_back("maxWindowAhead");
    allowedKeys.push_back("maxWindowBehind");
    allowedKeys.push_back("neighbourhoodSize");
    allowedKeys.push_back("percenttosample");
    allowedKeys.push_back("separate");
    allowedKeys.push_back("stepsize");
    allowedKeys.push_back("useThreshold");
    allowedKeys.push_back("weightingScheme");
    allowedKeys.push_back("metric");
    allowedKeys.push_back("saveGCM");
    allowedKeys.push_back("normalization");
    allowedKeys.push_back("windowLenAhead");
    allowedKeys.push_back("windowLenBehind");
    allowedKeys.push_back("wordlistfilename");
    allowedKeys.push_back("wordlistsize");
    allowedKeys.push_back("outputpath");
    allowedKeys.push_back("dbpath");
    allowedKeys.push_back("multipleFiles");
    allowedKeys.push_back("normCase");
    allowedKeys.push_back("englishContractions");
    allowedKeys.push_back("useVariance");
    ConfigFile cfg(configString, allowedKeys, multi);
    settings.dbname = cfg.get("dbname");
    settings.dictFilename = cfg.get("dictFilename");
    settings.contextSize = cfg.getPositiveInt("contextSize", 9000);
    settings.corpusFilename = cfg.get("corpusFilename", "");
    settings.eod = cfg.get("eod", "---END.OF.DOCUMENT---");
    settings.maxWindowAhead = cfg.getPositiveInt("maxWindowAhead", 10);
    settings.maxWindowBehind = cfg.getPositiveInt("maxWindowBehind", 10);
    settings.neighbourhoodSize = cfg.getPositiveInt("neighbourhoodSize", 20);
    settings.percenttosample = cfg.getDouble("percenttosample", 0.05);
    settings.separate = cfg.getBool("separate", false);
    settings.stepsize = cfg.getPositiveInt("stepsize", 800);
    settings.usezscore = cfg.getBool("useThreshold", false);
    settings.weightingScheme = cfg.getPositiveInt("weightingScheme", 1);
    settings.metric = cfg.get("metric", "");
    settings.saveGCM = cfg.getBool("saveGCM", "");
    settings.normalization = cfg.get("normalization", "");
    settings.windowLenAhead = cfg.getPositiveInt("windowLenAhead", 10);
    settings.windowLenBehind = cfg.getPositiveInt("windowLenBehind", 10);
    settings.wordlistfilename = cfg.get("wordlistfilename", "");
    settings.outputpath = cfg.get("outputpath", "");
    settings.dbpath = cfg.get("dbpath", "");
    settings.wordlistsize = cfg.getPositiveInt("wordlistsize", 300000000);
    settings.multipleFiles = cfg.get("multipleFiles", "");
    settings.normCase = cfg.getBool("normCase", true);
    settings.englishContractions = cfg.getBool("englishContractions", true);
    settings.useVariance = cfg.getBool("useVariance", false);
}


/** Create a new SDDB. */
void create()
{
    cerr << "Database Name = " << settings.dbname << '\n'
         << "Dictionary File Name = " << settings.dictFilename << '\n'
         << "Max Window length Behind = " << settings.maxWindowBehind << '\n'
         << "Max Window length Ahead = " << settings.maxWindowAhead << '\n'
         << "Case Normalization = " << settings.normCase << '\n'
         << "English Contraction Normalization = " << settings.englishContractions << '\n'
         << "EOD = " << settings.eod << '\n';
    if(settings.dbpath != "") {
      cerr << "\nSetting Alternate Working Directory Path = " << settings.dbpath << '\n';
      if (!dir_exists(settings.dbpath)) {
        ostringstream buffer;
        buffer << "  A working directory called " << settings.dbpath << " defined in your configuration file does not exist. Please change settings and try again. Exiting." ;
        throw Exception(buffer.str());
      }   
    }
    SDDB db(settings.dbname, settings.dbpath);
    db.setOptions(settings);
    db.initialize(settings.dictFilename, settings.maxWindowBehind, settings.maxWindowAhead, settings.eod);
}


/** Update the co-occurance counts of an SDDB with data a corpus files.
*/
void update()
{
    cerr << "Entering update processing....\n";
    if (settings.corpusFilename == "")
        throw Exception("  No Corpus file specified in config file");
    if (settings.stepsize < 1000)
        throw Exception("  A minimum stepsize of 1000 is needed");
    if (settings.dbpath != "") {
      if (!dir_exists(settings.dbpath)) {
        ostringstream buffer;
        buffer << "A working directory called " << settings.dbpath << " does not exist. Cannot update database. Exiting." ;
        throw Exception(buffer.str());
      }
    }   
    // Set up corpus filename.
    string filename = settings.dbpath + settings.corpusFilename;
    
    cerr << "stepsize = " << settings.stepsize << '\n'
         << "EOD = " << settings.eod << '\n'
         << "corpus = " << filename << '\n'
         << "db = " << settings.dbname << '\n'
         << "Case Normalization = " << settings.normCase << '\n'
         << "Max Memory = " << settings.maxMemory << '\n'
         << "Calculate Word Vector Variances = " << settings.useVariance << '\n'
         << "English Contraction Normalization = " << settings.englishContractions << '\n';
    SDDB db(settings.dbname, settings.dbpath);
    db.load(settings.eod, settings.maxMemory);
    db.setOptions(settings);
    db.setCurrentStep(0);
    
    // reprocess the corpus file in steps until there are no steps left.
    // this is due to the memory limitations in most environments.
    // 
    int stepcount = 1;
    while(true)
    {
	cerr << "Opening corpus file." << endl;
        ifstream in(filename.c_str());
        if (! in.is_open())
            throw Exception("Could not open corpus file: " + filename + " . Exiting.");
        else
            cerr << "Opened corpus file. Processing corpus..." << endl;
        db.update(in,0);
	cerr << "Closing corpus file." << endl;
        in.close();
        cerr << "Flushing data to disk." << endl;
        db.flushDB();
        if(! db.stepUp())
            break;
        cerr << "Finished step number " << stepcount << "\n";
        stepcount++;
    }
    cerr << "Finished updating database. Closing database files.\n";
    db.close();
}


/** Print out neighborhoods for a set of specified words. */
void printSDs()
{
    int errorcode;
    if (settings.dbpath != "") {
      if (!dir_exists(settings.dbpath)) {
        ostringstream buffer;
        buffer << "A directory called " << settings.dbpath << " does not exist. Cannot read database. Exiting." ;
        throw Exception(buffer.str());
      }
    }
    if (settings.wordlistfilename == "")
        throw Exception("No wordlistfilename specified in configfile. Exiting.");
    string wordlistfilename = settings.dbpath + settings.wordlistfilename;
    ifstream wordlistfile(wordlistfilename.c_str());
    if (! wordlistfile.is_open())
        throw Exception("  Could not open word list file \""
                        + wordlistfilename + "\"");
    if (settings.dbpath != "")
      cerr << "Working Directory = " << settings.dbpath << "\n";
    ostringstream buffer;
    buffer << "Database Name = " << settings.dbname << "\n"
         << "Wordlist FileName = " << settings.wordlistfilename << '\n'
         << "Context size = " << settings.contextSize << '\n'
         << "Window Length Behind " << settings.windowLenBehind << '\n'
         << "Window Length Ahead = " << settings.windowLenAhead << '\n'
         << "Weighting Scheme = " << settings.weightingScheme << '\n'
         << "Similarity Metric = " << settings.metric << '\n'
         << "Save GCM = " << settings.saveGCM << '\n'
         << "Normalization Method = " << settings.normalization << '\n'
         << "Max Neighbourhood Size = " << settings.neighbourhoodSize << '\n'
         << "Use Zscore Thresholds = " << settings.usezscore << '\n'
         << "Percent to sample for Zscore Thresholds = " << settings.percenttosample << '\n'
         << "Separate = " << settings.separate << '\n'
         << "Case Normalization = " << settings.normCase << '\n'
         << "English Contraction Normalization = " << settings.englishContractions << '\n'
         << "Using Variance for context selection = " << settings.useVariance << '\n'
         << "Word List Size = " << settings.wordlistsize << '\n';
    cerr << buffer.str();
    cerr << "Opening Database...\n";
    SDDB db(settings.dbname, settings.dbpath);
    db.setOptions(settings);
    db.load(settings.eod, settings.maxMemory);    
    cerr << "Starting to calculate neighborhoods ....\n";
    errorcode = db.printSDs(wordlistfile, settings.contextSize, settings.weightingScheme,
			    settings.metric, settings.normalization,
			    settings.windowLenBehind, settings.windowLenAhead,
			    settings.neighbourhoodSize, settings.usezscore,
			    settings.separate,
			    settings.percenttosample, settings.wordlistsize, 
			    settings.outputpath, 
			    settings.saveGCM, 
			    buffer.str()
			    );
    if (errorcode != 0)
        cerr << "ErrorCode was: " << errorcode << endl;
}

/** Print out vectors for a set of specified words. */
void printVectors()
{
    int errorcode;
    if (settings.dbpath != "") {
      if (!dir_exists(settings.dbpath)) {
        ostringstream buffer;
        buffer << "  A directory called " << settings.dbpath << " does not exist. Cannot read database. Exiting." ;
        throw Exception(buffer.str());
      }
    }
    if (settings.wordlistfilename == "")
        throw Exception("  No wordlistfilename specified");
    string wordlistfilename = settings.dbpath + settings.wordlistfilename;
    ifstream wordlistfile(wordlistfilename.c_str());
    if (! wordlistfile.is_open())
        throw Exception("  Could not open word list file \""
                        + wordlistfilename + "\"");
    cerr << "Database Name = " << settings.dbname << "\n"
	 << "Database Path = " << settings.dbpath << "\n"
         << "Wordlist FileName = " << settings.wordlistfilename << '\n'
         << "Context size = " << settings.contextSize << '\n'
         << "Window length behind = " << settings.windowLenBehind << '\n'
         << "Window length ahead = " << settings.windowLenAhead << '\n'
         << "Weighting scheme = " << settings.weightingScheme << '\n'
         << "Save GCM = " << settings.saveGCM << '\n'
         << "Neighbourhood Size = " << settings.neighbourhoodSize << '\n'
         << "separate = " << settings.separate << '\n'
         << "Case Normalization = " << settings.normCase << '\n'
         << "English Contraction Normalization = " << settings.englishContractions << '\n'
         << "wordlistsize = " << settings.wordlistsize << '\n'
         << "Using Variance for context selection = " << settings.useVariance << '\n'
         << "percenttosample = " << settings.percenttosample << '\n';
    cerr << "Opening Database...\n";
    SDDB db(settings.dbname, settings.dbpath);
    db.load(settings.eod, settings.maxMemory);    
    cerr << "Printing Vectors ....\n";
    errorcode = db.printVects(wordlistfile, settings.contextSize, settings.weightingScheme,
			      settings.windowLenBehind, settings.windowLenAhead,
                              settings.wordlistsize, 
                              settings.separate,
                              settings.outputpath,
			      settings.normalization,
			      settings.saveGCM
			      );
    if (errorcode != 0)
        cerr << "ErrorCode was: " << errorcode << endl;
}


/** Print out distances between wordpairs for a  set of word pairs. */
void printPairs()
{
    if (settings.wordlistfilename == "")
        throw Exception("No wordlistfilename specified");
    string wordlistfilename = settings.dbpath + settings.wordlistfilename;
    ifstream wordlistfile(wordlistfilename.c_str());
    if (! wordlistfile.is_open())
        throw Exception("  Could not open word list file \""
                        + wordlistfilename + "\"");
    cerr << "Database Name = " << settings.dbname << "\n"
	 << "Database Path = " << settings.dbpath << "\n"
         << "Wordlist FileName = " << settings.wordlistfilename << '\n'
         << "Context size = " << settings.contextSize << '\n'
         << "Window length Behind = " << settings.windowLenBehind << '\n'
         << "Window length Ahead = " << settings.windowLenAhead << '\n'
         << "Weighting scheme = " << settings.weightingScheme << '\n'
         << "Similarity Metric = " << settings.metric << '\n'
         << "Normalization Method = " << settings.normalization << '\n'
         << "Separate F/B Windows= " << settings.separate << '\n'
         << "Using Variance for context selection = " << settings.useVariance << '\n'
         << "Save GCM = " << settings.saveGCM << '\n';
    cerr << "Opening Database...\n";
    SDDB db(settings.dbname, settings.dbpath);
    db.setOptions(settings);
    db.load(settings.eod, settings.maxMemory);    
    cerr << "Starting to calculate word pair similarity ....\n";
    db.printPairs(wordlistfile, settings.contextSize, settings.weightingScheme,
                  settings.windowLenBehind, settings.windowLenAhead, 
                  settings.separate, settings.outputpath,
		  settings.metric,
		  settings.normalization, 
		  settings.saveGCM
		  );
}

void printmulti(const string& configFilename)
{
    fstream in(configFilename.c_str());
    if (! in.is_open())
        throw Exception("  Could not open " + configFilename);

    vector<string> configurations;
    while (! in.eof())
    {
        string configline;
        getline(in, configline);
        if (configline.size() > 1) 
            configurations.push_back(configline);
    }

    for (vector<string>::iterator fn = configurations.begin(); fn != configurations.end(); fn++) {
        cerr << "Starting proccessing with settings: " << *fn << endl;
        readSettings("", *fn);
        printSDs();
    }
    cerr << "Fishined Multi Batch from file " << configFilename << endl;
}

void PrintMultiPairs(const string& configFilename)
{
    fstream in(configFilename.c_str());
    if (! in.is_open())
        throw Exception("Could not open " + configFilename);

    vector<string> configurations;
    while (! in.eof())
    {
        string configline;
        getline(in, configline);
        if (configline.size() > 1) 
            configurations.push_back(configline);
    }

    for (vector<string>::iterator fn = configurations.begin(); fn != configurations.end(); fn++) {
        cerr << "Starting proccessing with settings: " << *fn << endl;
        readSettings("", *fn);
        printPairs();
    }
    cerr << "Fishined Multi Batch from file " << configFilename << endl;
}

/** Update the co-occurance counts of an SDDB with data a corpus files.*/
void remove()
{
    cout << "Are you sure you want to remove the database called " << settings.dbname << "?" << endl;
    cout << "This operation in not reversible." << endl;
    cout << "If you are sure, enter the letter (y): ";
    string answer;
    getline (cin, answer);
    if (answer == "y") {
        removeAllFiles(settings.dbname, settings.dbpath);
    } else {
        cerr << "Aborted removal." << endl;
    }
}

void usage(ostream& out, const char* programName)
{
    out << "usage: " << programName
        << " [-f <configfile>] -m <mode> \n"  
	<< "     If no configfile is specified, config.txt is used.\n" 
	<< "     Some popular modes are: create, update, remove, getneighbors, getsimilarity, and getvectors"
	<< "     For the list of all possible modes, please read the documentation included.\n";
}

void warranty()
{
  // Rember to change this when changing the version number
  Float versionNumber = 0.08;

  cerr <<  "**HiDEx Copyright (C) 2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014  Cyrus Shaoul and Geoff Hollis \n**This is HiDEx Version " << versionNumber << "\n**This program comes with ABSOLUTELY NO WARRANTY\n**This is free software, and you are welcome to redistribute it\n**under certain conditions; see the included file COPYING.txt for details.\n**(This is program is covered by the GPLv3).\n\nThe authors kindly request that if you report any findings that involved\nthe use of this software, please make sure to cite this\nsoftware in your report. For the correct way to cite this software,\nplease see the documentation.\nBeginning execution.\n\n" << endl; 
  
}

} // namespace



int main(int argc, char *argv[])
{
  
  // Start timing execution.
  time_t start = time(NULL);
  
  try
    {
      warranty();
      
      // Test for Open MP settings.
      char * numthreads = getenv("OMP_NUM_THREADS");
      if (numthreads) {
	cerr << "Parallel processing enabled. Number of OpenMP threads to be used is: " << numthreads << endl;
      } else {
	cerr << "OpenMP not active. Continuing without parallel processing.\nPlease set OMP_NUM_THREADS in your environment to use OpenMP.\n " << endl;
      }

      string mode;
      string configFilename;
      string multi = "";
      int ch;  
      
      // process command line arguments
      while ((ch = getopt(argc, argv, "m:f:?")) != -1)
        {
	  switch (ch)
            {  
            case 'f':
	      configFilename = optarg;
	      break;
            case 'm':
	      mode = optarg;
	      break; 
            case '?':
	      usage(cout, argv[0]);
	      return 0;
            default:
	      usage(cerr, argv[0]);
	      throw Exception("  No options specified Exiting.");
            }
        }
      if (mode == "") {
	usage(cout, argv[0]);
	throw Exception("  No mode specified. Exiting.");
      }
      if (configFilename == "")
	{
	  configFilename = "config.txt";
	  cerr << "Using default config file: \"config.txt\"\n";
        }
      if (mode == "create"){
	readSettings(configFilename, multi);
	create();
      } else if (mode == "getneighbors")
        {
	  readSettings(configFilename, multi);
	  printSDs();
        } else if (mode == "getvectors")
	  {
	    readSettings(configFilename, multi);
	    printVectors();
	  } else if (mode == "getsimilarity")
	    {
	      readSettings(configFilename, multi);
	      printPairs();
	    } else if (mode == "batchneighbors")
	      {
		printmulti(configFilename);
	      } else if (mode == "batchsimilarity")
		{
		  PrintMultiPairs(configFilename);
		} else if (mode == "update")
		  {
		    readSettings(configFilename, multi);
		    update();
		  } else if (mode == "remove")
		    {
		      readSettings(configFilename, multi);
		      remove();
		    } else {
		      usage(cout, argv[0]);
		      throw Exception("  Unknown mode: " + mode);
		    }
    }
  catch (const bad_alloc& x)
    {
      cerr << "  Out of memory error: " << x.what() << endl;
      abort();
    }
  catch (const exception& e)
    {
      cerr << endl <<"ERROR: " << e.what() << endl;
      return 1;
	}
  cerr << "Execution time: " <<  (time(NULL) - start) / 60.0 << " minutes. (Walltime)"<<  endl;
  cerr << "Goodbye." << endl;
  return 0;
}
