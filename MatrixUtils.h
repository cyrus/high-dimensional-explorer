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

#ifndef MatrixUtils_H
#define MatrixUtils_H

#include "SDDB.h"
#include "Exception.h"

// Why are these two functions here? Because they use templates, and must be in a .h file!


/** Collapse matrix.
    @param M Sparse matrix of co-occurrence values that has dimensions  x and
    y, where x is the window size and y is the size of the lexicon.
    @param weights Array of integers that are used to multiply the 
    values in each position of the window before summing them.
    @param realBehind The actual size of the backward window that is stored in
    the db.
    @param behind The desired size of the backward window
    @param realAhead see above
    @param ahead see above
    @param context Array of integers that contains all the word-IDs of the 
    words that are part of the context (in this case, the N words with the 
    highest raw frequencies in the corpus)
    @param context_size is the size of the context array
    @param separate Flag for  "separate mode". If it is TRUE, then 
    forward and backward values are summed separately. If it is FALSE, then 
    the forward and backward values are summed together. If true, it doubles 
    the size of the final matrix.
    @param num_dimensions if separate is FALSE, it is the context size. If 
    separate is TRUE, it is 2 * context_size.
*/



/* void normalizeWeightedCooccurenceVector(Float cooccurenceVector[], int size, Float wordfrequency) */
/* { */
  
/*   // normalize vectors by dividing by ofreq. */
/*   // make reciprocal of 1 over non-zero frequncy */
/*   // cerr << "Normalizing Word Frequency = " << wordfrequency << endl; */
/*   Float wordfrequency = static_cast<Float>(_frequency[i]); */
/*   Float adjustedfreq = (1.0 / (wordfrequency + 1.0)); */
/*   for(int i = 0; i < size; i++) { */
/*     if (cooccurenceVector[i]) { */
/*       //      cerr << "Adjusted real freq " << cooccurenceVector[i] << " to "; */
/*       cooccurenceVector[i] *= adjustedfreq; */
/*       //      cerr << cooccurenceVector[i] << endl; */
/*     }  */
/*   } */
/* } */


template <class T>
Float*computeWeightedCooccurenceVectorUnified(Matrix<T> *M, 
					      const vector<int>& weights, int realBehind,
					      int behind, int realAhead, 
					      int ahead, vector<int>& frequentWords)
{
  
  //    cerr << "Computing aggregate" << endl;
  
  // Sanity check for window ranges.
  //  bool temp = (ahead <= realAhead);
  assert(ahead <= realAhead);
  //  temp = (behind <= realBehind);
  assert(behind <= realBehind);
  
  int windowLen = behind+ahead;
  int startpoint = realBehind-behind;
  
  size_t contextsize = frequentWords.size();

  //    cerr << "converting array"<< endl;
  int **array = M->toArray();
  
  //    cerr << "Make new array of Floats."<< endl;
  Float *vect = new Float[contextsize];
  
  //    cerr << "Start aggregating vectors."<< endl;
  for(size_t i = 0; i < contextsize; i++) {
    // branch if using separate forward and backward contexts
    Float val = 0.0;
    // collapse all of the columns into a single cell
    //        cerr << "Check for NULLs."<< endl;
    if (array[frequentWords[i]] == NULL) {
      val = 0.0;
    }
    else {
      //            cerr << "Multiply values"<< endl;
      for (int k = 0; k < windowLen; k++) {
	//                cerr << " Value : " << i << "," << k << " = " << array[frequentWords[i]][k+startpoint]<< endl;
	val += static_cast<Float>((array[frequentWords[i]][k+startpoint] * weights[k]));
      }
    } 
    vect[i] = val;
    //        cerr << "Set aggregated value: " << val << endl;        
  }
  //  normalizeWeightedCooccurenceVector(vect, frequentWords.size(), wordfrequency, normalization);
  //    cerr << "Delete array"<< endl;    
  delete [] array;
  return vect;
}

template <class T>
Float*computeWeightedCooccurenceVectorSeparate(Matrix<T> *M, const vector<int>& weights, int realBehind,
                      int behind, int realAhead, int ahead, vector<int>& frequentWords)

{
  //  bool temp = (ahead <= realAhead);
  assert(ahead <= realAhead);
  //  temp = (behind <= realBehind);
  assert(behind <= realBehind);
  
  size_t num_dimensions = 2 * frequentWords.size();
  size_t context_size = frequentWords.size();
  int windowLen = behind+ahead;
  int startpoint = realBehind-behind;
  
  int **array = M->toArray();
  
  Float *vect = new Float[num_dimensions];
  
  for (size_t p=0; p < num_dimensions; p++)
    vect[p]=0.0;
  
  // branch if using separate forward and backward contexts
  //    if (windowLen % 2 != 0)
  //    {
  //        throw Exception("Window length was not an even number!  exiting!");
  //    }
  
  // debugging... printing data....
  /*     int fullw = realAhead + realBehind; */
  /*     for (size_t x = 0; x < 100; x++) { */
  /*         if (array[x] != NULL) { */
  /*             cerr << "Row #" << x << " :"; */
  /*             for (int y = 0; y < fullw; y++) { */
  /*                 cerr << array[x][y] << " "; */
  /*             } */
  /*         cerr << endl; */
  /*         } */
  /*     } */
  
  for(size_t i = 0; i < context_size; i++)
    {
      Float valforw = 0.0;
      Float valback = 0.0;
      if (array[frequentWords[i]] == NULL)
        {
	  valforw = 0.0;
	  valforw = 0.0;
        }
      else
        {
	  //	  cerr << "Full Window Width = " << windowLen << endl;
	  //	  cerr << "Window Start Point  = " << startpoint << endl;
	  //	  cerr << "For Word: " << frequentWords[i] << " -- " ; 
            for (int k = 0; k < (windowLen - ahead); k++)
	      {
                valback += static_cast<Float>((array[frequentWords[i]][k+startpoint] * weights[k]));
		//		cerr << "Position = " << k << ", Valback = " << valback << " , ";
	      }
            for (int k = (windowLen - ahead); k < windowLen; k++)
	      {
                valforw += static_cast<Float>((array[frequentWords[i]][k+startpoint] * weights[k]));
		//		cerr << "Position = " << k << " Valforw = " << valforw << " "; 
	      }
	    //	    cerr << endl;
        }
      vect[i] = valback;
      assert((i+ context_size) < num_dimensions);
      vect[i + context_size] = valforw;
    }

  //debug
  //    for (int x = 0; x < num_dimensions; x++) {
  //        if (vect[x] != 0.0)
  //            cerr << "At: " << x << " value: " << vect[x] << endl; 
    //    }
  
  //  normalizeWeightedCooccurenceVector(vect, num_dimensions, wordfrequency);
  
  //    cerr << "Num Dimensions = " << num_dimensions << endl;
  //    for (int x = 0; x < num_dimensions; x++) {
  //        if (vect[x] != 0.0)
  //            cerr << "At: " << x << " value: " << vect[x] << endl; 
  //    }

  delete [] array;
  return vect;
  
}


template <class T>
Float computeVariance(Matrix<T> *M, 
		      int realBehind,
		      int realAhead,
		      size_t numdimensions
		      )		    
{
  //see http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance  
  // for info on algorithm.

  //  cerr << "Computing Variance." << endl;
  
  //    cerr << "converting array"<< endl;
  int **array = M->toArray();
  
  //    initialize variables.
  Float variance = 0.0;
  Float n = 0.0;
  Float mean = 0.0;
  Float M2 = 0.0;
  Float delta = 0.0;
  Float x = 0.0;
  
  int windowLen = realBehind+realAhead;

  //    cerr << "Start aggregating vectors."<< endl;
  for(size_t i = 0; i < numdimensions; i++) {
    // collapse all of the columns into a single cell
    for (int k = 0; k < windowLen; k++) {
      n++;
      x = static_cast<Float>(array[i][k]);
      delta = x - mean;
      mean = mean + (delta/n);
      M2 = M2 + (delta * (x - mean));
      }
  }
  variance = M2/(n-1);
  cerr << "Variance = " << variance << endl;        
  }
  //    cerr << "Delete array"<< endl;    
  delete [] array;
  return variance;
}


#endif // MatrixUtils_H
