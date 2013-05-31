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

#ifndef Matrix_H
#define Matrix_H

#include <cassert>
#include <iostream>
#include <fstream>

using namespace std;


//*****************************************************************************
//
// Structures assocciated with storing data in the matrix
//
//*****************************************************************************

//
// A structure to hold data about rows, implemented in a
// binary tree fashion. Such an implementation saves quite a
// bit of space if we know many rows will be empty
//
template<class T>
class MatrixData
{
public:
    /** @todo Make private */
    int row;

    /** @todo Make private */
    MatrixData<T> *higher;

    /** @todo Make private */
    MatrixData<T> *lower;

    /** @todo Make private */
    T *data;

    MatrixData(int row, int cols);

    MatrixData(int row, T *rowData);

    ~MatrixData();

    size_t size(int cols) const;
};


template<class T>
class Matrix
{
public:
    /** @todo Make private */
    int rows;

    /** @todo Make private */
    int cols;

    /** @todo Make private */
    MatrixData <T> *data;

    Matrix(int rows, int cols);

    ~Matrix();

    /** Add value to row,col */
    void add(int row, int col, const T& val);

    /** Add all of the vals in from to to */
    void add(const Matrix<T>* from);

    /** Set to [row][col] val */
    void set(int row, int col, const T& val);

    /** Get the val at [row][col] */
    T get(int row, int col) const;

    /** Apply weights.
        Go through every Row in the matrix, and apply the weights
        contained in the *weights array to each row. In the matrix.
        It is assumed that the length of weights is equal to the
        length of a column in this matrix
    */
    void applyWeights(const T *weights);

    /** Collapses each row in the matrix into a single cell.
        After the transformation, the matrix will only have a column
        width of 1
    */
    void collapseColumns();

    /** Erase the contents of a row */
    void eraseRow(int row);

    /** Returns the sum of all the cells in the matrix */
    T sum() const;

    /** Returns the average of all the cells in the matrix.
        @todo Better return a double (and use internally) to avoid accumulated
        rounding errors if T is an integral type
    */
    T average() const;

    /** Returns the variance of all the cells within the matrix
        @todo Better return a double (and use internally) to avoid accumulated
        rounding errors if T is an integral type
    */
    T variance() const;

    /** Get the distance between a cell in the matrix and val. */
    T cellDistance(const T& val, int row, int col) const;

    /** Get the euclidean distance between a cell in the matrix and val. */
    T cellEuclideanDistance(const T& val, int row, int col) const;


    /** Get the euclidean distance to other matrix */
    T euclideanDistance(const Matrix<T>* m) const;

    /** Turn the matrix from a tree form into an array form */
    T** toArray() const;

    /** Write all rows with data in them to a file, along with info 
        about what rows the data belong to
    */
    void toFile(const char *filename) const;

    /** Open up the specified file and load the matrix contained inside */
    static Matrix<T>* fromFile(const char *filename, int rows, int cols);

    /** Print a matrix data to the stream in human-readable form */
    void show(ostream& out) const;

    size_t size() const;

private:
    /** Check if valid row index with assertion in debug mode. */
    void checkRowIndex(int rowIndex) const;

    /** Check if valid col index with assertion in debug mode. */
    void checkColIndex(int colIndex) const;
};


//*****************************************************************************
//
// Structures and functions assocciated with storing matrix data.
// THESE FUNCTIONS SHOULD NEVER BE USED DIRECTLY!! ONLY INTERACT WITH MATRICES
// THROUGH THE matrix*() FUNCTIONS!!
//
//*****************************************************************************


//constructor with empty rows
template<class T>
MatrixData<T>::MatrixData(int row, int cols)
{
    this->row = row;
    higher = 0;
    lower  = 0;
    data = new T[cols];
    for(int i = 0; i < cols; i++)
        data[i] = 0;
}

//constructor adding row data
template<class T>
MatrixData<T>::MatrixData(int row, T *rowData)
{
    assert(rowData != 0);
    this->row = row;
    data = rowData;
    higher = 0;
    lower = 0;
}


//destructor
template<class T>
MatrixData<T>::~MatrixData()
{
  if(higher)
      delete higher;
  if(this->lower)
      delete lower;
  delete[] data;
}

template<class T>
size_t MatrixData<T>::size(int cols) const
{
    size_t result = sizeof(MatrixData<T>);
    if (higher != 0)
        result += higher->size(cols);
    if (lower != 0)
        result += lower->size(cols);
    if (data != 0)
        result += cols * sizeof(T); // check if data always Cols elements
    //    assert(false);
    return result;

}

// @todo Make member function of MatrixData
template<class T>
void matrixDataEraseRow(MatrixData<T> *D, const int row, const int maxCols)
{
  if(D->row == row) {
    for(int i = 0; i < maxCols; i++)
      D->data[i] = 0;
  }
  else if(D->row > row && D->lower)
    matrixDataEraseRow(D->lower, row, maxCols);
  else if(D->row < row && D->higher)
    matrixDataEraseRow(D->higher, row, maxCols);
}


/*
template<class T>
void matrixDataEraseRow(MatrixData<T> *D, const int row) {
  if(D->higher && D->higher->row == row) {
    MatrixData<T> *lower = D->higher->lower;
    MatrixData<T> *higher = D->higher->higher;
    D->higher->lower = NULL;
    D->higher->higher = NULL;
    deleteMatrixData(D->higher);
    D->higher = NULL;
    
    if(!(lower || higher))
      return;
    else if(lower && higher) {
      D->higher = higher;
      matrixDataInsertRow(higher, lower);
    }
    else if(lower && !higher)
      D->higher = lower;
    else// if(higher && !lower)
      D->higher = higher;
  }

  else if(D->lower && D->lower->row == row) {
    MatrixData<T> *lower = D->lower->lower;
    MatrixData<T> *higher = D->lower->higher;
    D->lower->lower = NULL;
    D->lower->higher = NULL;
    deleteMatrixData(D->lower);
    D->lower = NULL;
    
    if(!(lower || higher))
      return;
    else if(lower && higher) {
      D->lower = higher;
      matrixDataInsertRow(higher, lower);
    }
    else if(lower && !higher)
      D->lower = lower;
    else// if(higher && !lower)
      D->lower = higher;
  }

  else {
    if(row > D->row && D->higher)
      matrixDataEraseRow(D->higher, row);
    else if(row < D->row && D->lower)
      matrixDataEraseRow(D->lower, row);
  }
}
*/


// @todo Make member function of MatrixData
template<class T>
void matrixDataAdd(MatrixData<T> *D, const int row, const T *rowData,
                   const int maxCol)
{
  if(D->row == row) {
    for(int i = 0; i < maxCol; i++)
      D->data[i] += rowData[i];
  }
  else if(D->row > row) {
    if(D->lower == 0) {
      D->lower = new MatrixData<T>(row, maxCol);
    }
    matrixDataAdd(D->lower, row, rowData, maxCol);
  }
  else {
    if(D->higher == 0) {
      D->higher = new MatrixData<T>(row, maxCol);
    }

    matrixDataAdd(D->higher, row, rowData, maxCol);
  }
}


// @todo Make member function of MatrixData
template<class T>
void matrixDataInsertRow(MatrixData<T> *D, MatrixData<T> *row)
{
  if(D->row > row->row) {
    if(D->higher == 0)
      D->higher = row;
    else
      matrixDataInsertRow(D->higher, row);
  }
  else {
    if(D->lower == NULL)
      D->lower = row;
    else
      matrixDataInsertRow(D->lower, row);
  }
}


// @todo Make member function of MatrixData
template<class T>
void matrixDataAdd(MatrixData<T> *to, const MatrixData<T> *from,
                   const int maxCol)
{
  if(from->higher) matrixDataAdd(to, from->higher, maxCol);
  if(from->lower)  matrixDataAdd(to, from->lower, maxCol);
  matrixDataAdd(to, from->row, from->data, maxCol);
}


// @todo Make member function of MatrixData
template<class T>
void matrixDataSet(MatrixData<T> *D, const int row, T *rowData)
{
  if(D->row == row) {
    if(D->data) delete[] D->data;
    D->data = rowData;
  }
  else if(D->row > row) {
    if(D->lower == 0) {
        D->lower = new MatrixData<T>(row, rowData) ;
    }
    else
      matrixDataSet(D->lower, row, rowData);
  }
  else {
    if(D->higher == 0) {
        D->higher = new MatrixData<T>(row, rowData);
    }
    else
      matrixDataSet(D->higher, row, rowData);
  }
}


// @todo Make member function of MatrixData
template<class T>
void matrixDataAdd(MatrixData<T> *D, const int row, const int col,
                   const int maxCol, const T val)
{
    int Drow = D->row;
    if(Drow == row) {
        D->data[col] += val;
    }
    else if(Drow > row) {
        if(D->lower == 0) {
            D->lower = new MatrixData<T>(row, maxCol);
        }
        matrixDataAdd(D->lower, row, col, maxCol, val);
    }
    else {
        if(D->higher == 0) {
            D->higher = new MatrixData<T>(row, maxCol);
        }
        matrixDataAdd(D->higher, row, col, maxCol, val);
    }
}


// @todo Make member function of MatrixData
template<class T>
void matrixDataSet(MatrixData<T> *D, const int row, const int col,
		   const int maxCol, const T val)
{
  if(D->row == row) {
    D->data[col] = val;
  }
  else if(D->row > row) {
    if(D->lower == 0) {
      D->lower = new MatrixData<T>(row, maxCol);
    }
    matrixDataSet(D->lower, row, col, maxCol, val);
  }
  else {
    if(D->higher == 0) {
      D->higher = new MatrixData<T>(row, maxCol);
    }
    matrixDataSet(D->higher, row, col, maxCol, val);
  }
}


// @todo Make member function of MatrixData
template<class T>
T matrixDataGet(const MatrixData<T> *D, const int row, const int col)
{
  if(D->row == row)
    return D->data[col];
  else if(D->row > row) {
    if(D->lower == 0) return 0;
    return matrixDataGet(D->lower, row, col);
  }
  else {
    if(D->higher == 0)  return 0;
    return matrixDataGet(D->higher, row, col);
  }
}


// @todo Make member function of MatrixData
template<class T>
void matrixDataApplyWeights(MatrixData<T> *D, const T *weights,
                            const int maxCol)
{
  if(D->lower) matrixDataApplyWeights(D->lower, weights, maxCol);
  if(D->higher) matrixDataApplyWeights(D->higher, weights, maxCol);

  for(int i = 0; i < maxCol; i++)
    D->data[i] *= weights[i];
}


// @todo Make member function of MatrixData
template<class T>
void matrixDataCollapseColumns(MatrixData<T> *D, const int maxCol)
{
  if(D->lower) matrixDataCollapseColumns(D->lower, maxCol);
  if(D->higher) matrixDataCollapseColumns(D->higher, maxCol);

  T newval = 0;
  for(int i = 0; i < maxCol; i++)
    newval += D->data[i];

  delete [] D->data;
  D->data = new T[1];
  D->data[0] = newval;
}


// @todo Make member function of MatrixData
template<class T>
T matrixDataSum(const MatrixData<T> *D, const int maxCol)
{
  T sum = 0;
  if(D->lower) sum += matrixDataSum(D->lower, maxCol);
  if(D->higher) sum += matrixDataSum(D->higher, maxCol);

  for(int i = 0; i < maxCol; i++)
    sum += D->data[i];

  return sum;
}


// @todo Make member function of MatrixData
template<class T>
T matrixDataSummedEuclideanDistance(const MatrixData<T> *D, const T val,
				    const int maxCol)
{
  T sum = 0;
  if(D->lower) sum += matrixDataSummedEuclideanDistance(D->lower, val, maxCol);
  if(D->higher) sum += matrixDataSummedEuclideanDistance(D->higher, val, maxCol);

  for(int i = 0; i < maxCol; i++)
    sum += (D->data[i] - val) * (D->data[i] - val);

  return sum;
}


// @todo Make member function of MatrixData
template<class T>
int matrixDataRows(const MatrixData<T> *D)
{
  int count = 1;
  if(D->lower) count += matrixDataRows(D->lower);
  if(D->higher) count += matrixDataRows(D->higher);

  return count;
}


// @todo Make member function of MatrixData
template<class T>
void matrixDataFillArray(MatrixData<T> *D, T ***array)
{
  if(D->higher) matrixDataFillArray(D->higher, array);
  if(D->lower) matrixDataFillArray(D->lower, array);

  (*array)[D->row] = D->data;
}


// @todo Make member function of MatrixData
// add assertion to check matrix bounds
template<class T, class N>
MatrixData<N> *matrixDataCopy(const MatrixData<T> *D, 
                              const int minRow, const int maxRow,
                              const int minCol, const int maxCol)
{
  if(D->row > maxRow || D->row < minRow)
    return 0;

  MatrixData<N> *newData = new MatrixData<N>(D->row, 1+maxCol-minCol);
  if(D->higher)
      newData->higher = matrixDataCopy<T, N>(D->higher, minRow, maxRow, minCol, maxCol);
  if(D->lower)
      newData->lower = matrixDataCopy<T, N>(D->lower, minRow, maxRow, minCol, maxCol);

  for(int i = minCol; i <= maxCol; i++)
    newData->data[i - minCol] = D->data[i];

  return newData;
}


template<class T>
void matrixDataToStreamBad(ostream& out, MatrixData<T> *D, const int maxCols)
{
  out << "Row: " << *(&(D->row)) << " Data: " << *(D->data) << endl;
  if(D->higher) matrixDataToStream(out, D->higher, maxCols);
  if(D->lower)  matrixDataToStream(out, D->lower,  maxCols);
}


// @todo Make member function of MatrixData
template<class T>
void matrixDataToStream(ostream& out, MatrixData<T> *D, const int maxCols)
{
  out.write(reinterpret_cast<char *>( &(D->row) ),
	    sizeof(int));
  out.write(reinterpret_cast<char *>( D->data ),
	    sizeof(T) * maxCols);
  if(D->higher) matrixDataToStream(out, D->higher, maxCols);
  if(D->lower)  matrixDataToStream(out, D->lower,  maxCols);
}


//
// read from a file that contains rows of data that are of length maxCols, 
// and are preceded by a row number. All cells in the row are of type T 
// (unsigned short, int, or perhaps even a type that is more complex like 
// another row of data, or another database!). the rows do not neccessarily 
// have to be stored in sequential order.
// @todo Make member function of MatrixData
template<class T>
MatrixData<T> *matrixDataFromStream(ifstream& in, const int maxCols)
{
  // we have never created a file w/ this name
  if(!in.good())
    return 0;

  MatrixData<T> *D = 0; // *D == our datastructure for storing the data

  while(in.good()) {
    int row;
    // reinterpret_cast<char *> == we are reading binary data
    in.read(reinterpret_cast<char *>( &row ), 
	    sizeof(row));

    // if we are at the end of file, will still be good until we
    // read one more row (which turns out to be absolutely nothing.
    // and we hit EOF). This is just a check to see if the last
    // entry we tried to read made us hit EOF. If we don't have
    // this check in place, we will screw up the last row we read...
    // its sad we have to check if in.good() is actually true twice per
    // loop, but at the moment I can't think of a better way to do this
    if(!in.good())
      break;

    // create an empty row of data
    T *rowData = new T[maxCols];
    // read the next maxCols of cells into this new row
    in.read(reinterpret_cast<char *>( rowData ),
	    sizeof(T) * maxCols );

    // set the appropriate row of the matrix to 
    // the values of the data we just read in
    if(D == 0)
        D = new MatrixData<T>(row, rowData);
    else
        matrixDataSet(D, row, rowData);
  }

  return D; 
}


// @todo Make member function of MatrixData
template<class T>
void showMatrixData(ostream &out, const MatrixData<T> *D, const int maxCol)
{
  if(D->lower) showMatrixData(out, D->lower, maxCol);
 
  out << "Row #" << D->row << ": ";
  for(int i = 0; i < maxCol; i++)
    out << D->data[i] << " ";
  out << endl;

  if(D->higher) showMatrixData(out, D->higher, maxCol);
}


//*****************************************************************************
//
// Functions used for setting values in a matrix, and retreiving info about
// the matrix. THESE SHOULD BE THE ONLY FUNCTIONS USED FOR INTERACTING WITH
// MATRICES!!!
//
//*****************************************************************************

template<class T>
Matrix<T>::Matrix(int rows, int cols)
{
    assert(rows >= 0);
    assert(cols >= 0);
    this->rows = rows;
    this->cols = cols;
    data = 0;
}


template<class T>
Matrix<T>::~Matrix()
{
    if(data != 0)
        delete data;
}

template<class T>
inline void Matrix<T>::checkRowIndex(int rowIndex) const
{
  //  bool temp = (rowIndex >= 0);
  assert(rowIndex >= 0);
  //  temp = (rowIndex < rows);
  assert(rowIndex < rows);
}

template<class T>
inline void Matrix<T>::checkColIndex(int colIndex) const
{
  //  bool temp = (colIndex >= 0); 
  assert(colIndex >= 0); 
  //  temp =(colIndex < cols); 
  assert(colIndex < cols); 
}

template<class T>
void Matrix<T>::add(int row, int col, const T& val)
{
    checkRowIndex(row);
    checkColIndex(col);
    if(data == 0)
        data = new MatrixData<T>(row, cols);
    matrixDataAdd(data, row, col, cols, val);
}

template<class T>
void Matrix<T>::add(const Matrix<T>* from)
{
    if(from->data == 0) {
        return;
    }
    else if(data == 0) {
        data = matrixDataCopy<T, T>(from->data, 0, rows, 0, cols);
    }
    else
        matrixDataAdd(data, from->data, cols);
}

template<class T>
void Matrix<T>::set(int row, int col, const T& val)
{
    checkRowIndex(row);
    checkColIndex(col);
    if(data == 0)
        data = new MatrixData<T>(row, cols);
    matrixDataSet(data, row, col, cols, val);
}

template<class T>
T Matrix<T>::get(int row, int col) const
{
    checkRowIndex(row);
    checkColIndex(col);
    if(data == 0)
        return 0;
    return matrixDataGet(data, row, col);
}

template<class T>
void Matrix<T>::applyWeights(const T *weights)
{
    if(data)
        matrixDataApplyWeights(data, weights, cols);
}

template<class T>
void Matrix<T>::collapseColumns()
{
  if(data)
      matrixDataCollapseColumns(data, cols);
  cols = 1;
}

template<class T>
void Matrix<T>::eraseRow(int row)
{
    checkRowIndex(row);
    if(data)
        matrixDataEraseRow(data, row, cols);
}

template<class T>
T Matrix<T>::sum() const
{
    return matrixDataSum(data, cols);
}

template<class T>
T Matrix<T>::average() const
{
    T sum = sum();
    int cells = cols * rows;
    return sum / cells;
}

template<class T>
T Matrix<T>::variance() const
{
    T avg = average();
    T variance = matrixDataSummedEuclideanDistance(data, avg, cols);
    // account for all of the zeros not stored
    int rows_counted = matrixDataRows(data);
    variance += avg * avg * (rows - rows_counted) * cols;
    return variance / (rows * cols - 1);
}

template<class T>
T Matrix<T>::cellDistance(const T& val, int row, int col) const
{
    checkRowIndex(row);
    checkColIndex(col);
    T cell_val = 0;
    if(data)
        cell_val = matrixDataGet(data, row, col);
    return (cell_val - val);
}

template<class T>
T Matrix<T>::cellEuclideanDistance(const T& val, int row, int col) const
{
    T distance = cellDistance(val, row, col);
    return distance * distance;
}

template<class T>
T Matrix<T>::euclideanDistance(const Matrix<T>* m) const
{
    T distance = 0;
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            distance += cellEuclideanDistance(m->get(i, j), i, j);
    return distance;
}

template<class T>
size_t Matrix<T>::size() const
{
    if (data != 0)
        return sizeof(Matrix<T>) + data->size(cols);
    else
        return sizeof(Matrix<T>);
}

template<class T>
T** Matrix<T>::toArray() const
{
    T **array = new T*[rows];
    for(int i = 0; i < rows; i++)
        array[i] = 0;
    if(data)
        matrixDataFillArray(data, &array);
    return array;
}

template<class T>
void Matrix<T>::toFile(const char* filename) const
{
    //    cerr << "Entering toFile. " << endl;
    if(data == 0)
    {
        //        cerr << "toFile: No data in Matrix. Not doing anything.";
        return;
    }
    ofstream out;
    out.open(filename, ios::out | ios::binary);
    matrixDataToStream(out, data, cols);
    out.close();
}

template<class T>
Matrix<T>* Matrix<T>::fromFile(const char *filename, int rows, int cols)
{
    assert(rows >= 0);
    assert(cols >= 0);
    ifstream in;
    in.open(filename, ios::in | ios::binary);
    MatrixData<T> *data = matrixDataFromStream<T>(in, cols);
    if(data == 0)
        return 0;
    Matrix<T>* m = new Matrix<T>(rows, cols);
    m->data = data;
    return m;
}

template<class T>
void Matrix<T>::show(ostream& out) const
{
    if(data)
        showMatrixData(out, data, cols);
    else
        out << "matrix contains no data to show\n";
}

/** Return a copy of the matrix, with the given dimensions.
    Can also be used to change the type of the matrix (e.g. from ints to
    doubles)
    @todo This is the last function that accesses members of Matrix directly
    and therefore keeps them from being made private. Since this function
    depends on a second template parameter, it cannot be a member function
    of Matrix. One option would be to replace it by a copy constructor and
    assignment operator in Matrix for the case that T==N (the conversion
    functionality for T!=N is not used elsewhere in the code anyway)
    
*/

// add assertion to check matrix bounds
template<class T, class N>
Matrix<N> *matrixCopy(const Matrix<T> *M, 
                      const int rowStart, const int rowEnd,
                      const int colStart, const int colEnd)
{
  Matrix<N> *Mnew = new Matrix<N>(1+rowEnd-rowStart, 1+colEnd-colStart);
  if(M->data)
      Mnew->data = matrixDataCopy<T, N>(M->data, rowStart, rowEnd, 
                                        colStart, colEnd);
  return Mnew;
}

#endif
