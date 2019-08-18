//
//  SparseMatrix.h
//  
//
//  Created by Ling Zou on 9/19/13.
//
//

#ifndef _SparseMatrix_h
#define _SparseMatrix_h

#include <vector>

/*
 * This class serves as an interface between our code and the solver.
 * We do not try to implement any additonal functions of sparse matrix.
 *
 * How it works:
 * ColumnEntry stores an entry: column number and its value
 * The row number is determined by its position in the vector.
 */
class SparseMatrix
{
public:
  SparseMatrix()
  {}
  
  ~SparseMatrix()
  {}

  typedef struct
  {
    int     col_number;
    double  value;
  } ColumnEntry;
  
  /*
   * Create empty rows for future ColumnEntry to be inserted
   */
  void SetRowNumber(unsigned int total_row_number);
  
  /*
   * Add an entry at given position (row, col)
   */
  void AddEntry(unsigned int row, unsigned int col, double value = 1);
  
  /*
   * Insert a SparseMatrix data structure to this SparseMatrix
   * To the position, i-th row and j-th column (index from 0)
   */
  void Insert(const SparseMatrix & sparse_matrix, unsigned int row, unsigned int col);
  
  const std::vector<std::vector<ColumnEntry> > & getEntryData() const { return EntryData; }

protected:
  std::vector<std::vector<ColumnEntry> > EntryData;
};

#endif
