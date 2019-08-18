//
//  SparseMatrix.C
//  
//
//  Created by Ling Zou on 9/19/13.
//
//

#include <stdio.h>
#include <iostream>
#include <sstream>
#include <cstdlib> // exit

#include "SparseMatrix.h"

void
SparseMatrix::SetRowNumber(unsigned int total_row_number)
{
  // Initialize the vector
  EntryData.resize(total_row_number);
}

void
SparseMatrix::AddEntry(unsigned int row, unsigned int col, double value)
{
  ColumnEntry entry;
  entry.col_number = col;
  entry.value = value;
  
  if(row < EntryData.size())
    EntryData[row].push_back(entry);
  else
  {
    std::cerr << "ERROR: The ColumnEntry being added is out of range. It's index, " << row << " (starts from 0), is beyond the EntryData size, " << EntryData.size() << "\n";
    exit(1);
  }
}

void
SparseMatrix::Insert(const SparseMatrix & sparse_matrix, unsigned int row, unsigned int col)
{
  const std::vector<std::vector<ColumnEntry> > & data_entries_being_inserted = sparse_matrix.getEntryData();
  
  if((EntryData.size() < data_entries_being_inserted.size() + row) || (EntryData.size() < data_entries_being_inserted.size() + col))
  {
    std::cerr << "Not enough room for the SparseMatrix being inserted.\n";
    exit(1);
  }
  
  for(int i = 0; i < data_entries_being_inserted.size(); i++)
    for(int j = 0; j < data_entries_being_inserted[i].size(); j++)
      AddEntry(i + row, data_entries_being_inserted[i][j].col_number + col, data_entries_being_inserted[i][j].value);
}