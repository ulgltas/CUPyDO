/*!
 * Header for CManager.
 * Authors : D. THOMAS.
 *
 * COPYRIGHT (C) University of Liège, 2017.
 */

#pragma once

#include <string>
#include <vector>

class CManager{

public:
  CManager();
  virtual ~CManager();
  int getGlobalIndex(std::string const& str_physics, int const& iProc, int const& iVertex);
  void setGlobalIndexing(std::string str_physics, std::vector<std::vector<int> > index_range);
  std::vector< std::vector< std::vector<int> > > globalIndexRange;
  int nPhyscis;
  int nIndex;
  int mpiSize;
};
