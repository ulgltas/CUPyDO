/*
 * Copyright 2018 University of Li�ge
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*!
 * Header for CManager.
 * Authors : D. THOMAS.
 */

#ifndef CMANAGER_H
#define CMANAGER_H

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

#endif //CMANAGER_H
