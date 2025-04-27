//
// Created by jianglimin on 23-12-19.
// https://stackoverflow.com/questions/9743485/natural-sort-of-directory-filenames-in-c
// https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
//

#ifndef LLVM_COMPARENAT_H
#define LLVM_COMPARENAT_H

#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

namespace llvm {

class CompareNat {
public:
  CompareNat() { }

  /// An algorithm that arranges strings into alphabetical order
  bool cmpNat(const std::string& a, const std::string& b);

  /// Return the indices change from alphabetical to natural order
  std::vector<size_t> Alpha2Natural(const int N);

  template <typename T>
  std::vector<size_t> sort_indices(const std::vector<T> &v) {

    // initialize original index locations
    std::vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values
    stable_sort(idx.begin(), idx.end(),
                [&v](size_t i1, size_t i2)
                { return stoi(v[i1]) < stoi(v[i2]); });

    return idx;
  }

};

} // end namespace llvm

#endif // LLVM_COMPARENAT_H
