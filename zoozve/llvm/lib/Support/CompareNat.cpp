//
// Created by jianglimin on 23-12-20.
//

#include "llvm/Support/CompareNat.h"
#include <algorithm>

using namespace llvm;

bool CompareNat::cmpNat(const std::string& a, const std::string& b) {
  if (a.empty()) return true;
  if (b.empty()) return false;
  if ( std::isdigit(a[0]) && !std::isdigit(b[0])) return true;
  if (!std::isdigit(a[0]) &&  std::isdigit(b[0])) return false;
  if (!std::isdigit(a[0]) && !std::isdigit(b[0])) {
    if (std::toupper(a[0]) == std::toupper(b[0]))
      return cmpNat(a.substr(1), b.substr(1));
    return (std::toupper(a[0]) < std::toupper(b[0]));
  }

  // Both strings begin with a digit --> parse both numbers
  std::istringstream issa(a);
  std::istringstream issb(b);
  int ia, ib;
  issa >> ia;
  issb >> ib;
  if (ia != ib) return ia < ib;

  // Numbers are the same --> remove numbers and recurse
  std::string anew, bnew;
  std::getline(issa, anew);
  std::getline(issb, bnew);
  return (cmpNat(anew, bnew));
}

std::vector<size_t> CompareNat::Alpha2Natural(const int N) {
  std::vector<std::string> str;
  std::vector<size_t> idx_list;

  for (int i = 1; i <= N; i++) {
    str.push_back(std::to_string(i));
  }

  // In alphabetical order
  std::sort(str.begin(), str.end());

  // In natural order
  for (auto i : this->sort_indices(str)) idx_list.push_back(i);

  return idx_list;
}
