#include "Allele.hpp"

  Allele::Allele(int const key):allele(key)
  {}

  Allele::~Allele()
  {}

  Allele::getAllele()
  {
    return allele;
  }

  setAllele(int const id)
  {
    allele=id;
  }
