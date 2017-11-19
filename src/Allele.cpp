#include "Allele.hpp"

  Allele::Allele(int const key):allele(key)
  {}

  Allele::~Allele()
  {}

  int Allele::getAllele() const
  {
    return allele;
  }

 void Allele::setAllele(int const id)
  {
    allele=id;
  }
