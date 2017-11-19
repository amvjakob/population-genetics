#ifndef ALLELE_H
#define ALLELE_H

class Allele {
  public:

  Allele(int const key);//!<Constructor
  ~Allele();//!<Destructor

  /**
  * @return int Allele
  */
  int getAllele() const;

  /**
  *@param int const id
  */
  void setAllele(int const id);

  private:
  int allele;//!< identifier specific to the allele
};

#endif
