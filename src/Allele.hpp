#ifndef ALLELE_H
#define ALLELE_H

class Allele {
  public:

  Allele(int const key);//!<Constructor
  ~Allele();//!<Destructor

  /**
  * @return int Allele
  */
  getAllele() const;

  /**
  *@param int const id
  */
  setAllele(int const id);

  private:
  int allele const;//!< identifier specific to the allele
};

#endif
