#ifndef ROKKO_BLACS_HPP
#define ROKKO_BLACS_HPP


extern "C" {
  /* BLACS */
  void blacs_pinfo_( int& mypnum, int& nprocs);
  void blacs_get_( const int& context, const int& request, int& value );
  void blacs_gridinfo_(const int& ictxt, int& nprow, int& npcol, int& myrow, int& mycol);
  void blacs_gridinit_(const int& ictxt, char* order, int& nprow, int& npcol );
  void blacs_gridexit_(int* ictxt);
  void blacs_exit_(const int& conti);
  void blacs_barrier_(const int& ictxt, const char* score);

  int numroc_(const int& n, const int& nb, const int& iproc, const int& isrcproc, const int& nprocs);

}

#endif //ROKKO_BLACS_HPP
