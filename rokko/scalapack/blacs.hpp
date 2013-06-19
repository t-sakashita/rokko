#ifndef ROKKO_BLACS_HPP
#define ROKKO_BLACS_HPP

#ifdef __cplusplus
extern "C" {
#endif
  
void descinit_(int* desc, const int* m, const int* n, const int* mb, const int* nb,
               const int* irsrc, const int* icsrc, const int* ixtxt, const int* lld, int* info);

inline void ROKKO_descinit(int* desc, int m, int n, int mb, int nb, int irsrc, int icsrc,
                           int ixtxt, int lld, int* info) {
  descinit_(desc, &m, &n, &mb, &nb, &irsrc, &icsrc, &ixtxt, &lld, info);
}

void blacs_pinfo_(int* mypnum, int* nprocs);

inline void ROKKO_blacs_pinfo(int* mypnum, int* nprocs) {
  blacs_pinfo_(mypnum, nprocs);
}
  
void blacs_get_(const int* context, const int* request, int* value);

inline void ROKKO_blacs_get(int context, int request, int* value) {
  blacs_get_(&context, &request, value);
}
  
void blacs_gridinfo_(const int* ictxt, int* nprow, int* npcol, int* myrow, int* mycol);

inline void ROKKO_blacs_gridinfo(int ictxt, int* nprow, int* npcol, int* myrow, int* mycol) {
  blacs_gridinfo_(&ictxt, nprow, npcol, myrow, mycol);
}
  
void blacs_gridinit_(const int* ictxt, const char* order, const int* nprow, const int* npcol );

inline void ROKKO_blacs_gridinit(int ictxt, char order, int nprow, int npcol) {
  blacs_gridinit_(&ictxt, &order, &nprow, &npcol);
}

void blacs_gridexit_(int* ictxt);

inline void ROKKO_blacs_gridexit(int* ictxt) { blacs_gridexit_(ictxt); }

void blacs_exit_(const int* conti);

inline void ROKKO_blacs_exit(int conti) { blacs_exit_(&conti); }

void blacs_barrier_(const int* ictxt, const char* score);

inline void ROKKO_blacs_barrier(int ictxt, char score) {
  blacs_barrier_(&ictxt, &score);
}

#ifdef __cplusplus
} // extern "C"
#endif

#endif // ROKKO_BLACS_HPP
