/*
 * cblas_f77.h
 * Written by Keita Teranishi
 *
 * Updated by Jeff Horner
 * Merged cblas_f77.h and cblas_fortran_header.h
 */

#ifndef CBLAS_F77_H
#define CBLAS_f77_H

#ifdef CRAY
   #include <fortran.h>
   #define F77_CHAR _fcd
   #define C2F_CHAR(a) ( _cptofcd( (a), 1 ) )
   #define C2F_STR(a, i) ( _cptofcd( (a), (i) ) )
   #define F77_STRLEN(a) (_fcdlen)
#endif

#ifdef WeirdNEC
   #define F77_INT long
#endif

#ifdef  F77_CHAR
   #define FCHAR F77_CHAR
#else
   #define FCHAR char *
#endif

#ifdef F77_INT
   #define FINT const F77_INT *
   #define FINT2 F77_INT *
#else
   #define FINT const int *
   #define FINT2 int *
#endif

#include "cblas_mangling.h"

#define F77_xerbla     BLAS_GLOBAL(xerbla,    XERBLA)
/*
 * Level 1 BLAS
 */
#define F77_srotg      BLAS_GLOBAL(srotg,     SROTG)
#define F77_srotmg     BLAS_GLOBAL(srotmg,    SROTMG)
#define F77_srot       BLAS_GLOBAL(srot,      SROT)
#define F77_srotm      BLAS_GLOBAL(srotm,     SROTM)
#define F77_drotg      BLAS_GLOBAL(drotg,     DROTG)
#define F77_drotmg     BLAS_GLOBAL(drotmg,    DROTMG)
#define F77_drot       BLAS_GLOBAL(drot,      DROT)
#define F77_drotm      BLAS_GLOBAL(drotm,     DROTM)
#define F77_sswap      BLAS_GLOBAL(sswap,     SSWAP)
#define F77_scopy      BLAS_GLOBAL(scopy,     SCOPY)
#define F77_saxpy      BLAS_GLOBAL(saxpy,     SAXPY)
#define F77_isamax_sub BLAS_GLOBAL(isamaxsub, ISAMAXSUB)
#define F77_dswap      BLAS_GLOBAL(dswap,     DSWAP)
#define F77_dcopy      BLAS_GLOBAL(dcopy,     DCOPY)
#define F77_daxpy      BLAS_GLOBAL(daxpy,     DAXPY)
#define F77_idamax_sub BLAS_GLOBAL(idamaxsub, IDAMAXSUB)
#define F77_cswap      BLAS_GLOBAL(cswap,     CSWAP)
#define F77_ccopy      BLAS_GLOBAL(ccopy,     CCOPY)
#define F77_caxpy      BLAS_GLOBAL(caxpy,     CAXPY)
#define F77_icamax_sub BLAS_GLOBAL(icamaxsub, ICAMAXSUB)
#define F77_zswap      BLAS_GLOBAL(zswap,     ZSWAP)
#define F77_zcopy      BLAS_GLOBAL(zcopy,     ZCOPY)
#define F77_zaxpy      BLAS_GLOBAL(zaxpy,     ZAXPY)
#define F77_izamax_sub BLAS_GLOBAL(izamaxsub, IZAMAXSUB)
#define F77_sdot_sub   BLAS_GLOBAL(sdotsub,   SDOTSUB)
#define F77_ddot_sub   BLAS_GLOBAL(ddotsub,   DDOTSUB)
#define F77_dsdot_sub  BLAS_GLOBAL(dsdotsub,  DSDOTSUB)
#define F77_sscal      BLAS_GLOBAL(sscal,     SSCAL)
#define F77_dscal      BLAS_GLOBAL(dscal,     DSCAL)
#define F77_cscal      BLAS_GLOBAL(cscal,     CSCAL)
#define F77_zscal      BLAS_GLOBAL(zscal,     ZSCAL)
#define F77_csscal     BLAS_GLOBAL(csscal,    CSSCAL)
#define F77_zdscal     BLAS_GLOBAL(zdscal,    ZDSCAL)
#define F77_cdotu_sub  BLAS_GLOBAL(cdotusub,  CDOTUSUB)
#define F77_cdotc_sub  BLAS_GLOBAL(cdotcsub,  CDOTCSUB)
#define F77_zdotu_sub  BLAS_GLOBAL(zdotusub,  ZDOTUSUB)
#define F77_zdotc_sub  BLAS_GLOBAL(zdotcsub,  ZDOTCSUB)
#define F77_snrm2_sub  BLAS_GLOBAL(snrm2sub,  SNRM2SUB)
#define F77_sasum_sub  BLAS_GLOBAL(sasumsub,  SASUMSUB)
#define F77_dnrm2_sub  BLAS_GLOBAL(dnrm2sub,  DNRM2SUB)
#define F77_dasum_sub  BLAS_GLOBAL(dasumsub,  DASUMSUB)
#define F77_scnrm2_sub BLAS_GLOBAL(scnrm2sub, SCNRM2SUB)
#define F77_scasum_sub BLAS_GLOBAL(scasumsub, SCASUMSUB)
#define F77_dznrm2_sub BLAS_GLOBAL(dznrm2sub, DZNRM2SUB)
#define F77_dzasum_sub BLAS_GLOBAL(dzasumsub, DZASUMSUB)
#define F77_sdsdot_sub BLAS_GLOBAL(sdsdotsub, SDSDOTSUB)
/*
 * Level 2 BLAS
 */
#define F77_ssymv      BLAS_GLOBAL(ssymv, SSYMV)
#define F77_ssbmv      BLAS_GLOBAL(ssbmv, SSBMV)
#define F77_sspmv      BLAS_GLOBAL(sspmv, SSPMV)
#define F77_sger       BLAS_GLOBAL(sger,  SGER)
#define F77_ssyr       BLAS_GLOBAL(ssyr,  SSYR)
#define F77_sspr       BLAS_GLOBAL(sspr,  SSPR)
#define F77_ssyr2      BLAS_GLOBAL(ssyr2, SSYR2)
#define F77_sspr2      BLAS_GLOBAL(sspr2, SSPR2)
#define F77_dsymv      BLAS_GLOBAL(dsymv, DSYMV)
#define F77_dsbmv      BLAS_GLOBAL(dsbmv, DSBMV)
#define F77_dspmv      BLAS_GLOBAL(dspmv, DSPMV)
#define F77_dger       BLAS_GLOBAL(dger,  DGER)
#define F77_dsyr       BLAS_GLOBAL(dsyr,  DSYR)
#define F77_dspr       BLAS_GLOBAL(dspr,  DSPR)
#define F77_dsyr2      BLAS_GLOBAL(dsyr2, DSYR2)
#define F77_dspr2      BLAS_GLOBAL(dspr2, DSPR2)
#define F77_chemv      BLAS_GLOBAL(chemv, CHEMV)
#define F77_chbmv      BLAS_GLOBAL(chbmv, CHBMV)
#define F77_chpmv      BLAS_GLOBAL(chpmv, CHPMV)
#define F77_cgeru      BLAS_GLOBAL(cgeru, CGERU)
#define F77_cgerc      BLAS_GLOBAL(cgerc, CGERC)
#define F77_cher       BLAS_GLOBAL(cher,  CHER)
#define F77_chpr       BLAS_GLOBAL(chpr,  CHPR)
#define F77_cher2      BLAS_GLOBAL(cher2, CHER2)
#define F77_chpr2      BLAS_GLOBAL(chpr2, CHPR2)
#define F77_zhemv      BLAS_GLOBAL(zhemv, ZHEMV)
#define F77_zhbmv      BLAS_GLOBAL(zhbmv, ZHBMV)
#define F77_zhpmv      BLAS_GLOBAL(zhpmv, ZHPMV)
#define F77_zgeru      BLAS_GLOBAL(zgeru, ZGERU)
#define F77_zgerc      BLAS_GLOBAL(zgerc, ZGERC)
#define F77_zher       BLAS_GLOBAL(zher,  ZHER)
#define F77_zhpr       BLAS_GLOBAL(zhpr,  ZHPR)
#define F77_zher2      BLAS_GLOBAL(zher2, ZHER2)
#define F77_zhpr2      BLAS_GLOBAL(zhpr2, ZHPR2)
#define F77_sgemv      BLAS_GLOBAL(sgemv, SGEMV)
#define F77_sgbmv      BLAS_GLOBAL(sgbmv, SGBMV)
#define F77_strmv      BLAS_GLOBAL(strmv, STRMV)
#define F77_stbmv      BLAS_GLOBAL(stbmv, STBMV)
#define F77_stpmv      BLAS_GLOBAL(stpmv, STPMV)
#define F77_strsv      BLAS_GLOBAL(strsv, STRSV)
#define F77_stbsv      BLAS_GLOBAL(stbsv, STBSV)
#define F77_stpsv      BLAS_GLOBAL(stpsv, STPSV)
#define F77_dgemv      BLAS_GLOBAL(dgemv, DGEMV)
#define F77_dgbmv      BLAS_GLOBAL(dgbmv, DGBMV)
#define F77_dtrmv      BLAS_GLOBAL(dtrmv, DTRMV)
#define F77_dtbmv      BLAS_GLOBAL(dtbmv, DTBMV)
#define F77_dtpmv      BLAS_GLOBAL(dtpmv, DTPMV)
#define F77_dtrsv      BLAS_GLOBAL(dtrsv, DTRSV)
#define F77_dtbsv      BLAS_GLOBAL(dtbsv, DTBSV)
#define F77_dtpsv      BLAS_GLOBAL(dtpsv, DTPSV)
#define F77_cgemv      BLAS_GLOBAL(cgemv, CGEMV)
#define F77_cgbmv      BLAS_GLOBAL(cgbmv, CGBMV)
#define F77_ctrmv      BLAS_GLOBAL(ctrmv, CTRMV)
#define F77_ctbmv      BLAS_GLOBAL(ctbmv, CTBMV)
#define F77_ctpmv      BLAS_GLOBAL(ctpmv, CTPMV)
#define F77_ctrsv      BLAS_GLOBAL(ctrsv, CTRSV)
#define F77_ctbsv      BLAS_GLOBAL(ctbsv, CTBSV)
#define F77_ctpsv      BLAS_GLOBAL(ctpsv, CTPSV)
#define F77_zgemv      BLAS_GLOBAL(zgemv, ZGEMV)
#define F77_zgbmv      BLAS_GLOBAL(zgbmv, ZGBMV)
#define F77_ztrmv      BLAS_GLOBAL(ztrmv, ZTRMV)
#define F77_ztbmv      BLAS_GLOBAL(ztbmv, ZTBMV)
#define F77_ztpmv      BLAS_GLOBAL(ztpmv, ZTPMV)
#define F77_ztrsv      BLAS_GLOBAL(ztrsv, ZTRSV)
#define F77_ztbsv      BLAS_GLOBAL(ztbsv, ZTBSV)
#define F77_ztpsv      BLAS_GLOBAL(ztpsv, ZTPSV)
/*
 * Level 3 BLAS
 */
#define F77_chemm      BLAS_GLOBAL(chemm,  CHEMM)
#define F77_cherk      BLAS_GLOBAL(cherk,  CHERK)
#define F77_cher2k     BLAS_GLOBAL(cher2k, CHER2K)
#define F77_zhemm      BLAS_GLOBAL(zhemm,  ZHEMM)
#define F77_zherk      BLAS_GLOBAL(zherk,  ZHERK)
#define F77_zher2k     BLAS_GLOBAL(zher2k, ZHER2K)
#define F77_sgemm      BLAS_GLOBAL(sgemm,  SGEMM)
#define F77_ssymm      BLAS_GLOBAL(ssymm,  SSYMM)
#define F77_ssyrk      BLAS_GLOBAL(ssyrk,  SSYRK)
#define F77_ssyr2k     BLAS_GLOBAL(ssyr2k, SSYR2K)
#define F77_strmm      BLAS_GLOBAL(strmm,  STRMM)
#define F77_strsm      BLAS_GLOBAL(strsm,  STRSM)
#define F77_dgemm      BLAS_GLOBAL(dgemm,  DGEMM)
#define F77_dsymm      BLAS_GLOBAL(dsymm,  DSYMM)
#define F77_dsyrk      BLAS_GLOBAL(dsyrk,  DSYRK)
#define F77_dsyr2k     BLAS_GLOBAL(dsyr2k, DSYR2K)
#define F77_dtrmm      BLAS_GLOBAL(dtrmm,  DTRMM)
#define F77_dtrsm      BLAS_GLOBAL(dtrsm,  DTRSM)
#define F77_cgemm      BLAS_GLOBAL(cgemm,  CGEMM)
#define F77_csymm      BLAS_GLOBAL(csymm,  CSYMM)
#define F77_csyrk      BLAS_GLOBAL(csyrk,  CSYRK)
#define F77_csyr2k     BLAS_GLOBAL(csyr2k, CSYR2K)
#define F77_ctrmm      BLAS_GLOBAL(ctrmm,  CTRMM)
#define F77_ctrsm      BLAS_GLOBAL(ctrsm,  CTRSM)
#define F77_zgemm      BLAS_GLOBAL(zgemm,  ZGEMM)
#define F77_zsymm      BLAS_GLOBAL(zsymm,  ZSYMM)
#define F77_zsyrk      BLAS_GLOBAL(zsyrk,  ZSYRK)
#define F77_zsyr2k     BLAS_GLOBAL(zsyr2k, ZSYR2K)
#define F77_ztrmm      BLAS_GLOBAL(ztrmm,  ZTRMM)
#define F77_ztrsm      BLAS_GLOBAL(ztrsm,  ZTRSM)

#ifdef __cplusplus
extern "C" {
#endif

   void F77_xerbla(FCHAR, void *);
/*
 * Level 1 Fortran Prototypes
 */

/* Single Precision */

   void F77_srot(FINT, float *, FINT, float *, FINT, const float *, const float *);
   void F77_srotg(float *,float *,float *,float *);    
   void F77_srotm( FINT, float *, FINT, float *, FINT, const float *);
   void F77_srotmg(float *,float *,float *,const float *, float *);
   void F77_sswap( FINT, float *, FINT, float *, FINT);
   void F77_scopy( FINT, const float *, FINT, float *, FINT);
   void F77_saxpy( FINT, const float *, const float *, FINT, float *, FINT);
   void F77_sdot_sub(FINT, const float *, FINT, const float *, FINT, float *);
   void F77_sdsdot_sub( FINT, const float *, const float *, FINT, const float *, FINT, float *);
   void F77_sscal( FINT, const float *, float *, FINT);
   void F77_snrm2_sub( FINT, const float *, FINT, float *);
   void F77_sasum_sub( FINT, const float *, FINT, float *);
   void F77_isamax_sub( FINT, const float * , FINT, FINT2);

/* Double Precision */

   void F77_drot(FINT, double *, FINT, double *, FINT, const double *, const double *);
   void F77_drotg(double *,double *,double *,double *);    
   void F77_drotm( FINT, double *, FINT, double *, FINT, const double *);
   void F77_drotmg(double *,double *,double *,const double *, double *);
   void F77_dswap( FINT, double *, FINT, double *, FINT);
   void F77_dcopy( FINT, const double *, FINT, double *, FINT);
   void F77_daxpy( FINT, const double *, const double *, FINT, double *, FINT);
   void F77_dswap( FINT, double *, FINT, double *, FINT);
   void F77_dsdot_sub(FINT, const float *, FINT, const float *, FINT, double *);
   void F77_ddot_sub( FINT, const double *, FINT, const double *, FINT, double *);
   void F77_dscal( FINT, const double *, double *, FINT);
   void F77_dnrm2_sub( FINT, const double *, FINT, double *);
   void F77_dasum_sub( FINT, const double *, FINT, double *);
   void F77_idamax_sub( FINT, const double * , FINT, FINT2);

/* Single Complex Precision */

   void F77_cswap( FINT, void *, FINT, void *, FINT);
   void F77_ccopy( FINT, const void *, FINT, void *, FINT);
   void F77_caxpy( FINT, const void *, const void *, FINT, void *, FINT);
   void F77_cswap( FINT, void *, FINT, void *, FINT);
   void F77_cdotc_sub( FINT, const void *, FINT, const void *, FINT, void *);
   void F77_cdotu_sub( FINT, const void *, FINT, const void *, FINT, void *);
   void F77_cscal( FINT, const void *, void *, FINT);
   void F77_icamax_sub( FINT, const void *, FINT, FINT2);
   void F77_csscal( FINT, const float *, void *, FINT);
   void F77_scnrm2_sub( FINT, const void *, FINT, float *);
   void F77_scasum_sub( FINT, const void *, FINT, float *);

/* Double Complex Precision */

   void F77_zswap( FINT, void *, FINT, void *, FINT);
   void F77_zcopy( FINT, const void *, FINT, void *, FINT);
   void F77_zaxpy( FINT, const void *, const void *, FINT, void *, FINT);
   void F77_zswap( FINT, void *, FINT, void *, FINT);
   void F77_zdotc_sub( FINT, const void *, FINT, const void *, FINT, void *);
   void F77_zdotu_sub( FINT, const void *, FINT, const void *, FINT, void *);
   void F77_zdscal( FINT, const double *, void *, FINT);
   void F77_zscal( FINT, const void *, void *, FINT);
   void F77_dznrm2_sub( FINT, const void *, FINT, double *);
   void F77_dzasum_sub( FINT, const void *, FINT, double *);
   void F77_izamax_sub( FINT, const void *, FINT, FINT2);

/*
 * Level 2 Fortran Prototypes
 */

/* Single Precision */

   void F77_sgemv(FCHAR, FINT, FINT, const float *, const float *, FINT, const float *, FINT, const float *, float *, FINT);
   void F77_sgbmv(FCHAR, FINT, FINT, FINT, FINT, const float *,  const float *, FINT, const float *, FINT, const float *, float *, FINT);
   void F77_ssymv(FCHAR, FINT, const float *, const float *, FINT, const float *,  FINT, const float *, float *, FINT);
   void F77_ssbmv(FCHAR, FINT, FINT, const float *, const float *, FINT, const float *, FINT, const float *, float *, FINT);
   void F77_sspmv(FCHAR, FINT, const float *, const float *, const float *, FINT, const float *, float *, FINT);
   void F77_strmv( FCHAR, FCHAR, FCHAR, FINT, const float *, FINT, float *, FINT);
   void F77_stbmv( FCHAR, FCHAR, FCHAR, FINT, FINT, const float *, FINT, float *, FINT);
   void F77_strsv( FCHAR, FCHAR, FCHAR, FINT, const float *, FINT, float *, FINT);
   void F77_stbsv( FCHAR, FCHAR, FCHAR, FINT, FINT, const float *, FINT, float *, FINT);
   void F77_stpmv( FCHAR, FCHAR, FCHAR, FINT, const float *, float *, FINT);
   void F77_stpsv( FCHAR, FCHAR, FCHAR, FINT, const float *, float *, FINT);
   void F77_sger( FINT, FINT, const float *, const float *, FINT, const float *, FINT, float *, FINT);
   void F77_ssyr(FCHAR, FINT, const float *, const float *, FINT, float *, FINT);
   void F77_sspr(FCHAR, FINT, const float *, const float *, FINT, float *); 
   void F77_sspr2(FCHAR, FINT, const float *, const float *, FINT, const float *, FINT,  float *); 
   void F77_ssyr2(FCHAR, FINT, const float *, const float *, FINT, const float *, FINT,  float *, FINT);

/* Double Precision */

   void F77_dgemv(FCHAR, FINT, FINT, const double *, const double *, FINT, const double *, FINT, const double *, double *, FINT);
   void F77_dgbmv(FCHAR, FINT, FINT, FINT, FINT, const double *,  const double *, FINT, const double *, FINT, const double *, double *, FINT);
   void F77_dsymv(FCHAR, FINT, const double *, const double *, FINT, const double *,  FINT, const double *, double *, FINT);
   void F77_dsbmv(FCHAR, FINT, FINT, const double *, const double *, FINT, const double *, FINT, const double *, double *, FINT);
   void F77_dspmv(FCHAR, FINT, const double *, const double *, const double *, FINT, const double *, double *, FINT);
   void F77_dtrmv( FCHAR, FCHAR, FCHAR, FINT, const double *, FINT, double *, FINT);
   void F77_dtbmv( FCHAR, FCHAR, FCHAR, FINT, FINT, const double *, FINT, double *, FINT);
   void F77_dtrsv( FCHAR, FCHAR, FCHAR, FINT, const double *, FINT, double *, FINT);
   void F77_dtbsv( FCHAR, FCHAR, FCHAR, FINT, FINT, const double *, FINT, double *, FINT);
   void F77_dtpmv( FCHAR, FCHAR, FCHAR, FINT, const double *, double *, FINT);
   void F77_dtpsv( FCHAR, FCHAR, FCHAR, FINT, const double *, double *, FINT);
   void F77_dger( FINT, FINT, const double *, const double *, FINT, const double *, FINT, double *, FINT);
   void F77_dsyr(FCHAR, FINT, const double *, const double *, FINT, double *, FINT);
   void F77_dspr(FCHAR, FINT, const double *, const double *, FINT, double *); 
   void F77_dspr2(FCHAR, FINT, const double *, const double *, FINT, const double *, FINT,  double *); 
   void F77_dsyr2(FCHAR, FINT, const double *, const double *, FINT, const double *, FINT,  double *, FINT);

/* Single Complex Precision */

   void F77_cgemv(FCHAR, FINT, FINT, const void *, const void *, FINT, const void *, FINT, const void *, void *, FINT);
   void F77_cgbmv(FCHAR, FINT, FINT, FINT, FINT, const void *,  const void *, FINT, const void *, FINT, const void *, void *, FINT);
   void F77_chemv(FCHAR, FINT, const void *, const void *, FINT, const void *, FINT, const void *, void *, FINT);
   void F77_chbmv(FCHAR, FINT, FINT, const void *, const void *, FINT, const void *, FINT, const void *, void *, FINT);
   void F77_chpmv(FCHAR, FINT, const void *, const void *, const void *, FINT, const void *, void *, FINT);
   void F77_ctrmv( FCHAR, FCHAR, FCHAR, FINT, const void *, FINT, void *, FINT);
   void F77_ctbmv( FCHAR, FCHAR, FCHAR, FINT, FINT, const void *, FINT, void *, FINT);
   void F77_ctpmv( FCHAR, FCHAR, FCHAR, FINT, const void *, void *, FINT);
   void F77_ctrsv( FCHAR, FCHAR, FCHAR, FINT, const void *, FINT, void *, FINT);
   void F77_ctbsv( FCHAR, FCHAR, FCHAR, FINT, FINT, const void *, FINT, void *, FINT);
   void F77_ctpsv( FCHAR, FCHAR, FCHAR, FINT, const void *, void *,FINT);
   void F77_cgerc( FINT, FINT, const void *, const void *, FINT, const void *, FINT, void *, FINT);
   void F77_cgeru( FINT, FINT, const void *, const void *, FINT, const void *, FINT, void *,  FINT);
   void F77_cher(FCHAR, FINT, const float *, const void *, FINT, void *, FINT);
   void F77_cher2(FCHAR, FINT, const void *, const void *, FINT, const void *, FINT, void *, FINT);
   void F77_chpr(FCHAR, FINT, const float *, const void *, FINT, void *);
   void F77_chpr2(FCHAR, FINT, const float *, const void *, FINT, const void *, FINT, void *);

/* Double Complex Precision */

   void F77_zgemv(FCHAR, FINT, FINT, const void *, const void *, FINT, const void *, FINT, const void *, void *, FINT);
   void F77_zgbmv(FCHAR, FINT, FINT, FINT, FINT, const void *,  const void *, FINT, const void *, FINT, const void *, void *, FINT);
   void F77_zhemv(FCHAR, FINT, const void *, const void *, FINT, const void *, FINT, const void *, void *, FINT);
   void F77_zhbmv(FCHAR, FINT, FINT, const void *, const void *, FINT, const void *, FINT, const void *, void *, FINT);
   void F77_zhpmv(FCHAR, FINT, const void *, const void *, const void *, FINT, const void *, void *, FINT);
   void F77_ztrmv( FCHAR, FCHAR, FCHAR, FINT, const void *, FINT, void *, FINT);
   void F77_ztbmv( FCHAR, FCHAR, FCHAR, FINT, FINT, const void *, FINT, void *, FINT);
   void F77_ztpmv( FCHAR, FCHAR, FCHAR, FINT, const void *, void *, FINT);
   void F77_ztrsv( FCHAR, FCHAR, FCHAR, FINT, const void *, FINT, void *, FINT);
   void F77_ztbsv( FCHAR, FCHAR, FCHAR, FINT, FINT, const void *, FINT, void *, FINT);
   void F77_ztpsv( FCHAR, FCHAR, FCHAR, FINT, const void *, void *,FINT);
   void F77_zgerc( FINT, FINT, const void *, const void *, FINT, const void *, FINT, void *, FINT);
   void F77_zgeru( FINT, FINT, const void *, const void *, FINT, const void *, FINT, void *,  FINT);
   void F77_zher(FCHAR, FINT, const double *, const void *, FINT, void *, FINT);
   void F77_zher2(FCHAR, FINT, const void *, const void *, FINT, const void *, FINT, void *, FINT);
   void F77_zhpr(FCHAR, FINT, const double *, const void *, FINT, void *);
   void F77_zhpr2(FCHAR, FINT, const double *, const void *, FINT, const void *, FINT, void *);

/*
 * Level 3 Fortran Prototypes
 */

/* Single Precision */

   void F77_sgemm(FCHAR, FCHAR, FINT, FINT, FINT, const float *, const float *, FINT, const float *, FINT, const float *, float *, FINT);
   void F77_ssymm(FCHAR, FCHAR, FINT, FINT, const float *, const float *, FINT, const float *, FINT, const float *, float *, FINT);
   void F77_ssyrk(FCHAR, FCHAR, FINT, FINT, const float *, const float *, FINT, const float *, float *, FINT);
   void F77_ssyr2k(FCHAR, FCHAR, FINT, FINT, const float *, const float *, FINT, const float *, FINT, const float *, float *, FINT);
   void F77_strmm(FCHAR, FCHAR, FCHAR, FCHAR, FINT, FINT, const float *, const float *, FINT, float *, FINT);
   void F77_strsm(FCHAR, FCHAR, FCHAR, FCHAR, FINT, FINT, const float *, const float *, FINT, float *, FINT);

/* Double Precision */

   void F77_dgemm(FCHAR, FCHAR, FINT, FINT, FINT, const double *, const double *, FINT, const double *, FINT, const double *, double *, FINT);
   void F77_dsymm(FCHAR, FCHAR, FINT, FINT, const double *, const double *, FINT, const double *, FINT, const double *, double *, FINT);
   void F77_dsyrk(FCHAR, FCHAR, FINT, FINT, const double *, const double *, FINT, const double *, double *, FINT);
   void F77_dsyr2k(FCHAR, FCHAR, FINT, FINT, const double *, const double *, FINT, const double *, FINT, const double *, double *, FINT);
   void F77_dtrmm(FCHAR, FCHAR, FCHAR, FCHAR, FINT, FINT, const double *, const double *, FINT, double *, FINT);
   void F77_dtrsm(FCHAR, FCHAR, FCHAR, FCHAR, FINT, FINT, const double *, const double *, FINT, double *, FINT);

/* Single Complex Precision */

   void F77_cgemm(FCHAR, FCHAR, FINT, FINT, FINT, const float *, const float *, FINT, const float *, FINT, const float *, float *, FINT);
   void F77_csymm(FCHAR, FCHAR, FINT, FINT, const float *, const float *, FINT, const float *, FINT, const float *, float *, FINT);
   void F77_chemm(FCHAR, FCHAR, FINT, FINT, const float *, const float *, FINT, const float *, FINT, const float *, float *, FINT);
   void F77_csyrk(FCHAR, FCHAR, FINT, FINT, const float *, const float *, FINT, const float *, float *, FINT);
   void F77_cherk(FCHAR, FCHAR, FINT, FINT, const float *, const float *, FINT, const float *, float *, FINT);
   void F77_csyr2k(FCHAR, FCHAR, FINT, FINT, const float *, const float *, FINT, const float *, FINT, const float *, float *, FINT);
   void F77_cher2k(FCHAR, FCHAR, FINT, FINT, const float *, const float *, FINT, const float *, FINT, const float *, float *, FINT);
   void F77_ctrmm(FCHAR, FCHAR, FCHAR, FCHAR, FINT, FINT, const float *, const float *, FINT, float *, FINT);
   void F77_ctrsm(FCHAR, FCHAR, FCHAR, FCHAR, FINT, FINT, const float *, const float *, FINT, float *, FINT);

/* Double Complex Precision */

   void F77_zgemm(FCHAR, FCHAR, FINT, FINT, FINT, const double *, const double *, FINT, const double *, FINT, const double *, double *, FINT);
   void F77_zsymm(FCHAR, FCHAR, FINT, FINT, const double *, const double *, FINT, const double *, FINT, const double *, double *, FINT);
   void F77_zhemm(FCHAR, FCHAR, FINT, FINT, const double *, const double *, FINT, const double *, FINT, const double *, double *, FINT);
   void F77_zsyrk(FCHAR, FCHAR, FINT, FINT, const double *, const double *, FINT, const double *, double *, FINT);
   void F77_zherk(FCHAR, FCHAR, FINT, FINT, const double *, const double *, FINT, const double *, double *, FINT);
   void F77_zsyr2k(FCHAR, FCHAR, FINT, FINT, const double *, const double *, FINT, const double *, FINT, const double *, double *, FINT);
   void F77_zher2k(FCHAR, FCHAR, FINT, FINT, const double *, const double *, FINT, const double *, FINT, const double *, double *, FINT);
   void F77_ztrmm(FCHAR, FCHAR, FCHAR, FCHAR, FINT, FINT, const double *, const double *, FINT, double *, FINT);
   void F77_ztrsm(FCHAR, FCHAR, FCHAR, FCHAR, FINT, FINT, const double *, const double *, FINT, double *, FINT);

#ifdef __cplusplus
}
#endif

#endif /*  CBLAS_F77_H */
