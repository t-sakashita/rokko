/*
 * cblas_test.h
 * Written by Keita Teranishi
 */
#ifndef CBLAS_TEST_H
#define CBLAS_TEST_H
#include "cblas.h"

#define  TRUE           1
#define  PASSED         1
#define  TEST_ROW_MJR	1

#define  FALSE          0
#define  FAILED         0
#define  TEST_COL_MJR	0

#define  INVALID       -1
#define  UNDEFINED     -1

typedef struct { float real; float imag; } CBLAS_TEST_COMPLEX;
typedef struct { double real; double imag; } CBLAS_TEST_ZOMPLEX;

#include "cblas_mangling.h"

#define F77_xerbla     BLAS_GLOBAL(xerbla, XERBLA)
/*
 * Level 1 BLAS
 */
#define F77_srotg      BLAS_GLOBAL(srotgtest,  SROTG)
#define F77_srotmg     BLAS_GLOBAL(srotmgtest, SROTMG)
#define F77_srot       BLAS_GLOBAL(srottest,   SROT)
#define F77_srotm      BLAS_GLOBAL(srotmtest,  SROTM)
#define F77_drotg      BLAS_GLOBAL(drotgtest,  DROTG)
#define F77_drotmg     BLAS_GLOBAL(drotmgtest, DROTMG)
#define F77_drot       BLAS_GLOBAL(drottest,   DROT)
#define F77_drotm      BLAS_GLOBAL(drotmtest,  DROTM)
#define F77_sswap      BLAS_GLOBAL(sswaptest,  SSWAP)
#define F77_scopy      BLAS_GLOBAL(scopytest,  SCOPY)
#define F77_saxpy      BLAS_GLOBAL(saxpytest,  SAXPY)
#define F77_isamax     BLAS_GLOBAL(isamaxtest, ISAMAX)
#define F77_dswap      BLAS_GLOBAL(dswaptest,  DSWAP)
#define F77_dcopy      BLAS_GLOBAL(dcopytest,  DCOPY)
#define F77_daxpy      BLAS_GLOBAL(daxpytest,  DAXPY)
#define F77_idamax     BLAS_GLOBAL(idamaxtest, IDAMAX)
#define F77_cswap      BLAS_GLOBAL(cswaptest,  CSWAP)
#define F77_ccopy      BLAS_GLOBAL(ccopytest,  CCOPY)
#define F77_caxpy      BLAS_GLOBAL(caxpytest,  CAXPY)
#define F77_icamax     BLAS_GLOBAL(icamaxtest, ICAMAX)
#define F77_zswap      BLAS_GLOBAL(zswaptest,  ZSWAP)
#define F77_zcopy      BLAS_GLOBAL(zcopytest,  ZCOPY)
#define F77_zaxpy      BLAS_GLOBAL(zaxpytest,  ZAXPY)
#define F77_izamax     BLAS_GLOBAL(izamaxtest, IZAMAX)
#define F77_sdot       BLAS_GLOBAL(sdottest,   SDOT)
#define F77_ddot       BLAS_GLOBAL(ddottest,   DDOT)
#define F77_dsdot      BLAS_GLOBAL(dsdottest,  DSDOT)
#define F77_sscal      BLAS_GLOBAL(sscaltest,  SSCAL)
#define F77_dscal      BLAS_GLOBAL(dscaltest,  DSCAL)
#define F77_cscal      BLAS_GLOBAL(cscaltest,  CSCAL)
#define F77_zscal      BLAS_GLOBAL(zscaltest,  ZSCAL)
#define F77_csscal     BLAS_GLOBAL(csscaltest, CSSCAL)
#define F77_zdscal     BLAS_GLOBAL(zdscaltest, ZDSCAL)
#define F77_cdotu      BLAS_GLOBAL(cdotutest,  CDOTU)
#define F77_cdotc      BLAS_GLOBAL(cdotctest,  CDOTC)
#define F77_zdotu      BLAS_GLOBAL(zdotutest,  ZDOTU)
#define F77_zdotc      BLAS_GLOBAL(zdotctest,  ZDOTC)
#define F77_snrm2      BLAS_GLOBAL(snrm2test,  SNRM2)
#define F77_sasum      BLAS_GLOBAL(sasumtest,  SASUM)
#define F77_dnrm2      BLAS_GLOBAL(dnrm2test,  DNRM2)
#define F77_dasum      BLAS_GLOBAL(dasumtest,  DASUM)
#define F77_scnrm2     BLAS_GLOBAL(scnrm2test, SCNRM2)
#define F77_scasum     BLAS_GLOBAL(scasumtest, SCASUM)
#define F77_dznrm2     BLAS_GLOBAL(dznrm2test, DZNRM2)
#define F77_dzasum     BLAS_GLOBAL(dzasumtest, DZASUM)
#define F77_sdsdot     BLAS_GLOBAL(sdsdottest, SDSDOT)
/*
 * Level 2 BLAS
 */
#define F77_s2chke     BLAS_GLOBAL(cs2chke, CS2CHKE)
#define F77_d2chke     BLAS_GLOBAL(cd2chke, CD2CHKE)
#define F77_c2chke     BLAS_GLOBAL(cc2chke, CC2CHKE)
#define F77_z2chke     BLAS_GLOBAL(cz2chke, CZ2CHKE)
#define F77_ssymv      BLAS_GLOBAL(cssymv,  CSSYMV)
#define F77_ssbmv      BLAS_GLOBAL(cssbmv,  CSSBMV)
#define F77_sspmv      BLAS_GLOBAL(csspmv,  CSSPMV)
#define F77_sger       BLAS_GLOBAL(csger,   CSGER)
#define F77_ssyr       BLAS_GLOBAL(cssyr,   CSSYR)
#define F77_sspr       BLAS_GLOBAL(csspr,   CSSPR)
#define F77_ssyr2      BLAS_GLOBAL(cssyr2,  CSSYR2)
#define F77_sspr2      BLAS_GLOBAL(csspr2,  CSSPR2)
#define F77_dsymv      BLAS_GLOBAL(cdsymv,  CDSYMV)
#define F77_dsbmv      BLAS_GLOBAL(cdsbmv,  CDSBMV)
#define F77_dspmv      BLAS_GLOBAL(cdspmv,  CDSPMV)
#define F77_dger       BLAS_GLOBAL(cdger,   CDGER)
#define F77_dsyr       BLAS_GLOBAL(cdsyr,   CDSYR)
#define F77_dspr       BLAS_GLOBAL(cdspr,   CDSPR)
#define F77_dsyr2      BLAS_GLOBAL(cdsyr2,  CDSYR2)
#define F77_dspr2      BLAS_GLOBAL(cdspr2,  CDSPR2)
#define F77_chemv      BLAS_GLOBAL(cchemv,  CCHEMV)
#define F77_chbmv      BLAS_GLOBAL(cchbmv,  CCHBMV)
#define F77_chpmv      BLAS_GLOBAL(cchpmv,  CCHPMV)
#define F77_cgeru      BLAS_GLOBAL(ccgeru,  CCGERU)
#define F77_cgerc      BLAS_GLOBAL(ccgerc,  CCGERC)
#define F77_cher       BLAS_GLOBAL(ccher,   CCHER)
#define F77_chpr       BLAS_GLOBAL(cchpr,   CCHPR)
#define F77_cher2      BLAS_GLOBAL(ccher2,  CCHER2)
#define F77_chpr2      BLAS_GLOBAL(cchpr2,  CCHPR2)
#define F77_zhemv      BLAS_GLOBAL(czhemv,  CZHEMV)
#define F77_zhbmv      BLAS_GLOBAL(czhbmv,  CZHBMV)
#define F77_zhpmv      BLAS_GLOBAL(czhpmv,  CZHPMV)
#define F77_zgeru      BLAS_GLOBAL(czgeru,  CZGERU)
#define F77_zgerc      BLAS_GLOBAL(czgerc,  CZGERC)
#define F77_zher       BLAS_GLOBAL(czher,   CZHER)
#define F77_zhpr       BLAS_GLOBAL(czhpr,   CZHPR)
#define F77_zher2      BLAS_GLOBAL(czher2,  CZHER2)
#define F77_zhpr2      BLAS_GLOBAL(czhpr2,  CZHPR2)
#define F77_sgemv      BLAS_GLOBAL(csgemv,  CSGEMV)
#define F77_sgbmv      BLAS_GLOBAL(csgbmv,  CSGBMV)
#define F77_strmv      BLAS_GLOBAL(cstrmv,  CSTRMV)
#define F77_stbmv      BLAS_GLOBAL(cstbmv,  CSTBMV)
#define F77_stpmv      BLAS_GLOBAL(cstpmv,  CSTPMV)
#define F77_strsv      BLAS_GLOBAL(cstrsv,  CSTRSV)
#define F77_stbsv      BLAS_GLOBAL(cstbsv,  CSTBSV)
#define F77_stpsv      BLAS_GLOBAL(cstpsv,  CSTPSV)
#define F77_dgemv      BLAS_GLOBAL(cdgemv,  CDGEMV)
#define F77_dgbmv      BLAS_GLOBAL(cdgbmv,  CDGBMV)
#define F77_dtrmv      BLAS_GLOBAL(cdtrmv,  CDTRMV)
#define F77_dtbmv      BLAS_GLOBAL(cdtbmv,  CDTBMV)
#define F77_dtpmv      BLAS_GLOBAL(cdtpmv,  CDTPMV)
#define F77_dtrsv      BLAS_GLOBAL(cdtrsv,  CDTRSV)
#define F77_dtbsv      BLAS_GLOBAL(cdtbsv,  CDTBSV)
#define F77_dtpsv      BLAS_GLOBAL(cdtpsv,  CDTPSV)
#define F77_cgemv      BLAS_GLOBAL(ccgemv,  CCGEMV)
#define F77_cgbmv      BLAS_GLOBAL(ccgbmv,  CCGBMV)
#define F77_ctrmv      BLAS_GLOBAL(cctrmv,  CCTRMV)
#define F77_ctbmv      BLAS_GLOBAL(cctbmv,  CCTBMV)
#define F77_ctpmv      BLAS_GLOBAL(cctpmv,  CCTPMV)
#define F77_ctrsv      BLAS_GLOBAL(cctrsv,  CCTRSV)
#define F77_ctbsv      BLAS_GLOBAL(cctbsv,  CCTBSV)
#define F77_ctpsv      BLAS_GLOBAL(cctpsv,  CCTPSV)
#define F77_zgemv      BLAS_GLOBAL(czgemv,  CZGEMV)
#define F77_zgbmv      BLAS_GLOBAL(czgbmv,  CZGBMV)
#define F77_ztrmv      BLAS_GLOBAL(cztrmv,  CZTRMV)
#define F77_ztbmv      BLAS_GLOBAL(cztbmv,  CZTBMV)
#define F77_ztpmv      BLAS_GLOBAL(cztpmv,  CZTPMV)
#define F77_ztrsv      BLAS_GLOBAL(cztrsv,  CZTRSV)
#define F77_ztbsv      BLAS_GLOBAL(cztbsv,  CZTBSV)
#define F77_ztpsv      BLAS_GLOBAL(cztpsv,  CZTPSV)
/*
 * Level 3 BLAS
 */
#define F77_s3chke     BLAS_GLOBAL(cs3chke, CS3CHKE)
#define F77_d3chke     BLAS_GLOBAL(cd3chke, CD3CHKE)
#define F77_c3chke     BLAS_GLOBAL(cc3chke, CC3CHKE)
#define F77_z3chke     BLAS_GLOBAL(cz3chke, CZ3CHKE)
#define F77_chemm      BLAS_GLOBAL(cchemm,  CCHEMM)
#define F77_cherk      BLAS_GLOBAL(ccherk,  CCHERK)
#define F77_cher2k     BLAS_GLOBAL(ccher2k, CCHER2K)
#define F77_zhemm      BLAS_GLOBAL(czhemm,  CZHEMM)
#define F77_zherk      BLAS_GLOBAL(czherk,  CZHERK)
#define F77_zher2k     BLAS_GLOBAL(czher2k, CZHER2K)
#define F77_sgemm      BLAS_GLOBAL(csgemm,  CSGEMM)
#define F77_ssymm      BLAS_GLOBAL(cssymm,  CSSYMM)
#define F77_ssyrk      BLAS_GLOBAL(cssyrk,  CSSYRK)
#define F77_ssyr2k     BLAS_GLOBAL(cssyr2k, CSSYR2K)
#define F77_strmm      BLAS_GLOBAL(cstrmm,  CSTRMM)
#define F77_strsm      BLAS_GLOBAL(cstrsm,  CSTRSM)
#define F77_dgemm      BLAS_GLOBAL(cdgemm,  CDGEMM)
#define F77_dsymm      BLAS_GLOBAL(cdsymm,  CDSYMM)
#define F77_dsyrk      BLAS_GLOBAL(cdsyrk,  CDSYRK)
#define F77_dsyr2k     BLAS_GLOBAL(cdsyr2k, CDSYR2K)
#define F77_dtrmm      BLAS_GLOBAL(cdtrmm,  CDTRMM)
#define F77_dtrsm      BLAS_GLOBAL(cdtrsm,  CDTRSM)
#define F77_cgemm      BLAS_GLOBAL(ccgemm,  CCGEMM)
#define F77_csymm      BLAS_GLOBAL(ccsymm,  CCSYMM)
#define F77_csyrk      BLAS_GLOBAL(ccsyrk,  CCSYRK)
#define F77_csyr2k     BLAS_GLOBAL(ccsyr2k, CCSYR2K)
#define F77_ctrmm      BLAS_GLOBAL(cctrmm,  CCTRMM)
#define F77_ctrsm      BLAS_GLOBAL(cctrsm,  CCTRSM)
#define F77_zgemm      BLAS_GLOBAL(czgemm,  CZGEMM)
#define F77_zsymm      BLAS_GLOBAL(czsymm,  CZSYMM)
#define F77_zsyrk      BLAS_GLOBAL(czsyrk,  CZSYRK)
#define F77_zsyr2k     BLAS_GLOBAL(czsyr2k, CZSYR2K)
#define F77_ztrmm      BLAS_GLOBAL(cztrmm,  CZTRMM)
#define F77_ztrsm      BLAS_GLOBAL(cztrsm,  CZTRSM)

void get_transpose_type(char *type, enum CBLAS_TRANSPOSE *trans);
void get_uplo_type(char *type, enum CBLAS_UPLO *uplo);
void get_diag_type(char *type, enum CBLAS_DIAG *diag);
void get_side_type(char *type, enum CBLAS_SIDE *side);

#endif /* CBLAS_TEST_H */
