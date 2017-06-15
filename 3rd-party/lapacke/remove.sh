#!/bin/sh

FUNC3r="\
geqrfp getrf2 gges3 ggev3 gghd3 lagge lagsy latms orcsd2by1 \
"
FUNC3c="\
geqrfp getrf2 gges3 ggev3 gghd3 heswapr hetri2 hetri2x lagge laghe lagsy latms sysv_aa sytrs2 sytrs_aa uncsd2by1 \
"

FUNC32r=" \
gesvxx posvxx gerfsx porfsx sysvxx gbsvxx syrfsx gbrfsx \
geequb syequb poequb gbequb \
"
FUNC32c=" \
gesvxx posvxx gerfsx porfsx sysvxx gbsvxx syrfsx gbrfsx hesvxx herfsx \
geequb syequb poequb gbequb heequb \
"

FUNC33r=" \
orcsd orbdb bbcsd lapmr lartgp lartgs \
sytrs2 syconv syswapr sytri2 sytri2x \
"
FUNC33c=" \
uncsd unbdb bbcsd lapmr \
hetrs2 syconv syswapr sytri2 sytri2x \
"

FUNC34r=" \
geqrt gemqrt geqrt2 \
geqrt3 \
tpmqrt tpqrt tpqrt2 tprfb \
"
FUNC34c=" \
geqrt gemqrt geqrt2 \
geqrt3 \
tpmqrt tpqrt tpqrt2 tprfb \
"

FUNC35r=" \
lasyf_rook sycon_rook sysv_rook sytf2_rook sytrf_rook sytri_rook sytrs_rook \
"
FUNC35c=" \
hecon_rook sytrf_rook hesv_rook sytri_rook hetf2_rook sytrs_rook hetrf_rook hetri_rook hetrs_rook lahef_rook lasyf_rook sycon_rook sysv_rook sytf2_rook \
"

FUNC36r=" \
ggsvd3 ggsvp3 \
bdsvdx gesvdx \
potrf2 \
"
FUNC36c=" \
ggsvd3 ggsvp3 \
gesvdx \
gejsv gesvj gsvj0 gsvj1 \
potrf2 \
"

FUNC37r="
gelq gelqt gelqt3 gemlq gemlqt gemqr geqr getsls lamswlq lamtsqr laswlq latsqr tplqt tplqt2 tpmlqt \
sysv_aa sytrf_aa sytrs_aa lasyf_aa \
sytf2_rk sytf2_rk lasyf_rk lasyf_rk sytrf_rk sytrf_rk sytrs_3 sytrs_3 sycon_3 sycon_3 sytri_3 sytri_3 sytri_3x sytri_3x sysv_rk sysv_rk \
sb2st_kernels sbev_2stage sbevd_2stage sbevx_2stage syev_2stage syevd_2stage syevr_2stage syevx_2stage sygv_2stage sytrd_2stage sytrd_sb2st sytrd_sy2sb larfy \
"
FUNC37c="\
gelq gelqt gelqt3 gemlq gemlqt gemqr geqr getsls lamswlq lamtsqr laswlq latsqr tplqt tplqt2 tpmlqt \
hesv_aa hetrf_aa hetrs_aa lahef_aa \
sytf2_rk hetf2_rk lasyf_rk lahef_rk sytrf_rk hetrf_rk sytrs_3 hetrs_3 sycon_3 hecon_3 sytri_3 hetri_3 sytri_3x hetri_3x sysv_rk hesv_rk \
hb2st_kernels hbev_2stage hbevd_2stage hbevx_2stage heev_2stage heevd_2stage heevr_2stage heevx_2stage hegv_2stage hetrd_2stage hetrd_hb2st hetrd_he2hb larfy \
gejsv gesvj gsvj0 gsvj1 \
"

for func in $FUNC3r $FUNC32r $FUNC33r $FUNC34r $FUNC35r $FUNC36r $FUNC37r; do
  for p in s d; do
    rm -f src/lapacke_${p}${func}.c src/lapacke_${p}${func}_work.c
  done
done
for func in $FUNC3c $FUNC32c $FUNC33c $FUNC34c $FUNC35c $FUNC36c $FUNC37c; do
  for p in c z; do
    rm -f src/lapacke_${p}${func}.c src/lapacke_${p}${func}_work.c
  done
done
