!#include "f2c.h"
!#include <malloc.h>
!----
!
!extern void
!dsyevd_(char *jobz, char *uplo, integer *n, doublereal *
!a, integer *lda, doublereal *w, doublereal *work, integer *lwork, 
!integer *iwork, integer *liwork, integer *info);
!
!int lapack_eigen2_(int *n, int *n, int *hbw, int *id, double *d,
!double *e, int *lde, double *q, int *ldq) {

      subroutine lapack_eigen2(n_global, n, hbw, id, d, e, lde, q, ldq)

      integer :: n_global, n, hbw, id, lde, ldq
      real(8) :: d(*), e(lde,*), q(ldq,*)

      integer :: i, j, info, lwork, liwork
      real(8) :: temp
      real(8), pointer :: work(:)
      integer, pointer :: iwork(:)
      character*1 :: jobu, jobvt

      do j=1,n
         do i=1,n
            q(i,j)=0.0d0
         enddo
      enddo

      do i=1,n
         q(i,i)=d(i+id-1)
      enddo

      do j=1,hbw
         do i=1,n-j
            q(i,i+j)=e(i+id-1,j)
            q(i+j,i)=e(i+id-1,j)
         enddo
      enddo

      ldqs = ldq

      jobu = 'v'; jobvt= 'u'

      lwork = -1
      liwork = -1

      call dsyevd(jobu, jobvt, n, q, ldq, d(1+id-1),
     &            temp, lwork, i, liwork, info)

      lwork  = int(temp)
      liwork = i;
!---- print*,lwork, liwork

      allocate(work(lwork), iwork(liwork))

      call dsyevd(jobu, jobvt, n, q, ldq, d(1+id-1),
     &             work, lwork, iwork, liwork, info)

      deallocate(work,iwork)

      end subroutine

