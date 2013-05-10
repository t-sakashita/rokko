FF = mpifrt
FC = mpifrt
CC = mpifcc

# Optimization flags for Fujitsu compiler
CCFLAG = -Kfast,OMP,parallel -X9 -Am -Cpp

SCALAPACK = -L$(HOME)/lib -lscalapack -SSL2BLAMP
BLAS = 

INCPATH = -I/usr/local/include

INCFLAG = $(INCPATH)

LIB = Eigen_s
AR     = ar
RANLIB = ranlib
RM     = rm -f

# FOR PAPI AVAILABLE
LDFLAG  = $(CCFLAG) $(SCALAPACK) $(BLAS)
OBJS = cstab_get_optdim.o communication_s.o eigen_dc.o dlaed3x.o eigen_s.o \
       machine_epsilon.o matrix_adjust_s.o pdlaed2x.o pdlaed3x.o \
       eigen_trd_2update.o eigen_trd_Au.o eigen_trd_v.o eigen_trd_u.o \
       eigen_trd_2update_local.o eigen_trd_2update_local_special.o \
       eigen_trd_2update_v.o eigen_trd_d_preserve_restore.o \
       eigen_trd_init_finalize.o eigen_tbk.o eigen_tbk_main.o \
       eigen_tbk_step.o eigen_trd.o CSTAB_mem.o
OBJSM= cstab_get_optdim.o communication_s.o eigen_dc.o dlaed3x.o eigen_s.o \
       ev_test_2D.o machine_epsilon.o main.o matrix_set.o matrix_adjust_s.o \
       pdlaed2x.o pdlaed3x.o eigen_trd_2update.o eigen_trd_Au.o \
       eigen_trd_v.o eigen_trd_u.o eigen_trd_2update_local.o \
       eigen_trd_2update_local_special.o eigen_trd_2update_v.o \
       eigen_trd_d_preserve_restore.o eigen_trd_init_finalize.o eigen_tbk.o \
       eigen_tbk_main.o eigen_tbk_step.o eigen_trd.o \
       CSTAB_mem.o

all : lib$(LIB).a

lib$(LIB).a: $(OBJS)
	$(RM) $@
	$(AR) r $@ $(OBJS)
	$(RANLIB) $@

a.out : $(OBJSM)
	$(FF) -o a.out $(CCFLAG) $(OBJSM) $(LDFLAG) -L. -l$(LIB)

main.o : main.f
	$(FF) -c main.f $(CCFLAG) $(INCFLAG)
eigen_real_symmetric.o : eigen_real_symmetric.f
	$(FF) -c eigen_real_symmetric.f $(CCFLAG) $(INCFLAG)
eigen_s.o : eigen_s.f
	$(FF) -c eigen_s.f $(CCFLAG) $(INCFLAG)
eigen_trd.o : eigen_trd.f
	$(FF) -c eigen_trd.f $(CCFLAG) $(INCFLAG)
eigen_trd_2update.o : eigen_trd_2update.f
	$(FF) -c eigen_trd_2update.f $(CCFLAG) $(INCFLAG)
eigen_trd_Au.o : eigen_trd_Au.f
	$(FC) -c eigen_trd_Au.f $(CCFLAG) $(INCFLAG)
eigen_trd_v.o : eigen_trd_v.f
	$(FC) -c eigen_trd_v.f $(CCFLAG) $(INCFLAG)
eigen_trd_u.o : eigen_trd_u.f
	$(FC) -c eigen_trd_u.f $(CCFLAG) $(INCFLAG)
eigen_trd_2update_local.o : eigen_trd_2update_local.f
	$(FC) -c eigen_trd_2update_local.f $(CCFLAG) $(INCFLAG)
eigen_trd_2update_local_special.o : eigen_trd_2update_local_special.f
	$(FC) -c eigen_trd_2update_local_special.f $(CCFLAG) $(INCFLAG)
eigen_trd_2update_v.o : eigen_trd_2update_v.f
	$(FC) -c eigen_trd_2update_v.f $(CCFLAG) $(INCFLAG)
eigen_trd_d_preserve_restore.o : eigen_trd_d_preserve_restore.f
	$(FC) -c eigen_trd_d_preserve_restore.f $(CCFLAG) $(INCFLAG)
eigen_trd_init_finalize.o : eigen_trd_init_finalize.f
	$(FC) -c eigen_trd_init_finalize.f $(CCFLAG) $(INCFLAG)
eigen_tbk.o : eigen_tbk.f
	$(FF) -c eigen_tbk.f $(CCFLAG)  $(INCFLAG)
eigen_tbk_main.o : eigen_tbk_main.f
	$(FF) -c eigen_tbk_main.f $(CCFLAG)  $(INCFLAG)
eigen_tbk_step.o : eigen_tbk_step.f
	$(FF) -c eigen_tbk_step.f $(CCFLAG)  $(INCFLAG)
machine_epsilon.o : machine_epsilon.f
	$(FF) -c machine_epsilon.f -O0 $(CCFLAG) -O0 $(INCFLAG)
matrix_set.o: matrix_set.f
	$(FF) -c matrix_set.f $(CCFLAG) $(INCFLAG)
matrix_adjust_s.o: matrix_adjust_s.f
	$(FF) -c matrix_adjust_s.f $(CCFLAG) $(INCFLAG)
ev_test_2D.o: ev_test_2D.f
	$(FF) -c ev_test_2D.f $(CCFLAG) $(INCFLAG)
communication_s.o : communication_s.f
	$(FF) -c communication_s.f $(CCFLAG) $(INCFLAG)
eigen_dc.o : eigen_dc.f
	$(FC) -c eigen_dc.f $(CCFLAG) $(INCFLAG)
dlaed3x.o : dlaed3x.f
	$(FC) -c dlaed3x.f $(CCFLAG) $(INCFLAG)
pdlaed2x.o : pdlaed2x.f
	$(FC) -c pdlaed2x.f $(CCFLAG) $(INCFLAG)
pdlaed3x.o : pdlaed3x.f
	$(FC) -c pdlaed3x.f $(CCFLAG) $(INCFLAG)

clean :
	-\rm -f *.o *.mod a.out lib$(LIB).a

eigen_tbk.o: CSTAB.h trd.h
eigen_tbk_main.o: CSTAB.h trd.h
eigen_tbk_step.o: CSTAB.h trd.h
eigen_trd.o: CSTAB.h trd.h
communication_s.o: trd.h
eigen_dc.o: trd.h
ev_test_2D.o: trd.h
main.o: trd.h
matrix_adjust_s.o: trd.h
eigen_trd_2update.o: trd.h
eigen_trd_Au.o: trd.h trd_au.h
eigen_trd_v.o: trd.h
eigen_trd_u.o: trd.h
eigen_trd_2update_local.o: trd.h
eigen_trd_2update_local_special.f.o: trd.h
eigen_trd_2update_v.o: trd.h
eigen_trd_d_preserve_restore.o: trd.h
eigen_trd_init_finalize.o: trd.h trd_au.h

