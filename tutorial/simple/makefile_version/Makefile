# Makefile for rokko/tutorial/simple

ROKKO_ROOT_DIR = ${HOME}/rokko_install
include $(ROKKO_ROOT_DIR)/share/rokko/include.mk

TARGET  = frank_mpi
SRCS    = frank_mpi.cpp
OBJS    = $(SRCS:.cpp=.o)

CC  = $(CMAKE_CXX_COMPILER)
LIBS    = $(ROKKO_LIBS)
CPPFLAGS = $(ROKKO_CPPFLAGS)
LDFLAGS = $(ROKKO_LDFLAGS)

# リンク時にライブラリを追加
$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS) $(LIBDIR) $(LIBS)
 
# ソースをコンパイルするときにインクルードディレクトリを指定
$(OBJS): $(SRCS)
	$(CC) $(INCDIR) $(CPPFLAGS) -c $(SRCS)



