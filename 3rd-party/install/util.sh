#!/bin/sh

set_prefix() {
  PREFIX_OPT_DEF="$HOME/opt"
  PREFIX_ROKKO_DEF="$HOME/rokko"
  SUDO=

  hostname -f > /dev/null 2>&1
  if [ $? = 0 ]; then
    HOSTNAME=$(hostname -f)
  else
    HOSTNAME=$(hostname)
  fi

  if [ -d /opt/rokko ]; then
    SUDO="sudo LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
    PREFIX_ROKKO_DEF="/opt/rokko"
    if [ -d /opt/local ]; then
      PREFIX_OPT_DEF="/opt/local"
    else
      PREFIX_OPT_DEF="/opt/alps"
    fi
  elif [ -d /opt/nano/rokko ]; then
    SUDO="sudo LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
    PREFIX_ROKKO_DEF="/opt/nano/rokko"
    if [ -d /opt/local ]; then
      PREFIX_OPT_DEF="/opt/local"
    else
      PREFIX_OPT_DEF="/opt/nano/alps"
    fi
  elif [ -d /opt/nano/alps/rokko ]; then
    SUDO=
    PREFIX_ROKKO_DEF="/opt/nano/alps/rokko"
    PREFIX_OPT_DEF="/opt/nano/alps"
  fi

  # for Mac OS X
  if [ $(uname) = Darwin ]; then
    if [ -d /opt/rokko ]; then
      PREFIX_OPT_DEF="/opt/alps"
      PREFIX_ROKKO_DEF="/opt/rokko"
    fi
  fi

  # for camphor.kudpc.kyoto-u.ac.jp
  if [ -d /LARGE0/hp120237 ]; then
    PREFIX_OPT_DEF="/LARGE0/hp120237/opt"
    PREFIX_ROKKO_DEF="/LARGE0/hp120237/rokko"
  fi

  # for k.aics.riken.jp
  if [ -d /opt/spire/alps ]; then
    PREFIX_OPT_DEF="/opt/spire/alps"
    PREFIX_ROKKO_DEF="/opt/spire/alps/rokko"
  fi

  # for kashiwa.issp.u-tokyo.ac.jp
  if [ -d /home/issp/materiapps ]; then
    PREFIX_OPT_DEF="/home/issp/materiapps/opt"
    PREFIX_ROKKO_DEF="/home/issp/materiapps/rokko"
  fi

  # for maki.issp.u-tokyo.ac.jp
  if [[ ! -z `echo "$HOSTNAME" | egrep "^maki.\.fx10hpc$"` ]]; then
    PREFIX_OPT_DEF="/global/app/materiapps/opt"
    PREFIX_ROKKO_DEF="/global/app/materiapps/rokko"
  fi

  # for oakleaf-fx.cc.u-tokyo.ac.jp
  if [[ ! -z `echo "$HOSTNAME" | egrep "^oakleaf-fx.*$"` ]]; then
    PREFIX_OPT_DEF="/group/gc25/share/opt"
    PREFIX_ROKKO_DEF="/group/gc25/share/rokko"
  fi

  if [ -z "$PREFIX_OPT" ]; then
    PREFIX_OPT="$PREFIX_OPT_DEF"
  fi
  if [ -z "$PREFIX_ROKKO" ]; then
    PREFIX_ROKKO="$PREFIX_ROKKO_DEF"
  fi
  echo "PREFIX_OPT = $PREFIX_OPT"
  echo "PREFIX_ROKKO = $PREFIX_ROKKO"
  echo "SUDO = $SUDO"
  return 0
}

set_build_dir() {
  if [ -z "$BUILD_DIR" ]; then
    BUILD_DIR="$HOME/build"
  fi
  mkdir -p "$BUILD_DIR"
}

check() {
  "$@"
  result=$?
  if [ $result -ne 0 ]; then
    echo "Failed: $@" >&2
    exit $result
  fi
  return 0
}

start_info() {
  echo "Start: $(date) on $(hostname)"
}

finish_info() {
  echo "Finish: $(date)"
}
