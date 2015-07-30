#!/bin/sh

set_prefix() {
  PREFIX_DEF="$HOME/rokko"
  BUILD_DIR_DEF="$HOME/build"
  SOURCE_DIR_DEF="$HOME/source"
  SUDO_DEF="/usr/bin/sudo"
  SUDO_DEF="/usr/bin/sudo"

  if [ -f "$HOME/.rokkoinstaller" ]; then
    source $HOME/.rokkoinstaller
  fi

  if [ -z "$PREFIX_ROKKO" ]; then
    PREFIX_ROKKO="$PREFIX_DEF"
  fi
  echo "PREFIX_ROKKO=$PREFIX_ROKKO"
  if [ -d "$PREFIX_ROKKO" ]; then :; else
    echo "Fatal: target directory $PREFIX_ROKKO does not exist!"
    exit 127
  fi

  if [ -z "$SUDO" ]; then
    RES=$(touch $PREFIX_ROKKO/.rokkoinstaller.tmp > /dev/null 2>&1; echo $?; rm -f $PREFIX_ROKKO/.rokkoinstaller.tmp)
    if [ $RES = 0 ]; then
      SUDO=
    else
      SUDO="$SUDO_DEF"
    fi
  fi
  echo "SUDO=$SUDO"

  if [ -z "$BUILD_DIR" ]; then
    BUILD_DIR="$BUILD_DIR_DEF"
  fi
  echo "BUILD_DIR=$BUILD_DIR"
  if [ -d "$BUILD_DIR" ]; then :; else
    echo "Fatal: target directory $BUILD_DIR does not exist!"
    exit 127
  fi
  RES=$(touch $BUILD_DIR/.rokkoinstaller.tmp > /dev/null 2>&1; echo $?; rm -f $BUILD_DIR/.rokkoinstaller.tmp)
  if [ $RES = 0 ]; then :; else
    echo "Fatal: have no permission to write in build directory $BUILD_DIR"
    exit 127
  fi

  if [ -z "$SOURCE_DIR" ]; then
    SOURCE_DIR="$SOURCE_DIR_DEF"
  fi
  echo "SOURCE_DIR=$SOURCE_DIR"

  return 0
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
