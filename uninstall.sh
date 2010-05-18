#!/bin/bash
# This is a very simplistic uninstall scipt. Use with care!

if [ -n $1 ] && [ "$1" = "--system" ]; then
  rm -vr /usr/local/lib/python*/site-packages/tamkin
else
  rm -vr $HOME/lib/python/tamkin
fi
