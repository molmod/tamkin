#!/bin/bash
# This is a very simplistic install scipt. Use with care!

if [ -n $1 ] && [ "$1" = "--system" ]; then
  ./uninstall.sh --system
  python setup.py install --prefix=/usr/local
  ./cleanfiles.sh
else
  ./uninstall.sh
  python setup.py install --home=$HOME
  ./cleanfiles.sh
  echo "Don't forget to add 'export PYTHONPATH=\$HOME/lib/python' to your .bashrc file."
fi

