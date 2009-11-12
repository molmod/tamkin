#!/bin/bash
echo "Cleaning python code in $(pwd) and subdirectories"
for file in `find lib test examples | egrep "(\.py$)|(\.f90$)|(^scripts/tr-)"`; do
  echo Cleaning ${file}
  sed -i -e $'s/\t/    /' ${file}
  sed -i -e $'s/[ \t]\+$//' ${file}
done
