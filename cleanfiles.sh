#!/bin/bash
for i in $(find tamkin test | egrep "\.pyc$|\.py~$|\.pyc~$|\.bak$") ; do rm -v ${i}; done

rm -vr python-build-stamp-* 

rm -v MANIFEST
rm -vr dist
rm -vr build
rm -vr doc/_build
rm -vr doctrees
rm -v scripts/tamkin-driverc
rm -vr tamkin.egg-info
