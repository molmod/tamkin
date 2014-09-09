#!/usr/bin/env python

import imp
import codecs

script = imp.load_source('tamkin-driver', '../scripts/tamkin-driver.py')
with codecs.open('reference/tamkin-driver.rst', 'w', 'utf-8') as f:
    f.write(script.__doc__)
