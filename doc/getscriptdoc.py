#!/usr/bin/env python

import imp

script = imp.load_source('tamkin-driver', '../scripts/tamkin-driver.py')
with open('reference/tamkin-driver.rst', 'w') as f:
    f.write(script.__doc__)
