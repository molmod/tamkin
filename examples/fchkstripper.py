#!/usr/bin/env python

# This script takes a list of Gaussian fchk files as arguments and will strip
# all sections that are irrelevant for TAMkin. Also the trailing white space on
# each line is stripped. This is used to reduce the size of the example fchk
# files in TAMkin.


keep_fields = set([
    "Multiplicity", "Total Energy", "Atomic numbers",
    "Current cartesian coordinates", "Real atomic weights",
    "Cartesian Gradient", "Cartesian Force Constants",
])

def strip(fn):
    # a list with all lines that we'd like to keep
    lines = []
    # load the good lines
    f = file(fn, "r")
    busy = False
    keep = False
    for line in f:
        line = line.rstrip()
        if len(lines) < 2:
            # keep the first two lines
            lines.append(line)
        else:
            if busy:
                if not line.startswith(" "):
                    keep = False
                    busy = False
                if keep:
                    lines.append(line)
            if not busy:
                busy = line[47:49] == "N="
                field = line[:43].strip()
                if field in keep_fields:
                    lines.append(line)
                    keep = True
                else:
                    keep = False
    f.close()
    # print stuff back into the same file
    f = file(fn, "w")
    for line in lines:
        print >> f, line
    f.close()


if __name__ == "__main__":
    import sys
    fns_fchk = sys.argv[1:]
    for fn in fns_fchk:
        if fn.endswith(".fchk"):
            print "Stripping", fn
            strip(fn)
        else:
            print "Skipping", fn, "(wrong extension)"
