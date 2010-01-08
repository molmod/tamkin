"""
This script finds the overlap between the 7th mode and the 7th mode
and writes the result to a file delta.overlaps.7-7.csv.
First line/column contains the frequencies.
Next column contains the overlap (still to be squared).

Script does this for every delta.*overlap.csv file present in the directory.

usage:  python %progr
"""


def get_freqs_and_overlap_7_7(filename):
    f = file(filename)
    print "-"*20
    print "Reading file...", filename

    # read 7th freq on first line (first element is 0, so is 8th freq)
    for line in f:
        words = line.split(";")
        #freq1 = float(words[7])  # !!!!!!!!!!!!!
        freq1 = float(words[-1])
        break

    # get 7th freq in first column (we already skipped first line)
    # and overlap 7th-7th
    count = 0
    for line in f:
        count +=1
        if count==7:
            words = line.split(";")
            freq2 = float(words[0])
            #overlap = float(words[7])  # !!!!!!!!!
            overlap = float(words[-1])
            break

    f.close()
    print "freq1: ", freq1
    print "freq2: ", freq2
    print "overlap: ", overlap

    return freq1, freq2, overlap


import sys, glob

freqs1 = []
freqs2 = []
overlaps = []
filenames = glob.glob("delta.*overlap.csv")

filenames.sort()
print "Considered files: ", filenames
for filename in filenames:
    #y = str(x)
    #if x < 10:
    #    y = "0"+str(x)
    freq1, freq2, overlap = get_freqs_and_overlap_7_7(filename)
    freqs1.append(freq1)
    freqs2.append(freq2)
    overlaps.append(overlap)

conv = 7.25163277859109E-07  # conversion factor such that frequencies are printed in 1/cm

filename_out = "delta.overlaps.7-7.csv"
f = file(filename_out,"w+")
for i in range(len(freqs1)):
    print >> f, filenames[i]+";"+str(freqs1[i]/conv)+";"+str(freqs2[i]/conv)+";"+str(overlaps[i])
f.close()
print "file written:", filename_out


