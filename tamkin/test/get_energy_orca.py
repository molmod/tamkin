# 
#    get_energy_orca.py: a python function to read and return the optimized energy
#    from ORCA output file
#    Usage: energy=get_energy(orcf)
#    Copyright (C) 2020  Mark E. Fuller
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#    Mark E. Fuller: mark.e.fuller@gmx.de


#setup terminal output later:
#    get_energy_orca.py  Copyright (C) 2020  Mark E. Fuller
#    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
#    This is free software, and you are welcome to redistribute it
#    under certain conditions; type `show c' for details.


################################################################################

def get_energy_orca(orcf):

    #extract energy from orca output file
    E0 = [] # Fill this list with E0 lines.

    # Get the data from the single-point calculation.
    outfile = open(orcf,'r')
    outlines = outfile.readlines()
    outfile.close()

    for (l,line) in enumerate(outlines):
        if line.startswith('FINAL SINGLE POINT ENERGY'):
            E0.append(line)

    #print (float(E0[-1].split()[-1]))
    
    q = float(E0[-1].split()[-1])
    
    return q
