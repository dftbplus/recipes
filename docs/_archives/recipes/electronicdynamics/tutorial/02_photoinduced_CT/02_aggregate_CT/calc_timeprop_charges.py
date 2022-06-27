#!/usr/bin/env python3

'''
Charge per fragment of DFTB+ TD data
'''

import sys
import optparse
import numpy as np

USAGE = """usage: %prog -l ii:jj,ll:mm 

Reads output from TD calculation with external laser and produces net charges per fragment
(subtracting value at time = 0).

Needs qsvst.dat file present in working directory."""

def main():
    parser = optparse.OptionParser(usage=USAGE)
    parser.add_option("-l", "--list", type='str', action="store", dest="at_list",
            help="list of atom indices starting from 1 (initial and final index separated by colons, ranges separated by commas)") 

    (options, args) = parser.parse_args()

    print(options.at_list)
    list_ranges = options.at_list.split(',')
    nfrag = len(list_ranges)          #number of fragments   
    idx_ini = np.zeros(nfrag, dtype=int)
    idx_end = np.zeros(nfrag, dtype=int)
    for ii, ir in enumerate(list_ranges):
        ini_end = ir.split(':')
        idx_ini[ii] = int(ini_end[0]) - 1 # because first index is 0
        idx_end[ii] = int(ini_end[1]) - 1

    qsdata = np.genfromtxt('qsvst.dat')
    time = qsdata[:,0]
    natoms = qsdata.shape[1] - 2
    print ('Number of atoms = {}'.format(natoms))
    if any(idx_ini > natoms) or any(idx_end) > natoms+1:
        print('List ranges above total number of atoms = {}'.format(natoms))
        sys.exit()
    
    for ii in range(nfrag):       #sum of the atomic charges inside each fragment
        print ('Fragment {}: from atom {} till {}'.format(ii, idx_ini[ii]+1, idx_end[ii]+1))
        q_fragment = qsdata[:,2+idx_ini[ii]:2+idx_end[ii]+1].sum(axis=1)
        np.savetxt("charge_frag%s.dat"%(ii+1), np.column_stack((time, q_fragment-q_fragment[0])))

if __name__ == "__main__":
    main()

