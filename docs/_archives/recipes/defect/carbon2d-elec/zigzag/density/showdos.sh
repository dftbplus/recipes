grep "Fermi level" detailed.out
plotxy --xlabel "Energy [eV]" --ylabel "DOS" dos.dat pdos.C1.dat pdos.C2.dat pdos.C3.dat pdos.H.dat &
