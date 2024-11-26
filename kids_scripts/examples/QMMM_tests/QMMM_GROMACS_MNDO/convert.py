import parmed as pmd

amber = pmd.load_file('amber.prmtop', 'amber.inpcrd')
#Save a GROMACS topology and GRO files
amber.save('gromacs.top')
amber.save('gromacs.gro')

