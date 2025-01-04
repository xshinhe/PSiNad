from pyscf import gto, scf, mcscf, nac
'''
mol = gto.M(
        atom = atoms_string,
        basis = basis,
        verbose=0
)

mf = scf.ROHF()
mf.kernel()
'''

mol = gto.M(atom='N 0 0 0; N 0 0 1', basis='ccpvdz')
mf = scf.RHF(mol).run()
mc = mcscf.CASSCF(mf, 2, 2).state_average([0.5, 0.5]).run()
mc_nac = nac.sacasscf.NonAdiabaticCouplings(mc)
# mc_nac = mc.nac_method() # Also valid
mc_nac.kernel(state=(0,1), use_etfs=False)

