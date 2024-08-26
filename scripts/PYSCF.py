from pyscf import gto, scf, mcscf
'''
mol = gto.M(
        atom = atoms_string,
        basis = basis,
        verbose=0
)

mf = scf.ROHF()
mf.kernel()
'''

mol = gto.M(atom='N 0 0 0; N 0 0 1.1', verbose=0)
mc_grad_scanner1 = mcscf.CASSCF(scf.RHF(mol), 4, 4).nuc_grad_method().as_scanner(1)
mc_grad_scanner2 = mcscf.CASSCF(scf.RHF(mol), 4, 4).nuc_grad_method().as_scanner(2)
etot, grad = mc_grad_scanner1(gto.M(atom='N 0 0 0; N 0 0 1.1'))
print(etot, grad)
etot, grad = mc_grad_scanner2(gto.M(atom='N 0 0 0; N 0 0 1.1'))
print(etot, grad)
#etot, grad = mc_grad_scanner(gto.M(atom='N 0 0 0; N 0 0 1.5'))
#print(etot, grad)

