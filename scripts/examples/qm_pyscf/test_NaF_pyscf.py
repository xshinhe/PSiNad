import numpy as np
from pyscf import lib
from pyscf.lib import logger
from pyscf.fci import direct_spin1
from pyscf.mcscf import newton_casscf
from pyscf.grad import casscf as casscf_grad
from pyscf.grad import sacasscf as sacasscf_grad
from pyscf.nac.sacasscf import NonAdiabaticCouplings
from functools import reduce
from pyscf import gto, scf, mcscf
from scipy import linalg


def calc_fix(R:float):
# build molecule
    mol = gto.M(atom = '''
    Na 0 0 0; F 0 0 %.8f;
    '''%(R)
    , basis='cc-pVTZ'
    , output='mol.log')

    mf = scf.RHF(mol).run()
    mc = mcscf.CASSCF(mf, 8, 8).fix_spin_(ss=0).state_average ([0.5,0.5]).run(conv_tol=1e-10)

#openmolcas_energies = np.array ([-7.85629118, -7.72175252])
#print ("energies:",mc.e_states)
#print ("disagreement w openmolcas:", np.around (mc.e_states-openmolcas_energies, 8))

    for i in range(len(mc.e_states)):
        print("ROOT%i = %12.8f"%(i, mc.e_states[i]))


    mc_grad = mc.nuc_grad_method ()
    for i in range(len(mc.e_states)):
        de = mc_grad.kernel (state=i)
        #print('GRAD ON STATE %d'%i)
        #print('A               X           Y           Z')
        for k in range(len(de)):
            print('%s  %12.8f  %12.8f  %12.8f'%(mol.atom_symbol(k), de[k][0], de[k][1], de[k][2] ))

    mc_nacs = NonAdiabaticCouplings(mc)
    first_ind, second_ind = np.triu_indices(len(mc.e_states), k=1)
    for i, j in zip(first_ind, second_ind):
        
        nac = mc_nacs.kernel (state=(i,j), use_etfs=True)
        #print('NAC ON STATE %d %d'%(i,j))
        #print('A               X           Y           Z')
        for k in range(len(nac)):
            print('%s  %12.8f  %12.8f  %12.8f'%(mol.atom_symbol(k), nac[k][0], nac[k][1], nac[k][2] ))

if __name__ == '__main__':
    calc_fix(1.5)

