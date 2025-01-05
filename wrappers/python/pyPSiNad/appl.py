from pyPSiNad import Kernel, Model, Param, DataSet, Status, data
import copy 

class PopulationTransfer(Kernel):
	def __init__(self, M: Model):
		self._name = 'PopulationTransfer'
		Kernel.__init__(self, self._name)

		self._model = M
	
	def getName(self):
		return self._name

	def setpInputParam_impl(self, PM: Param):
		self._occ0 = PM.get_int('solver.occ', 0)
		return
	
	def setInputDataSet_impl(self, DS: DataSet):
		self._K0 = DS.numpy_tpl_complex(data.integrator.K0)
		self._K0_init = copy.deepcopy(self._K0)
		M,F,_ = np.shape(self._K0)
		self._pop = DS._def('appl.pop', (F), 'int', '')
		print(f'K0 = {K0}!')
		return
	
	def executeKernel_impl(self, stat: Status):
		self._pop[:] = np.einsum('m,mii->i', self._K0_init[:, self._occ, self._occ], self._K0[:,:,:]) / M
		return 0

