from __future__ import absolute_import
from pyPSiNad import Kernel, Model, Param, DataSet, Status, data
import typing
import numpy as np

class PyModel_SBtest(Model):
	def __init__(self, name: str):
		Model.__init__(self, name)
		self._name = name
	
	def getName(self):
		return self._name

	def setInputParam_impl(self, PM: Param):
		print(PM)
		self._P = PM.get_int(['solver.P'])
		self._N = PM.get_int(['model.N'])
		self._F = PM.get_int(['model.F'])

		omegac = PM.get_real(['model.omegac'])
		lamb = PM.get_real(['model.lamb'])

		z = np.arange(self._N) + 0.5
		self._w = - omegac * np.log( z / self._N)
		self._c = np.sqrt(2.0 * lamb /  self._N) * self._w

		# read parameters from Param
		return
	
	def setInputDataSet_impl(self, DS: DataSet):
		DS._def("model.vpes", (self._P,), 'real', '')
		DS._def("model.grad", (self._P, self._N), 'real', '')
		DS._def("model.V", (self._P, self._F, self._F), 'real', '')
		DS._def("model.dV", (self._P, self._N, self._F, self._F), 'real', '')

		self._x = DS.numpy_tpl_real(data.integrator.x)
		self._vpes = DS.numpy_tpl_real(data.model.vpes)
		self._grad = DS.numpy_tpl_real(data.model.grad)
		self._V = DS.numpy_tpl_real(data.model.V)
		self._dV  = DS.numpy_tpl_real(data.model.dV)
		return
	
	def executeKernel_impl(self, stat: Status):
		vpes = np.zeros((self._P))
		grad = np.zeros((self._P, self._N))
		V = np.zeros((self._P, self._F, self._F))
		dV = np.zeros((self._P, self._N, self._F, self._F))

		vpes = 0.5 * np.sum( self._w **2 * self._x**2 )
		grad = self._w**2 * self._x
		V[:,0,0] = epsilon + self._c * self._x
		V[:,1,1] =-epsilon - self._c * self._x
		V[:,0,1] = delta
		V[:,1,0] = delta

		dV[:,0,0] = self._c 
		dV[:,1,1] =-self._c

		# perfrom some algorithm with PM & DS
		return 0
