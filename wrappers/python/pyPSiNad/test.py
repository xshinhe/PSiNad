from __future__ import absolute_import
from pyPSiNad import Kernel, Param, DataSet, Status
import typing

class PyKernel_Customized(Kernel):
	def __init__(self, name: str):
		Kernel.__init__(self, name)
		self._name = name
		print(self.getName())
	
	def getName(self):
		return self._name

	def setpInputParam_impl(self, PM: Param):
		# read parameters from Param
		return
	
	def setInputDataSet_impl(self, DS: DataSet):
		# fetch some data from DataSet
		print(f'read param {DS}!')
		return
	
	def executeKernel_impl(self, stat: Status):
		# perfrom some algorithm with PM & DS
		return 0

if __name__ == '__main__':
	print("test")
	