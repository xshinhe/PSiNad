import pyPSiNad

class Kernel_Example_Test(pyPSiNad.Kernel):
    def __init__(self, name: str):
        pyPSiNad.Kernel.__init__(self, name)
        self._name = name
        print(self.getName())
        
    def getName(self):
        return self._name
    
    def setpInputParam_impl(self, PM):
        print(f'read param {PM}!')
        return
    
    def setInputDataSet_impl(self, DS):
        print(f'read param {DS}!')
        return
    
    def executeKernel_impl(self, stat:int):
        print(f'hello here {self.name()}!')
        return 0