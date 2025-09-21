# Kernel_Customized.py
import pypsnd
import pypsnd.libpypsnd_v1.data as DA

class Kernel_Customized_X(pypsnd.lib.Kernel):
    # ......
    
    def __init__(self, name: str):
        pypsnd.lib.Kernel.__init__(self, name)
        self._name = name
        print(self.getName())
        
    def getName(self):
        return self._name
    
    def setInputParam_impl(self, PM: pypsnd.lib.Param):
        return
    
    # import pypsnd.lib.data as DA
    def setInputDataSet_impl(self, DS: pypsnd.lib.DataSet):
        self.dt = DS.numpy_tpl_real(DA.iter.dt)
        self.x  = DS.numpy_tpl_real(DA.integrator.x)
        self.p  = DS.numpy_tpl_real(DA.integrator.p)
        self.m  = DS.numpy_tpl_real(DA.integrator.m)
        self.minv_dt = 1 / self.m * self.dt[0]
        return
    
    def initializeKernel_impl(self, stat: pypsnd.lib.Status):
        return stat
    
    def executeKernel_impl(self, stat: pypsnd.lib.Status):
        self.x[:] += self.p[:] * self.minv_dt
        return stat
    
    