# Kernel_Customized.py
import pykids
import pykids.libpykids_v1.data as DA

class Kernel_Customized_X(pykids.lib.Kernel):
    # ......
    
    def __init__(self, name: str):
        pykids.lib.Kernel.__init__(self, name)
        self._name = name
        print(self.getName())
        
    def getName(self):
        return self._name
    
    def setInputParam_impl(self, PM: pykids.lib.Param):
        return
    
    # import pykids.lib.data as DA
    def setInputDataSet_impl(self, DS: pykids.lib.DataSet):
        self.dt = DS.numpy_tpl_real(DA.iter.dt)
        self.x  = DS.numpy_tpl_real(DA.integrator.x)
        self.p  = DS.numpy_tpl_real(DA.integrator.p)
        self.m  = DS.numpy_tpl_real(DA.integrator.m)
        self.minv_dt = 1 / self.m * self.dt[0]
        return
    
    def initializeKernel_impl(self, stat: pykids.lib.Status):
        return stat
    
    def executeKernel_impl(self, stat: pykids.lib.Status):
        self.x[:] += self.p[:] * self.minv_dt
        return stat
    
    