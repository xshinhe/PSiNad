# Develop Models {#dev_models}

User can add extensions to Kernel class with customized implementation.
```python
import libpykids
import typing

class PyKernel_Test(libpykids.Kernel):
    def __init__(self, name: str):
        libpykids.Kernel.__init__(self, name)
    
    def read_param_impl(self, PM: libpykids.Param):
        self._ndof = PM.get_int('N')
        self._edof = PM.get_int('F')
        self._kcoeff = PM.get_double('k')
        # ...
        return
    
    def init_data_impl(self, DS: libpykids.DataSet):
        self._x_init = DS.numpy('init.x')
        self._p_init = DS.numpy('init.p')
        self._x = DS.numpy('integrator.x')
        self._p = DS.numpy('integrator.p')

        self._pes = DS.numpy('model.vpes')
        self._grad = DS.numpy('model.grad')
        self._hess = DS.numpy('model.hess')
        self._V = DS.numpy('model.V')
        self._dV = DS.numpy('model.dV')

        return

    def init_calc_impl(self, stat: int):
        self._x_init[:] = np.zeros((self._ndof))
        self._p_init[:] = np.zeros((self._ndof))
        self._x[:] = self._x_init[:]
        self._p[:] = self._p_init[:]
        return
    
    def exec_kernel_impl(self, stat: int):
        print(f'hello here {self.name()}!')
        return 0

```