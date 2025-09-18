# Python APIs {#apis_python}


## Installation

```
pip install pybind11 Cython
cd PSiNaD/build
make
make PythonInstall
```

## Construction of Param

```python
import binder.libpypsnd as psnd

dir(psnd)
dir(psnd.Param)
print(dir(psnd.DataSet))

PM = psnd.Param('''
{
    \"model\": \"MNDO_Intef\",
}
''', psnd.Param.fromString)
```

## Construction of DataSet

```python
import libpypsnd as psnd

print(dir(psnd))
print(dir(psnd.Param))
print(dir(psnd.DataSet))

DS = psnd.DataSet()

# define variables in DataSet
DS._def('a.a', (10,), dtype='int', doc='this is int array')
DS._def('a.b', (10,), dtype='real', doc='this is real array')
DS._def('a.b', (10,), dtype='real', doc='this is real2 array') # it only allows to define once
DS._def('a.c', (10,), dtype='complex', doc='this is complex array')
DS._def('a.d.e', (10,), dtype='bool', doc='this is bool array')

print(DS)

DS._undef('a.a')
b = DS.numpy('a.b')
b[1] = 3.14

print(DS)

print(DS.help())
print(DS.help('a.b'))
print(DS.help('a.d'))
```

## Construction of Kernel


## Embeded PSND in Thirdpart Library



