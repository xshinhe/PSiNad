# Python APIs {#apis_python}

## Construction of Param

```python
import binder.libpykids as kids

dir(kids)
dir(kids.Param)
print(dir(kids.DataSet))

PM = kids.Param('''
{
    \"model\": \"MNDO_Intef\",
}
''', kids.Param.fromString)
```

## Construction of DataSet

```python
import libpykids as kids

print(dir(kids))
print(dir(kids.Param))
print(dir(kids.DataSet))

DS = kids.DataSet()

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


## Embeded KIDS in Thirdpart Library



