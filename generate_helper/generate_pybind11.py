#!/bin/python

import json
import glob
import re
import sys
import os

root=sys.argv[1]
config_file = sys.argv[2]

if os.path.exists(root+'/config.json'):
    config_file = root + '/config.json'


with open(config_file, 'r', encoding='utf-8') as load_f:
    global data
    data = json.load(load_f)

decro_list = ['const', 'static', 'inline', 'virtual', 'constexpr']
type_list = ['void', 'bool', 'int', 'double', 'kids_real', 'kids_complex',
             'std::string', 'std::map', 'enum', 'class', 'namespace']

def parse_decro_type_flag(line):
    decro1, decro2 = '', ''

    # more compact
    res = re.search('<(.*?)>', line, re.S)
    if res:
        i1, i2 = res.span()
        old = line[i1:i2]
        line = re.sub(old, '<__>', line)

    tmpsx = re.split('\(|;|\{', line)
    tmps0 = tmpsx[0].split(' ')
    while '' in tmps0: tmps0.remove('')

    if '(' not in line:
        se_type = '_var'
    else:
        if len(tmps0) == 1:
            se_type = '_init'
        else:
            se_type = '_fun'
    xxx = ''
    for t in tmps0:
        if t in decro_list:
            if xxx == '':
                decro1 += t + ' '
            else:
                decro2 += t + ' '
        else:
            xxx += t + ' '
    xxx2 = re.split('=', xxx)[0]
    tmps0 = re.split(' |,|;', xxx2)
    while '' in tmps0: tmps0.remove('')

    typeinfo = ''
    if len(tmps0) > 1:
        typeinfo = tmps0[0]
        flags = tmps0[1:]
    else:
        flags = tmps0[:]

    if res:
        typeinfo = re.sub('<__>', old, typeinfo)

    if se_type == '_fun' and typeinfo == '': #bugs: assume that ~final all with virtual!
        se_type = '_final'

    return decro1, typeinfo, decro2, flags, se_type

def remove_comment(txt):
    '''
        this function remove c-style comment of in the text
    '''
    pn1 = re.compile('//.*')
    pn2 = re.compile('\/\*(?:[^\*]|\*+[^\/\*])*\*+\/')
    txt = re.sub(pn1, '', txt)
    txt = re.sub(pn2, '', txt)
    return txt

def safe_split_comma(line_in):
    '''
        this function split comma, such as
        <double, std::string>, <std::map<int, double>>, <<int>, 2>
    '''
    line = line_in
    _dd = []
    res = re.search('(<[^<]*?>)', line, re.S)
    while res:
        i1, i2 = res.span()
        _dd += [ line[i1:i2] ]

        line = line[:i1] + "I__%d__I"%len(_dd) + line[i2:]
        res = re.search('(<[^<]*?>)', line, re.S)

    line = re.sub(',', '@', line)

    while _dd:
        line = line.replace("I__%d__I"%len(_dd), _dd[-1])
        _dd = _dd[:-1]

    return line.split('@')

def unique_type(typi):
    typi = re.sub(' +?', ' ', typi)
    typi = re.sub('(\w) (\w)', r'\1I__0__I\2', typi)
    typi = re.sub(' ', '', typi)
    typi = re.sub('I__0__I', ' ', typi)
    return typi

def parse_type_and_name(term):
    '''
        this function split (unique) typeinfo & variable name.
        passed for:
        double* abc=0;
        double* abc_d;
        double* abc_d[];
        double* abc_d[2];
        double*&abc_d[];
        const std::sting &abc_d[];
    '''

    # remove by symbol '=' (assignment) & ';' (end sentence)
    _tmps = re.split(';|=', term)[0].strip()

    # split array literal, such as 'A[2]' to 'A' and '[2]'
    liter = ''
    res =  re.search('\[.*?$', _tmps, re.S)
    if res:
        i1, i2 = res.span()
        _tmps, liter = _tmps[:i1], _tmps[i1:i2]

    res =  re.search('[\w]*?$', _tmps, re.S)
    i1, i2 = res.span()
    name = _tmps[i1:i2]
    typi = _tmps[:i1] + liter

    # make type to be unique
    typi = unique_type(typi)

    return typi, name

def parse_type_and_name_multiple(term):
    '''
        this function split multiple typeinfo & variable name.
        passed for:
        "const int *a, b, **c;"
    '''

    _tmps = safe_split_comma(re.split(';', term)[0]);

    type0, name0 = parse_type_and_name(_tmps[0])
    arg_type = [type0]
    arg_name = [name0]

    # find base type for 0-th variable
    res =  re.search('(\W)*?$', type0, re.S)
    if res:
        i1, i2 = res.span()
        type0 = type0[:i1]

    # generate derived type for other variable
    for i in range(1, len(_tmps)):
        res =  re.search('^(\W)*?\w', _tmps[i], re.S)
        typex = type0
        namex = _tmps[i]
        if res:
            i1, i2 = res.span()
            typex += _tmps[i][i1:i2-1]
            namex = _tmps[i][i2-1:]
            typex = unique_type(typex)
        arg_type += [typex]
        arg_name += [namex]

    return arg_type, arg_name


def dict_append(d, keys, c):
    ks = keys.split('.')
    L = len(ks)
    ds = [{}]*L
    ds[0] = d
    for i in range(L-1):
        if ks[i] in ds[i]:
            ds[i+1] = ds[i][ks[i]]
        else:
            ds[i][ks[i]] = {}
            ds[i+1] = ds[i][ks[i]]

    if ks[-1] not in ds[-1]:
        ds[-1][ks[-1]] = c
    else:
        if(isinstance(c, list)):
            ds[-1][ks[-1]] += c
        if(isinstance(c, set)):
            ds[-1][ks[-1]] |= c
        if(isinstance(c, str)):
            ds[-1][ks[-1]] += '\n' + c
        if(isinstance(c, dict)):
            ds[-1][ks[-1]] = dict( ds[-1][ks[-1]].items() + c.items())

def sweep_scope(lines, istart):
    cnt_1 = 0 # count for ( and )
    cnt_2 = 0 # count for { and }
    max_1 = 0
    max_2 = 0
    for i in range(istart, len(lines)):
        cnt_1 += lines[i].count('(')
        cnt_2 += lines[i].count('{')
        max_1 = max(max_1, cnt_1)
        max_2 = max(max_2, cnt_2)
        cnt_1 -= lines[i].count(')')
        cnt_2 -= lines[i].count('}')
        if cnt_1 == 0 and cnt_2 == 0: # and ';' in lines[i]: # may not
            return i, max_1, max_2
    return istart, 0, 0

def parse_argument(content):
    tokens = safe_split_comma( (content.strip().split(')')[0]).split('(')[-1] )
    arg_type = []
    arg_name = []
    for term in tokens:
        _1, _2 = parse_type_and_name(term)
        arg_type += [_1]
        arg_name += [_2]
    return arg_type, arg_name

def parse_scope(info, kstring, lines, init_attr):
    # print('\n'.join(lines))
    # exit(-1)

    il = 0
    scope_attr = init_attr
    for _ in lines:
        if il >= len(lines): break

        line = lines[il]
        if line.strip() == '' or line[0] == '#':
            il += 1
            continue

        # check attribute
        if 'public:' in line:
            scope_attr = 'public'
            il += 1
            continue
        if 'protected:' in line:
            scope_attr = 'protected'
            il += 1
            continue
        if 'private:' in line:
            scope_attr = 'private'
            il += 1
            continue

        if 'DEFINE_POINTER' in line: # bind for container (at first token off| not confused with fun!)
            tmps = re.split(' |,|\(|\)', line)
            while '' in tmps: tmps.remove('')

            flag = tmps[2]
            dict_append(info, kstring+'._bind', [
                {flag: { '_type':'EigMX<%s>'%tmps[1], '_raw': line.strip(), '_attr': scope_attr}}
            ])
            il += 1
            continue

        # sweep a scope and judge it to be _var, _fun or _init
        il2, m1, m2 = sweep_scope(lines, il)
        content = '\n'.join(lines[il:il2+1])

        de1, typi, de2, flags, se_type = parse_decro_type_flag(line)

        if m1 == 0 and se_type != '_var':
            print('Error'*10)
        else:
            arg_type, arg_name = [], []
            if se_type == '_fun' or se_type == '_init':
                arg_type, arg_name = parse_argument(content)
                # print(arg_type, arg_name)
            for flag in flags:
                dict_append(info, kstring+'.'+ se_type, [{
                    flag : {'_type': typi, '_raw': content.strip(), '_attr': scope_attr
                    ,'_arg_type': ', '.join(arg_type)
                    ,'_arg_name': ', '.join(arg_name)
                    }}
                ])
        il = il2 + 1

def file_parse(fn):

    f = open(fn, 'r',encoding='utf-8')
    txt = ''.join(f.readlines())
    txt = remove_comment(txt)

    info = {
        "name": fn, # fn.split('/')[-1],
        "deps": re.findall('#include "(.*?)"', txt, re.S)
    }
    kstring =''
    scope_attr = ''

    lines = txt.split('\n'); il = 0
    for _ in lines:
        if il >= len(lines): break

        line = lines[il]
        if line.strip() == '' or line[0] == '#':
            il += 1
            continue

        terms = line.split(' ')
        if terms[0] == 'using':
            while ';' not in lines[il]:
                il += 1
            il += 1
            continue

        if '{' in terms:
            if terms[0] == 'namespace':
                flag = terms[1]
                kstring = kstring + '.' + flag
                dict_append(info, kstring, {'_type': 'namespace', '_var': []})
                il2, _, _ = sweep_scope(lines, il)
                parse_scope(info, kstring, lines[il+1:il2], 'public')
                kstring = '.'.join(kstring.split('.')[:-1])
                il = il2 + 1
                continue

            if terms[0] == 'class' or terms[0] == 'struct' :
                flag = terms[1]
                superclass = ''
                if terms[2] == ':':
                    superclass = terms[4]
                kstring = kstring + '.' + flag
                dict_append(info, kstring, {'_type': terms[0], '_superclass': superclass,
                                            '_var':[], '_fun':[], '_init':[], '_bind':[],
                                           })
                il2, _, _ = sweep_scope(lines, il)
                parse_scope(info, kstring, lines[il+1:il2], 'private')
                kstring = '.'.join(kstring.split('.')[:-1])
                il = il2 + 1
                continue

        # sweep a scope and judge it to be _var, _fun or _init
        il2, m1, m2 = sweep_scope(lines, il)
        content = '\n'.join(lines[il:il2+1])

        de1, typi, de2, flags, se_type = parse_decro_type_flag(line)

        if m1 == 0 and se_type != '_var':
            print('Error'*10)
        else:
            for flag in flags:
                dict_append(info, kstring+'.'+ se_type, [{
                    flag : {'_type': typi, '_raw': content.strip(), '_attr': scope_attr}}
                ])
        il = il2 + 1
    return info

# info = file_parse('/home/public/hexin/share/github/opendf/solvers/solvers_md/traj.h')
# print(json.dumps(info, indent=4))
# exit(-1)
######################################################################

objs = []
for i in data['models']:
    list1 = glob.glob(os.path.abspath(root+'/'+i))
    for j in list1:
        objs += [ '%s.h'%j[:-4] ]

for i in data['solvers']:
    list1 = glob.glob(os.path.abspath(root+'/'+i))
    for j in list1:
        objs += [ '%s.h'%j[:-4] ]

# print(objs)
# objs = [
# '/home/public/hexin/share/github/opendf/models/nad_forcefield/systembath.h',
# '/home/public/hexin/share/github/opendf/models/forcefieldbase.h',
# '/home/public/hexin/share/github/opendf/solvers/solvers_md/traj.h'
# '/home/public/hexin/share/github/opendf/solvers/solvers_nad/solver_mmd.h'
# ]
# exit(-1)

complete=False
infos = {'incl': [], 'deps':[], '': {}}
while not complete:
    infos1 = {'incl': [], 'deps':[], '': {}}
    for i in objs:
        info = file_parse(i)
        infos1['incl'] += [ info['name'] ]
        for j in info['deps']:
            infos1['deps'] += [ os.path.abspath(
                os.path.dirname(info['name']) + '/' + j
                )]
        infos1[''] = {**infos1[''], **info['']}

    infos['incl'] += infos1['incl']
    infos['deps'] += infos1['deps']
    infos[''] = {**infos[''], **infos1['']}

    s1 = set(infos['deps'])
    s2 = set(infos['incl'])
    s3 = set()
    for i in s1-s2:
        # print('checking deps: ', i)
        if 'opendf/models' in i or 'opendf/solvers' in i:
            s3 = s3 | {i}
    if s3:
        objs = list(s3)
    else:
        complete=True
        # print("****")

# print(json.dumps(infos, sort_keys=False, indent=4))
# print(infos[''])
# exit(-1)


#####################################################


def get_fathers(dd, child):
    if '_superclass' in dd[child] and dd[child]['_superclass'] != '':
        return child + ', ' + get_fathers(dd, dd[child]['_superclass'])
    return child

def creat_trampoline_fun(ns, ns0, dd, coll):
    for f in dd[ns]['_fun']:
        for i in f:
            if f[i]['_type'] not in type_list:
                continue
            if 'virtual' in f[i]['_raw']:
                virtual_in_name = 'PYBIND11_OVERRIDE'
                if ') = 0;' in f[i]['_raw']:
                    virtual_in_name = 'PYBIND11_OVERRIDE_PURE'

                tag = '%s(%s)'%(i, f[i]['_arg_type'])
                if tag in coll:
                    continue
                coll += [tag]

                token = f[i]['_raw'].split(')')[0] + ')'
                token = token.replace('virtual ', '')
                print('''
        %s override {
            %s(
            %s, // return type
            %s, // parent class
            %s, // func name
            %s
            );
        }'''%(token, virtual_in_name,
                    f[i]['_type'],
                    ns,
                    i,
                    f[i]['_arg_name']
                    )
                )
    return coll

def creat_trampoline(ns, dd):
    print('''
    class PyTrampoline_%s : public %s {
        public:
        using %s::%s;'''%(ns, ns, ns, ns)
        )
    coll = []
    coll = creat_trampoline_fun(ns, ns, dd, coll)
    ns1 = ns
    while dd[ns1]['_superclass'] != '':
        ns1 = dd[ns1]['_superclass']
        coll = creat_trampoline_fun(ns1, ns, dd, coll)
    print('    };\n')

# def creat_public_members(ns, ns0, dd, coll):
#     # for f in dd[ns]['_var']:
#     #     for i in f:
#     #         if f[i]['_attr'] == 'protected':
#     #             print('        using %s::%s;'%(ns, i))
#     for f in dd[ns]['_bind']:
#         for i in f:
#             if f[i]['_attr'] == 'protected':
#                 print('        using %s::%s;'%(ns, i))
#                 print('        using %s::%s_eigen_container;'%(ns, i))
# def creat_publicist(ns, dd):
#     print('''
#     class PyPublicist_%s : public %s {
#         public:'''%(ns, ns)
#         )
#     coll = []
#     coll = creat_public_members(ns, ns, dd, coll)
#     # ns1 = ns
#     # while dd[ns1]['_superclass'] != '':
#     #     ns1 = dd[ns1]['_superclass']
#     #     coll = creat_trampoline_fun(ns1, ns, dd, coll)
#     print('    };\n')

def creat_class_init(ns, dd):
    for f in dd[ns]['_init']:
        for i in f:
            print('\n    .def(py::init<%s>())'%f[i]['_arg_type'], end='')

def creat_class_var(ns, dd):
    for f in dd[ns]['_var']:
        for i in f:
            if f[i]['_attr'] == 'public':
                if f[i]['_type'] in ['int', 'double', 'kids_real', 'kids_complex',
                'std::string', 'bool'] and i[0]!='*':
                    print('\n    .def_readwrite("%s", &%s::%s)'%(
                        i, ns, i
                        ), end='')

def creat_class_bind(ns, dd):
    for f in dd[ns]['_bind']:
        for i in f:
            print('\n    .def("ref_%s", &%s::ref_%s, py::return_value_policy::reference_internal)'%(
                i, ns, i
                ), end='')

def creat_class_fun(ns, ns0, dd):
    for f in dd[ns]['_fun']:
        for i in f:
            if f[i]['_type'] not in type_list or 'kids_complex*' in f[i]['_raw']:
                continue

            if i == 'name' and ns != ns0:
                continue
            if 'static' in f[i]['_raw']:
                print('\n    .def_static("%s", &%s::%s)'%
                    (
                    i, ns, i
                    ), end='')
            elif '*' not in f[i]['_arg_type']:
                print('\n    .def("%s", static_cast<%s (%s::*)(%s)>(&%s::%s))'%
                    (
                    i, f[i]['_type'], ns, f[i]['_arg_type'], ns, i
                    ), end='')
            else:
                arg_type, arg_name = parse_argument(f[i]['_raw'])
                arg_name1 = []
                arg_name2 = []
                pairs = []
                for k in range(len(arg_type)):
                    if arg_type[k][-1] == '*':
                        arg_type[k] = 'py::array_t<%s, py::array::c_style | py::array::forcecast>'%arg_type[k][:-1]
                        arg_name1 += [arg_name[k]+'_arr']
                        arg_name2 += [arg_name[k]+'_arr.mutable_data()']
                    else:
                        arg_name1 += [arg_name[k]]
                        arg_name2 += [arg_name[k]]
                    pairs += [arg_type[k] + ' ' + arg_name1[k]]
                print('\n    .def("%s", [](%s& self, %s) {\n            return self.%s(%s); \n        }\n    )'%
                    (
                    i, ns, ', \n        '.join(pairs), i, ', '.join(arg_name2)
                    ), end='')

def creat_class_(ns, dd, field):
    tplname = ns
    if dd[ns]['_superclass'] != '':
        tplname = '%s, %s'%(ns, dd[ns]['_superclass'])
    tplname += ', PyTrampoline_%s'%(ns)
    # clsname = ns.lower()

    print('    py::class_<%s>(%s, "%s", py::dynamic_attr())'
        %(tplname, field, ns), end='')
    creat_class_init(ns, dd)
    creat_class_var(ns, dd)
    creat_class_bind(ns, dd)
    creat_class_fun(ns, ns, dd)

    ns1 = ns
    while dd[ns1]['_superclass'] != '':
        ns1 = dd[ns1]['_superclass']
        creat_class_fun(ns1, ns, dd)

    print(';', end='\n')

print('''
#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
''')
for i in infos['incl']:
    print('#include "%s"'%os.path.relpath(i, root+'/python'))
# exit(-1)

print('''

#include "../utils/definitions.h"

namespace py = pybind11;

// clang-format off
PYBIND11_MODULE(libopendf, m) {

    #include "opendf_phys.bind"

    py::module models_m = m.def_submodule("models");
    py::module solvers_m = m.def_submodule("solvers");

''')

dd = infos['']
created_class = []

def try_creat_class(ns):
    global created_class
    if (ns not in created_class
        ## following are classes with bugs
        and 'SCF' not in ns
        and 'Atomic_BasisSet' not in ns
        and 'PPIMD' not in ns
        ):
        if dd[ns]['_superclass'] != '':
            try_creat_class(dd[ns]['_superclass']) # father must previously created!!
        # get_fathers(dd, ns)
        creat_trampoline(ns, dd)
        if 'Solver' in ns:
            creat_class_(ns, dd, 'solvers_m')
        else:
            creat_class_(ns, dd, 'models_m')

        created_class += [ns]

for ns in dd.keys():
    if '_type' in dd[ns] and dd[ns]['_type'] == 'class':
        try_creat_class(ns)

# for i in infos['deps']:
#     if i not in infos['incl']:
#         print(i)
print('''}
// clang-format off
''')

