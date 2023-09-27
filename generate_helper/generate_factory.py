#!/bin/python 

import json
import glob
import re
import sys 
import os 

root=sys.argv[1]
config_file = sys.argv[2]
target = sys.argv[3]
factoryfile = sys.argv[4]

if os.path.exists(root+'/config.json'):
    config_file = root + '/config.json'

# print(root, config_file, target, factoryfile)

with open(config_file, 'r', encoding='utf-8') as load_f:
    global data
    data = json.load(load_f)

def creat_object_codes(header, name1, name2, name3):
    _txt = ''.join(open(header.replace('..', root)).readlines())
    codes = ''
    for i in re.findall('class (.*?) : public', _txt, re.S):
        codes += '''
    } else if (%s == %s::name()) {
        %s = new %s(%s);'''%(name1, i, name2, i, name3)
    return codes

objs = []
# print modelfactory.cpp
if target == 'models':
    for i in data['models']:
        list1 = glob.glob(root+'/'+i)
        for j in list1:
            objs += [ j.replace(root, "..") ]

    txt=''
    txt+='#include "modelfactory.h"\n\n'
    txt+='#include <string>\n\n'
    for i in objs:
        txt += '#include "%s.h"\n'%i[:-4]

    txt += '''
Model* init_model(const std::string& model_name, const Param& parm) {
    Model* mymodel = NULL;
    if (false) { // just do nothing (placeholder)'''

    for i in objs:
        txt += creat_object_codes(i[:-4]+'.h', 'model_name', 'mymodel', 'parm')

    txt += '''
    } else {
        LOG(FATAL) << "Cannot parse <forcefield> " << model_name << std::endl;
    }
    return mymodel;
}
'''
    f = open(factoryfile, 'w')
    f.write(txt)
    f.close()

    for i in objs:
        print(i.replace('../models/', ''))

elif target == 'solvers':
    for i in data['solvers']:
        list1 = glob.glob(root+'/'+i)
        for j in list1:
            objs += [ j.replace(root, "..") ]

    txt=''
    txt+='#include "solverfactory.h"\n\n'
    txt+='#include <string>\n\n'
    for i in objs:
        txt += '#include "%s.h"\n'%i[:-4]

    txt += '''
Solver* init_solver(const std::string& solver_name, const Param& parm, Model* pM) {
    Solver* mysolver = NULL;
    if (false) { // just do nothing (placeholder)'''

    for i in objs:
        txt += creat_object_codes(i[:-4]+'.h', 'solver_name', 'mysolver', 'parm, pM')
    txt += '''
    } else {
        LOG(FATAL) << "Cannot parse method name: " << solver_name << std::endl;
    }
    return mysolver;
}
'''
    f = open(factoryfile, 'w')
    f.write(txt)
    f.close()

    for i in objs:
        print(i.replace('../solvers/', '') + ';')


