#!/usr/bin/env python3

# standard library
import sys
import os
# external modules
import numpy as np

# hijack $PYTHONPATH (sys.path) to find cobramm local modules
sys.path.insert(0, os.path.join(os.getenv('COBRAM_PATH'),'cobramm'))
# now we can import output module
try:
    from output import Output
except ImportError:
    print( 'Error: cannot import cobramm output.py module.\n'
           'Please check if $COBRAM_PATH is defined in your environment, if not run:\n'
           'export COBRAM_PATH=/path/to/COBRAMM' )
    sys.exit(1)


filename = 'cobramm.xml'
# section 1: initial data
# read cobramm output
out = Output(parse=True)

# set default print config for numpy
np.set_printoptions(threshold=sys.maxsize,      # try to avoid truncation for long arrays
                    suppress=True,              # try to avoid scientific notation
                    sign=' ',                   # whitespace in front of positive values
                    floatmode='fixed')          # print all decimals for selected precision (default 8)


# geometry [6 decimals]
geometry = ' ' + np.array2string(np.reshape(out.geometry.cartesian,
                                            (out.geometry.atomNum * 3, ),
                                            order='F'),
                                 precision=6,
                                 separator='\t',
                                 max_line_width=40,
                                 )[1:-1]

# atom links BA and DC are merged in a single 2d array
_atom_links = np.concatenate((out.geometry.atomLink_BA,
                              out.geometry.atomLink_DC), axis=1)
_atom_links = _atom_links.astype('int64')
if len(_atom_links[0]) > 0:
    atom_links = '{0}\n{1}'.format('\t'.join([str(s) for s in _atom_links[1]]),
                                   '\t'.join([str(s) for s in _atom_links[0]]))
else:
    atom_links = None

# formatter handle spacing between mono and double chars labels and print unquoted strings
atom_labels = ' ' + np.array2string(out.geometry.atomLabel,
                                    separator='',
                                    max_line_width=50,
                                    formatter={'str_kind': lambda x: '{0}{1}'.format(x, '   ' if len(x)==1 else '  ')}
                                    )[1:-1]

atom_layers = ' ' + np.array2string(out.geometry.atomType,
                                    separator='  ',
                                    max_line_width=50,
                                    formatter={'str_kind': lambda x: x}
                                    )[1:-1]

# charges [8 decimals]
if np.any(out.chargeMM):
    chargeMM = ' ' + np.array2string(out.chargeMM,
                                     separator='\t',
                                     )[1:-1]
else:
    chargeMM = None

if np.any(out.chargeQM):
    chargeQM = ' ' + np.array2string(out.chargeQM,
                                     separator='\t',
                                     )[1:-1]
else:
    chargeQM = None

# list of unchanged attributes
i_list = ['version', 'calculation_type', 'optimization_type', 'qm_type',
          'modelH_top','real_top']
# list of crafted strings
i_mod_list = ['geometry', 'atom_labels', 'atom_layers', 'atom_links',
              'chargeMM', 'chargeQM']
# init empty dict
i_keys = {}
# populate dict
for k in i_list:
    exec( "i_keys['{0}'] = out.{0}".format(k) )
for k in i_mod_list:
    exec( "i_keys['{0}'] = {0}".format(k) )

# prepare text
text = '<INITIAL DATA>\n'
# loop over lists to preserve attributes ordering (dicts are NOT ordered in python < 3.7)
for k in i_list + i_mod_list:
    if i_keys[k] is not None:
        text += '    <{0}>\n{1}\n    </{0}>\n'.format(k, i_keys[k])
text += '</INITIAL DATA>\n'

# write initial data to file
with open('HR_{0}'.format(filename), 'w') as out_file:
    out_file.write(text)


# only for frequencies calculations, write FREQ section with high precision normal modes
if out.optimization_type in ['freqxg', 'freqxgp']:
    # prepare frequencies and normal modes text
    text = '<FREQ>\n'
    text += out.freq_to_string()
    text += '\n</FREQ>\n'

    with open('HR_{0}'.format(filename), 'a') as out_file:
        out_file.write(text)


# section 2: steps
# loop over out.steps + 1 to include steps "0 to last"
for n in range(out.steps + 1):
    # fetch Step object
    step = out.get_step(n)

    # geometry [6 decimals]
    geometry = ' ' + np.array2string(np.reshape(step.geometry.cartesian,
                                                (step.geometry.atomNum * 3, ),
                                                order='F'),
                                     precision=6,
                                     separator='\t',
                                     max_line_width=40,
                                     )[1:-1]

    # charges [8 decimals]
    if np.any(step.chargeMM):
        chargeMM = ' ' + np.array2string(step.chargeMM,
                                         separator='\t',
                                         )[1:-1]
    else:
        chargeMM = None

    if np.any(step.chargeQM):
        chargeQM = ' ' + np.array2string(step.chargeQM,
                                         separator='\t',
                                         )[1:-1]
    else:
        chargeQM = None

    # forces [8 decimals]
    if np.any(step.forces_QMMM_HM):
        forces_QMMM_HM = ' ' + np.array2string(np.reshape(step.forces_QMMM_HM,
                                                          (np.prod(np.shape(step.forces_QMMM_HM)), ),
                                                          order='F'),
                                               separator='\t',
                                               max_line_width=40,
                                               )[1:-1]

    else:
        forces_QMMM_HM = None

    # energy [8 decimals]
    if np.any(step.E_QM):
        E_QM = np.array2string(np.asarray(step.E_QM),
                               separator='\n',
                               )[1:-1]
    else:
        E_QM = None

    if np.any(step.E_QMMM):
        E_QMMM = np.array2string(np.asarray(step.E_QMMM),
                                 separator='\n',
                                 )[1:-1]
    else:
        E_QMMM = None

    # velocity [10 decimals]
    if np.any(step.velocity):
        velocity = ' ' + np.array2string(np.reshape(step.velocity,
                                                    (np.prod(np.shape(step.velocity)), ),
                                                    order='F'),
                                         precision=10,
                                         separator='\t',
                                         max_line_width=50,
                                         )[1:-1]

    else:
        velocity = None

    # all attributes list
    key_list = ['step', 'geometry', 'chargeMM', 'chargeQM', 'forces_QMMM_HM',
                'E_MM_real', 'E_MM_nocharge', 'E_MM_modelH', 'E_self_emb', 'E_QM', 'E_QMMM',
                'time', 'state', 'velocity',
                'MD_Etot', 'MD_Epot', 'MD_Ekin', 'MD_temp',
                'optstate', 'nIRCstep', 'Fmax', 'Frms', 'Dmax', 'Drms']

    # unchanged attributes
    s_unchanged_list = ['step', 'time', 'state', 'nIRCstep', 'optstate']
    # float attributes (6 decimals)
    s_float_list = ['E_MM_real', 'E_MM_nocharge', 'E_MM_modelH', 'E_self_emb',
                    'MD_Etot', 'MD_Epot', 'MD_Ekin', 'MD_temp',
                    'Fmax', 'Frms', 'Dmax', 'Drms']
    # hand crafted attributes
    s_mod_list = ['geometry', 'chargeMM', 'chargeQM',
                  'forces_QMMM_HM', 'E_QM', 'E_QMMM', 'velocity']

    step_keys = {}

    # prepare all strings
    for k in s_unchanged_list:
        exec( "step_keys['{0}'] = step.{0} if step.{0} is not None else None".format(k) )
    for k in s_float_list:
        result = step.get_attribute(k)
        if result:
            step_keys[k] = '{0:.8f}'.format(result)
        else:
            step_keys[k] = None
    for k in s_mod_list:
        exec( "step_keys['{0}'] = {0}".format(k) )

    # prepare step text
    text = '<STEP {0}>\n'.format(n)
    for k in key_list:
        if step_keys[k] is not None:
            text += '    <{0}>\n{1}\n    </{0}>\n'.format(k, step_keys[k])
    text += '</STEP {0}>\n'.format(n)

    # write text to file
    with open('HR_{0}'.format(filename), 'a') as out_file:
        out_file.write(text)
