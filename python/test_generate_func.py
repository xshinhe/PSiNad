import re
line = "<double, std::string>, <std::map<int, double>>, <<int>, 2>"


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
        double* abc=0
        double* abc_d
        double* abc_d[]
        double* abc_d[2]
        double*&abc_d[]
        const std::sting &abc_d[]
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

def safe_split_comma(line_in):
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

def parse_type_and_name_multiple(term):

    _tmps = safe_split_comma(re.split(';', term)[0]);

    type0, name0 = parse_type_and_name(_tmps[0])
    arg_type = [type0]
    arg_name = [name0]

    res =  re.search('(\W)*?$', type0, re.S)
    if res:
        i1, i2 = res.span()
        type0 = type0[:i1]

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

def test(term):
    '''
        this function split (unique) typeinfo & variable name.
        passed for:
        double* abc=0
        double* abc_d
        double* abc_d[]
        double* abc_d[2]
        double*&abc_d[]
        const std::sting &abc_d[]
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
    typi = re.sub(' +?', ' ', typi)
    typi = re.sub('(\w) (\w)', r'\1I__0__I\2', typi)
    typi = re.sub(' ', '', typi)
    typi = re.sub('I__0__I', ' ', typi)

    return typi, name


print(safe_split_comma(line))

print(test('double* abc=0'))
print(test('double* abc_d'))
print(test('double* abc_d[]'))
print(test('double* abc_d[2]'))
print(test('double*&abc_d[]'))
print(test('const std::sting &abc_d[]'))

# while res:
#     for i in res:
        # line.replace()
# if res:
# i1, i2 = res.span()
# old = line[i1:i2]
# line = re.sub(old, '<__>', line)

print(parse_type_and_name_multiple(
    "const int *a, b, **c;"
))