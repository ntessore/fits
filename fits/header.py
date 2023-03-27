'''FITS header reader'''

import re
from enum import Enum


class HDUError(TypeError):
    pass


class KeywordValue(str, Enum):
    '''Grammar for FITS keyword values.

    Based on Appendix A of the FITS Standard Document v4.0.

    '''

    __str__ = str.__str__

    STR = '(?: \' [ -~]* \' )+'

    BOOL = '[TF]'

    INT = '[-+]? [0-9]+'

    FLOAT = r'''
        [-+]? (?: [0-9]+ \.? | [0-9]* \. [0-9]+ )   # decimal
        (?: [ED] [-+]? [0-9]+ )?                    # exponent
    '''

    COMPLEX = fr'''
        \( \ * (?P<real> {FLOAT} ) \ * ,   # real part
        \ * (?P<imag> {FLOAT} ) \ * \)     # imaginary part
    '''

    VALUE = f'''
        (?P<value>
            (?P<string> {STR}) |
            (?P<bool> {BOOL}) |
            (?P<int> {INT}) |
            (?P<float> {FLOAT}) |
            (?P<complex> {COMPLEX})
        )
    '''

    COMMENT = '(?: / (?P<comment> [ -~]* ) )'

    PARSE = fr'\ * {VALUE}? \ * {COMMENT}?'


keyword_value_regex = re.compile(KeywordValue.PARSE, re.X)


def parse_keyword_value(expr):
    m = keyword_value_regex.fullmatch(expr)

    if m is None:
        raise ValueError('invalid keyword value')
    elif m.group('value') is None:
        value = None
    elif (s := m.group('string')) is not None:
        value = s[1:-1].replace("''", "'")
    elif (b := m.group('bool')) is not None:
        value = (b == 'T')
    elif (i := m.group('int')) is not None:
        value = int(i)
    elif (f := m.group('float')) is not None:
        value = float(f.replace('D', 'E'))
    elif m.group('complex') is not None:
        real, imag = m.group('real'), m.group('imag')
        value = complex(float(real), float(imag))
    else:
        # should never happen
        value = None

    if (co := m.group('comment')) is not None:
        comment = co.strip()
    else:
        comment = None

    return value, comment


def parse_keyword(record):
    keyword = record[:8].rstrip(' ')

    if keyword == 'END':
        keyword, value, comment = None, None, None

    elif (keyword == 'CONTINUE'
          or keyword not in ('COMMENT', 'HISTORY', '')
          and record[8:10] == '= '):

        value, comment = parse_keyword_value(record[10:])

    else:

        value, comment = record[8:].rstrip(' '), None

    return keyword, value, comment


def read_block(fp, buffer=None):
    if buffer is None or len(buffer) != 2880:
        buffer = bytearray(2880)
    size = fp.readinto(buffer)
    if size != 2880:
        buffer = None
    return buffer


def read_header(fp, buffer=None, primary=True):
    keywords = {}
    nrecords = 0

    ended = False
    while not ended:

        block = read_block(fp, buffer)

        if block is None:
            if nrecords == 0:
                raise HDUError
            else:
                raise TypeError('incomplete block in header')

        if min(block) < 32 or max(block) > 126:
            raise TypeError('invalid characters in header')

        for i in range(0, 2880, 80):
            record = block[i:i+80].decode('ascii')
            nrecords += 1

            try:
                keyword, value, comment = parse_keyword(record)
            except ValueError as exc:
                if nrecords == 1:
                    raise HDUError
                if (m := re.fullmatch(r'(\w+)\s*', record[:8])):
                    msg = f' [{m.group(1)}?]'
                else:
                    msg = ''
                raise ValueError(f'record {nrecords}{msg}: {exc}') from None

            if nrecords == 1:
                if primary:
                    if keyword != 'SIMPLE' or value is not True:
                        raise HDUError
                else:
                    if keyword != 'XTENSION' or not isinstance(value, str):
                        raise HDUError

            if keyword is None:
                ended = True
                break

            # TODO: keep comment
            if value is not None:
                ndup = 0
                name = keyword
                while name in keywords:
                    ndup += 1
                    name = f'{keyword}.{ndup}'
                keywords[name] = value

    return keywords
