'''FITS file and data reader'''

import mmap
import math
import re
from weakref import finalize

import numpy as np

from .header import HDUError, read_header


class HduArray(np.ndarray):

    def __new__(cls, data, header=None):
        obj = np.asarray(data).view(cls)
        obj.header = header
        return obj

    def __array_finalize__(self, obj):
        if obj is not None:
            self.header = getattr(obj, 'header', None)


def shape_from_header(header):
    naxis = header['NAXIS']
    if naxis > 0:
        shape = tuple(header[f'NAXIS{i+1}'] for i in range(naxis))
    else:
        shape = None
    return shape


tform_regex = re.compile(r'''
    (?P<count>[0-9]+)?          # repeat count
    (?P<value>[LXBIJKAEDCM])    # data type, excluding variable array
    .*                          # additional characters (unused)
''', re.X)


BITPIX_DTYPE_MAP = {
    8: 'u1',
    16: 'i2',
    32: 'i4',
    64: 'i8',
    -32: 'f4',
    -64: 'f8',
}


TFORM_DTYPE_MAP = {
    'L': '?',
    'X': 'u1',
    'B': 'u1',
    'I': 'i2',
    'J': 'i4',
    'K': 'i8',
    'A': 'U',
    'E': 'f4',
    'D': 'f8',
    'C': 'c8',
    'M': 'c16'
}


def image_dtype(header):
    bitpix = header['BITPIX']
    try:
        dtype = BITPIX_DTYPE_MAP[bitpix]
    except KeyError:
        raise ValueError(f'invalid BITPIX value: {bitpix}') from None
    return f'>{dtype}'


def bintable_dtype(header):
    dtype = []
    tfields = header['TFIELDS']
    for n in range(tfields):
        ttype = header.get(f'TTYPE{n+1}', f'f{n}')
        tform = header[f'TFORM{n+1}']

        name = ttype.rstrip()

        if (m := tform_regex.fullmatch(tform)) is None:
            raise ValueError(f'TFORM{n+1}: invalid format: {tform}')
        count, value = m.group('count'), m.group('value')

        fmt = TFORM_DTYPE_MAP[value]

        if count is not None:
            if value == 'X':
                c = int(count)
                if c % 8:
                    c += 1
                fmt = f'{c}{fmt}'
            elif value == 'A':
                fmt = f'{fmt}{count}'
            else:
                fmt = f'{count}{fmt}'

        dtype += [(name, f'>{fmt}')]

    return dtype


def dtype_from_header(header):
    kind = header.get('XTENSION', 'IMAGE').rstrip()

    if kind == 'IMAGE':
        dtype = image_dtype(header)
    elif kind == 'BINTABLE':
        dtype = bintable_dtype(header)
    else:
        raise ValueError(f'{kind} extension not supported')

    return dtype


def skip_data(fp, header):
    shape = shape_from_header(header)
    if shape is not None:
        bitpix = header['BITPIX']
        pcount = header.get('PCOUNT', 0)
        ndata = abs(bitpix)//8 * (pcount + math.prod(shape))
        nskip = ndata + (2880 - ndata % 2880) % 2880
        fp.seek(nskip, 1)


def hdu(fp, header):
    shape = shape_from_header(header)
    dtype = dtype_from_header(header)
    buffer = data = None

    if shape is not None:

        bitpix = header['BITPIX']
        pcount = header.get('PCOUNT', 0)
        gcount = header.get('GCOUNT', 1)

        nbytes = abs(bitpix)//8 * gcount * (pcount + math.prod(shape))

        if isinstance(dtype, list):
            shape = shape[1:]

        file_offset = fp.tell()

        array_offset = file_offset % mmap.ALLOCATIONGRANULARITY

        buffer_offset = file_offset - array_offset
        buffer_length = array_offset + nbytes

        try:

            buffer = mmap.mmap(fp.fileno(), buffer_length,
                               offset=buffer_offset, access=mmap.ACCESS_READ)

            data = np.ndarray.__new__(np.ndarray, shape, dtype,
                                      buffer=buffer, offset=array_offset,
                                      order='F')

        finally:

            if buffer is not None and data is None:
                buffer.close()

        finalize(data, buffer.close)

    else:

        data = np.empty(0, dtype)

    return HduArray(data, header)


def load(fp):
    hdus = []

    buffer = bytearray(2880)

    try:
        header = read_header(fp, buffer)
    except HDUError:
        raise TypeError('not a FITS file')

    hdus.append(hdu(fp, header))

    if header['EXTEND'] is True:

        while True:

            skip_data(fp, header)

            try:
                header = read_header(fp, buffer, False)
            except HDUError:
                break

            hdus.append(hdu(fp, header))

    return hdus
