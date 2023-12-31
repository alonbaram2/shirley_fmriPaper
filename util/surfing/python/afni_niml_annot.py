# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the PyMVPA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
'''
Experimental support for AFNI NIML annotation files 

Created on Feb 19, 2012

@author: Nikolaas. N. Oosterhof (nikolaas.oosterhof@unitn.it)
'''

import numpy as np

import afni_niml as niml
import afni_niml_dset as dset

def rawniml2annot(p):
    if type(p) is list:
        return map(rawniml2annot, p)

    r = dset.rawniml2dset(p)
    for node in p['nodes']:
        name = node.get('name', None)
        if name == 'AFNI_labeltable':
            t = _rawniml_labeltable2annot(node)
            r.update(t)
        elif name == 'AFNI_atr' and node['atr_name'].startswith('UNIQUE_VALS_'):
            data = node['data']
            unique_keys = set([data[i, 0] for i in xrange(len(data))])

    keys = r['keys']
    r['key2row'] = dict((keys[r], r) for r in xrange(len(keys)))

    if len(unique_keys - set(r['keys'])) > 0:
        raise ValueError('key mismatch')

    # clean up
    r.pop('stats')

    return r

def _rawniml_labeltable2annot(p):
    t = dset.rawniml2dset(p)

    r = dict()
    table = t['data']
    nrows = len(table[0])
    r['rgba'] = [(table[0][i],
                table[1][i],
                table[2][i],
                table[3][i]) for i in xrange(nrows)]
    r['keys'] = table[4]

    r['names'] = table[5]

    return r


def read(fn, itemifsingletonlist=True):
    return niml.read(fn, itemifsingletonlist, rawniml2annot)

