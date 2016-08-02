#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Command-line tool for simulating predictive datasets from MrBayes' output.
"""

import argparse
import math
import sys

import dendropy
import pandas


__authors__ = 'Markus Englund'
__license__ = 'MIT'
__version__ = '0.1.0'


def _get_skiprows(cnt):
    """
    Return a list with rows to skip when
    reading the p-file.

    Parameters
    ----------
    cnt : int
        Number of rows to be skipped.
    """
    skiprows = [0]
    if cnt > 0:
        skiprows.extend(range(2, cnt + 2))
    return skiprows


def _iterrecords(frame):
    """Iterate over the rows of a DataFrame as dicts."""
    for record in frame.to_dict('records'):
        yield record


def kappa_to_titv(kappa, piA, piC, piG, piT):
    """Calculate transistion/transversion ratio from kappa."""
    tot = piA + piC + piG + piT
    if math.fabs(tot - 1.0) > 1.e-6:
        piA = piA/tot
        piC = piC/tot
        piG = piG/tot
        piT = piT/tot
    titv = kappa * (piA * piG + piC * piT)/((piA + piG) * (piC + piT))
    return titv


def get_seqgen_params(mrbayes_params):
    """
    Adapt MrBayes parameter values for use with Seq-Gen.

    Paramters
    ---------
    mrbayes_prams : dict
        Parameter values from a single row in a MrBayes p-file.

    Returns
    -------
    seqgen_params : dict
    """
    seqgen_params = {}
    try:
        seqgen_params['state_freqs'] = (
            str(mrbayes_params['pi(A)']) + ',' +
            str(mrbayes_params['pi(C)']) + ',' +
            str(mrbayes_params['pi(G)']) + ',' +
            str(mrbayes_params['pi(T)']))
    except KeyError as ex:
        raise KeyError(
            'Could not find any base frequences:\n{ex}'.format(ex=str(ex)))
    try:
        seqgen_params['ti_tv'] = kappa_to_titv(
            float(mrbayes_params['kappa']),
            float(mrbayes_params['pi(A)']),
            float(mrbayes_params['pi(C)']),
            float(mrbayes_params['pi(G)']),
            float(mrbayes_params['pi(T)']))
    except KeyError:
        pass
    try:
        seqgen_params['general_rates'] = (
            str(mrbayes_params['r(A<->C)']) + ',' +
            str(mrbayes_params['r(A<->G)']) + ',' +
            str(mrbayes_params['r(A<->T)']) + ',' +
            str(mrbayes_params['r(C<->G)']) + ',' +
            str(mrbayes_params['r(C<->T)']) + ',' +
            str(mrbayes_params['r(G<->T)']))
    except KeyError:
        pass
    try:
        seqgen_params['gamma_shape'] = str(mrbayes_params['alpha'])
    except KeyError:
        pass
    try:
        seqgen_params['prop_invar'] = str(mrbayes_params['pinvar'])
    except KeyError:
        pass
    return seqgen_params


def simulate_matrix(
        seqgen_path, seqgen_params, tree, seq_len=1000, gamma_cats=None):
    """
    Simulate a predictive datasets with Seq-Gen.

    Parameters
    ----------
    seqgen_path : str
        Path to Seq-Gen executable.
    seqgen_params : dict
        Parameter inputs for Seq-Gen.
    tree : Dendropy.Tree
    seq_len : int (default: 1000)
        Lengt of sequences to simulate.
    gamma_cats : int (default: None)
        Number of discrete gamma rate categories
        (continous by default).

    Returns
    -------
    matrix, command : tuple (dendropy.CharacterMatrix, str)
        Simulated dataset and used Seq-Gen command.
    """
    s = dendropy.interop.seqgen.SeqGen()
    s.seqgen_path = seqgen_path
    s.char_model = 'GTR' if ('general_rates' in seqgen_params) else 'HKY'
    s.seq_len = seq_len
    try:
        s.gamma_shape = seqgen_params['gamma_shape']
        s.gamma_cats = gamma_cats
    except KeyError:
        pass
    try:
        s.prop_invar = seqgen_params['prop_invar']
    except KeyError:
        pass
    try:
        s.state_freqs = seqgen_params['state_freqs']
    except KeyError:
        pass
    if 'general_rates' in seqgen_params:
        s.general_rates = seqgen_params['general_rates']
    elif 'ti_tv' in seqgen_params:
        s.ti_tv = seqgen_params['ti_tv']
    matrix = s.generate(tree).char_matrices[0]
    tree_string = tree.as_string('newick', suppress_rooting=True).rstrip()
    command = ' '.join(s._compose_arguments()) + '\t' + tree_string
    return (matrix, command)


def simulate_multiple_matrices(
        seqgen_path, p_frame, tree_list, seq_len=1000, gamma_cats=None):
    """
    Simulate multiple predictive datasets with Seq-Gen.

    Parameters
    ----------
    seqgen_path : str
        Path to Seq-Gen executable.
    p_frame : pandas.DataFrame
        Parameter inputs for Seq-Gen.
    tree_list : Dendropy.TreeList
    seq_len : int (default: 1000)
        Lengt of sequences to simulate.
    gamma_cats : int (default: None)
        Number of discrete gamma rate categories
        (continous by default).

    Returns
    -------
    matrix, command : tuple (dendropy.CharacterMatrix, str)
        Simulated dataset and used Seq-Gen command.
    """
    if len(p_frame) == 0:
        raise ValueError('No parameter values found')
    if len(p_frame) != len(tree_list):
        raise ValueError(
            'Number of parameter values do not match the number of trees.')
    p_frame['tree'] = tree_list
    matrices = dendropy.DataSet()
    seqgen_commands = []
    for record in _iterrecords(p_frame):
        tree = record.pop('tree')
        seqgen_params = get_seqgen_params(record)
        matrix, command = simulate_matrix(
            seqgen_path, seqgen_params, tree, seq_len, gamma_cats)
        matrices.add(matrix)
        seqgen_commands.append(command)
    matrices.unify_taxon_namespaces()
    return (matrices, seqgen_commands)


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=(
            'A command-line utility that reads posterior '
            'output of MrBayes and simulates predictive '
            'datasets with Seq-Gen.'))
    parser.add_argument(
        '-V', '--version', action='version',
        version='predsim ' + __version__)
    parser.add_argument(
        '-l', '--length', type=int, action='store', default=1000,
        metavar='INT', dest='length',
        help='sequence lenght (default: 1000)')
    parser.add_argument(
        '-g', '--gamma-cats', type=int, action='store', metavar='INT',
        dest='gamma_cats',
        help='number of gamma rate categories (default: continuous)')
    parser.add_argument(
        '-c', '--commands-file',
        type=argparse.FileType('w'),
        dest='commands_file', metavar='PATH',
        help='path to output file with used Seq-Gen commands')
    parser.add_argument(
        '-s', '--skip', type=int, action='store', metavar='INT',
        dest='skip', default=0, help=(
            'number of records (trees) to skip at the beginning '
            'of the sample (default: 0)'))
    parser.add_argument(
        '-p', '--seqgen-path',
        type=str, default='seq-gen',
        dest='seqgen_path', metavar='PATH',
        help='path to a Seq-Gen executable (default: "seq-gen")')
    parser.add_argument(
        'pfile', type=argparse.FileType('rU'),
        help='path to a MrBayes p-file')
    parser.add_argument(
        'tfile', type=argparse.FileType('rU'),
        help='path to a MrBayes t-file')
    parser.add_argument(
        'outfile', nargs='?', type=argparse.FileType('w'),
        default=sys.stdout,
        help='path to output file (default: <stdout>)')
    return parser.parse_args(args)


def main(args=None):
    if args is None:
        args = sys.argv[1:]
    parser = parse_args(args)
    skiprows = _get_skiprows(parser.skip)
    p_frame = pandas.read_table(parser.pfile, skiprows=skiprows)
    tree_list = dendropy.TreeList.get_from_stream(
        parser.tfile, 'nexus', tree_offset=parser.skip)
    simulated_matrices, seqgen_commands = simulate_multiple_matrices(
        parser.seqgen_path,
        p_frame,
        tree_list,
        seq_len=parser.length,
        gamma_cats=parser.gamma_cats)
    if parser.commands_file:
        parser.commands_file.write('\n'.join(seqgen_commands))
    parser.outfile.write(simulated_matrices.as_string(schema='nexus'))


if __name__ == '__main__':  # pragma: no cover
    main()
