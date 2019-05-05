#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Command-line tool for simulating predictive datasets from MrBayes' output.
"""

import argparse
import csv
import os
import shutil
import sys

from collections import namedtuple
from contextlib import ExitStack
from itertools import islice
from math import fabs

import dendropy


__author__ = 'Markus Englund'
__license__ = 'MIT'
__version__ = '0.6.0'


def main(args=None):
    if args is None:
        args = sys.argv[1:]
    parser = parse_args(args)
    tree_list = read_tfile(parser.tfile_path, parser.skip, parser.num_records)
    p_dicts = read_pfile(parser.pfile_path, parser.skip, parser.num_records)
    if parser.seeds_filepath:
        with open(parser.seeds_filepath, 'r') as seeds_fo:
            lines = seeds_fo.readlines()
        rng_seeds = [line for line in lines if line.strip() != '']
    else:
        rng_seeds = None

    simulation_input = combine_simulation_input(tree_list, p_dicts, rng_seeds)

    result_iterator = iter_seqgen_results(
        simulation_input, seq_len=parser.length, gamma_cats=parser.gamma_cats,
        seqgen_path=parser.sg_filepath)

    if parser.out_format == 'nexus':
        schema_kwargs = {'schema': 'nexus', 'simple': False}
    elif parser.out_format == 'phylip':
        schema_kwargs = {'schema': 'phylip'}

    with ExitStack() as cm:  # write to multiple files simultaneously
        write_funcs = []

        if parser.commands_filepath is not None:
            commands_fo = cm.enter_context(open(parser.commands_filepath, 'w'))
            write_funcs.append(get_write_func(commands_fo, 'command'))

        if parser.trees_filepath is not None:
            trees_fo = cm.enter_context(open(parser.trees_filepath, 'w'))
            write_funcs.append(get_write_func(trees_fo, 'tree'))

        for result in result_iterator:
            sys.stdout.write(result.char_matrix.as_string(**schema_kwargs))
            for write_func in write_funcs:
                write_func(result)


def parse_args(args):
    parser = argparse.ArgumentParser(
        prog='predsim', description=(
            'A command-line utility that reads posterior output of MrBayes '
            'and simulates predictive datasets with Seq-Gen.'))
    parser.add_argument(
        '-V', '--version', action='version', version='%(prog)s ' + __version__)
    parser.add_argument(
        '-l', '--length', action='store', default=1000, type=int,
        help='sequence lenght (default: 1000)', metavar='N', dest='length')
    parser.add_argument(
        '-g', '--gamma-cats', action='store', type=int,
        help='number of gamma rate categories (default: continuous)',
        metavar='N', dest='gamma_cats')
    parser.add_argument(
        '-s', '--skip', action='store', default=0, type=int, help=(
            'number of records (trees) to skip at the beginning '
            'of the sample (default: 0)'), metavar='N', dest='skip')
    parser.add_argument(
        '-n', '--num-records', action='store', default=None, type=int,
        help='number of records (trees) to use in the simulation',
        metavar='N', dest='num_records')
    parser.add_argument(
        '-f', '--format', default='nexus', choices=['nexus', 'phylip'],
        help='output format (default: "nexus")', dest='out_format')
    parser.add_argument(
        '-p', '--seqgen-path', default='seq-gen', type=str,
        help='path to a Seq-Gen executable (default: "seq-gen")',
        metavar='FILE', dest='sg_filepath')
    parser.add_argument(
        '--seeds-file', action=StoreExpandedPath, type=is_file,
        help='path to file with seed numbers (e.g. for debugging purposes)',
        metavar='FILE', dest='seeds_filepath')
    parser.add_argument(
        '--commands-file', action=StoreExpandedPath, type=str,
        help='path to output file with commands used by Seq-Gen',
        metavar='FILE', dest='commands_filepath')
    parser.add_argument(
        '--trees-file', action=StoreExpandedPath, type=str,
        help='path to output file with trees used by Seq-Gen',
        metavar='FILE', dest='trees_filepath')
    parser.add_argument(
        'pfile_path', action=StoreExpandedPath, type=is_file,
        help='path to a MrBayes p-file', metavar='pfile')
    parser.add_argument(
        'tfile_path', action=StoreExpandedPath, type=is_file,
        help='path to a MrBayes t-file', metavar='tfile', )

    return parser.parse_args(args)


def read_tfile(filepath, skip=0, num_records=None):
    """
    Read MrBayes t-file into a dendropy TreeList.

    Parameters
    ----------
    filepath : str
    skip : int
        Number of records to skip in the beginning of the file.
    num_records : int
        Number of records to read after the skipped records.

    Returns
    -------
    tree_list : dendropy.TreeList
    """
    tree_list = dendropy.TreeList.get_from_path(
        filepath, 'nexus', tree_offset=skip)[:num_records]
    return tree_list


def read_pfile(filepath, skip=0, num_records=None):
    """
    Read MrBayes p-file into a list of dicts.

    Parameters
    ----------
    filepath : str
    skip : int
        Number of records to skip in the beginning of the file.
    num_records : int
        Number of records to read after the skipped records.

    Returns
    -------
    p_dicts : list
    """
    def process_file(fo, skip=0, stop=None):
        fo.seek(0)
        try:
            next(fo)
        except StopIteration:
            raise ValueError('No records to process in p-file.')
        reader = csv.DictReader(fo, delimiter='\t')
        sliced = islice(reader, skip, stop)
        p_dicts = list(sliced)
        return p_dicts

    stop = skip + num_records if num_records else None
    with open(filepath) as fo:
        p_dicts = process_file(fo, skip=skip, stop=stop)
    return p_dicts


class StoreExpandedPath(argparse.Action):
    """Invoke shell-like path expansion for user- and relative paths."""

    def __call__(self, parser, namespace, values, option_string=None):
        if values:
            filepath = os.path.abspath(os.path.expanduser(str(values)))
            setattr(namespace, self.dest, filepath)


def is_file(filename):
    """Check if a path is a file."""
    if not os.path.isfile(filename):
        msg = '{0} is not a file'.format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename


def kappa_to_titv(kappa, piA, piC, piG, piT):
    """Calculate transistion/transversion ratio from kappa."""
    tot = piA + piC + piG + piT
    if fabs(tot - 1.0) > 1.e-6:
        piA = piA / tot
        piC = piC / tot
        piG = piG / tot
        piT = piT / tot
    titv = kappa * (piA * piG + piC * piT) / ((piA + piG) * (piC + piT))
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
    except KeyError:
        seqgen_params['state_freqs'] = '0.25,0.25,0.25,0.25'
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


def combine_simulation_input(tree_list, p_dicts, rng_seeds=None):
    """Combine input for multiple simulations."""
    assert len(p_dicts) == len(tree_list), (
        'Number of trees does not match the number of records '
        'with parameter values.')
    if rng_seeds is not None:
        assert len(p_dicts) <= len(rng_seeds), (
            'There must be at least ' + str(len(p_dicts)) + ' seed numbers.')
    assert len(p_dicts) > 0, 'No records to process!'
    rng_seeds = rng_seeds if rng_seeds else [None] * len(p_dicts)
    zipped = zip(tree_list, p_dicts, rng_seeds)
    return zipped


def iter_seqgen_results(
        simulation_input, seq_len=1000, gamma_cats=None,
        seqgen_path='seq-gen'):
    """Iterate over multiple simulations."""
    for tree, p_dict, rng_seed in simulation_input:
        seqgen_params = get_seqgen_params(p_dict)
        result = simulate_matrix(
            tree, seq_len=seq_len, rng_seed=rng_seed,
            seqgen_path=seqgen_path, **seqgen_params)
        yield result


def simulate_matrix(
        tree, seq_len=1000, state_freqs=None, ti_tv=None, general_rates=None,
        gamma_shape=None, gamma_cats=None, prop_invar=None, rng_seed=None,
        seqgen_path='seq-gen'):
    """
    Simulate a dataset with Seq-Gen.

    Parameters
    ----------
    tree : dendropy.Tree
    seq_len : int (default: 1000)
        Lengt of sequences to simulate.
    state_freqs : str (default: None)
        State frequences. If `None`, use equal frequences.
    ti_tv : float (default: None)
        transition/transversion rate ratio.
    general_rates : str (default: None)
        General rate matrix.
    gamma_shape : float (default: None)
    gamma_cats : int (default: None)
        Number of gamma rate categories. If `None`,
        use a continous gamma distribution.
    prop_invar : float (default: None)
        Proportion of invariable sites.
    rng_seed : int (default: None)
        Seed for the random number generator. If `None`,
        a seed number will be generated automatically.
    seqgen_path : str (default: "seq-gen")
        Path to Seq-Gen executable.
    """
    if ti_tv and general_rates:
        raise ValueError(
            '"ti_tv" or "general_rates" or both must be set to "None"')
    s = dendropy.interop.seqgen.SeqGen()
    s.seqgen_path = shutil.which(seqgen_path)
    s.char_model = 'GTR' if general_rates else 'HKY'
    s.seq_len = seq_len
    s.state_freqs = state_freqs
    if general_rates:
        s.general_rates = general_rates
    elif ti_tv:
        s.ti_tv = ti_tv
    if gamma_shape:
        s.gamma_shape = gamma_shape
        s.gamma_cats = gamma_cats
    elif gamma_cats:
        raise ValueError(
            'If "gamma_cats" is specified, "gamma_shape" cannot be None')
    s.gamma_shape = gamma_shape
    s.gamma_cats = gamma_cats
    s.prop_invar = prop_invar
    s.rng_seed = rng_seed
    result = SeqGenResult(
        s.generate(tree).char_matrices[0],
        ' '.join(s._compose_arguments()) + '\n',
        tree.as_string(schema='newick'))
    return result


SeqGenResult = namedtuple(
    'SeqGenResult', ['char_matrix', 'command', 'tree'])


def get_write_func(fo, field):
    """
    Return a function for writing data to a file object.

    Parameters
    ----------
    fo : file object
    field : named tuple field name
        Field name in SeqGenResult.
    """
    def write_to_file(result):
        fo.write(getattr(result, field))
        fo.flush()
    return write_to_file


if __name__ == '__main__':  # pragma: no cover
    main()
