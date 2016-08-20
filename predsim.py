#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Command-line tool for simulating predictive datasets from MrBayes' output.
"""

import argparse
import csv
import itertools
import math
import shutil
import sys

import dendropy


__authors__ = 'Markus Englund'
__license__ = 'MIT'
__version__ = '0.2.0'


def _read_parameter_file(filepath, skip=0, num_records=None):
    """Read MrBayes p-file into a list of dicts.

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
    def process_file(p_file, skip=0, num_records=None):
        p_file.seek(0)
        try:
            next(p_file)
        except StopIteration:
            raise ValueError('No records to process in p-file.')
        reader = csv.DictReader(p_file, delimiter='\t')
        sliced = itertools.islice(reader, skip, stop)
        p_dicts = list(sliced)
        return p_dicts

    stop = skip + num_records if num_records else None
    if (sys.version_info >= (3, 0)):
        with open(filepath, newline='') as p_file:
            p_dicts = process_file(p_file, skip=skip, num_records=num_records)
    else:
        with open(filepath) as p_file:
            p_dicts = process_file(p_file, skip=skip, num_records=num_records)
    return p_dicts


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
        seqgen_path, tree, seq_len=1000, state_freqs=None,
        ti_tv=None, general_rates=None, gamma_shape=None,
        gamma_cats=None, prop_invar=None, rng_seed=None):
    """
    Simulate a dataset with Seq-Gen.

    Parameters
    ----------
    seqgen_path : str
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
    """
    if ti_tv and general_rates:
        raise ValueError(
            '"ti_tv" or "general_rates" or both must be set to "None"')
    s = dendropy.interop.seqgen.SeqGen()
    if (sys.version_info >= (3, 3)):
        full_path = shutil.which(seqgen_path)
        s.seqgen_path = full_path if full_path else seqgen_path
    else:
        s.seqgen_path = seqgen_path
    s.char_model = 'GTR' if general_rates else 'HKY'
    s.seq_len = seq_len
    s.state_freqs = state_freqs
    if general_rates:
        s.general_rates = general_rates
    elif ti_tv:
        s.ti_tv = ti_tv
    s.gamma_shape = gamma_shape
    s.gamma_cats = gamma_cats
    s.prop_invar = prop_invar
    s.rng_seed = rng_seed
    matrix = s.generate(tree).char_matrices[0]
    tree_string = tree.as_string('newick', suppress_rooting=True).rstrip()
    command = ' '.join(s._compose_arguments()) + '\t' + tree_string
    return (matrix, command)


def simulate_multiple_matrices(
        seqgen_path, tree_list, p_dicts, rng_seeds=None,
        seq_len=1000, gamma_cats=None):
    """
    Simulate multiple predictive datasets with Seq-Gen.

    Parameters
    ----------
    seqgen_path : str
        Path to Seq-Gen executable.
    tree_list : Dendropy.TreeList
    p_dicts : list
        Parameter values from a MrBayes' p-file.
    rng_seeds : list
        Random seed numbers passed to Seq-Gen.
    seq_len : int (default: 1000)
        Lengt of sequences to simulate.
    gamma_cats : int (default: None)
        Number of discrete gamma rate categories
        (continous by default).

    Returns
    -------
    matrices, commands : tuple (dendropy.DataSet, list)
        Simulated dataset and used Seq-Gen commands.
    """
    assert len(p_dicts) == len(tree_list), (
        'Number of trees does not match the number of parameter values.')
    if rng_seeds is not None:
        assert len(p_dicts) == len(rng_seeds), (
            'Number of seed numbers does not match '
            'the number of parameter values.')
    assert len(p_dicts) > 0, 'No records to process!'
    rng_seeds = rng_seeds if rng_seeds else [None] * len(p_dicts)
    zipped = zip(tree_list, p_dicts, rng_seeds)
    matrices = dendropy.DataSet()
    seqgen_commands = []
    for tree, p_dict, rng_seed in zipped:
        seqgen_params = get_seqgen_params(p_dict)
        matrix, command = simulate_matrix(
            seqgen_path, tree, seq_len=seq_len, rng_seed=rng_seed,
            **seqgen_params)
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
        metavar='N', dest='length',
        help='sequence lenght (default: 1000)')
    parser.add_argument(
        '-g', '--gamma-cats', type=int, action='store', metavar='N',
        dest='gamma_cats',
        help='number of gamma rate categories (default: continuous)')
    parser.add_argument(
        '-s', '--skip', type=int, action='store', metavar='N',
        dest='skip', default=0, help=(
            'number of records (trees) to skip at the beginning '
            'of the sample (default: 0)'))
    parser.add_argument(
        '-n', '--num-records', type=int, action='store', metavar='N',
        dest='num_records', default=None, help=(
            'number of records (trees) to use in the simulation'))
    parser.add_argument(
        '-p', '--seqgen-path',
        type=str, default='seq-gen',
        dest='seqgen_path', metavar='FILE',
        help='path to a Seq-Gen executable (default: "seq-gen")')
    parser.add_argument(
        '--random-seeds-file',
        type=argparse.FileType('rU'),
        dest='random_seeds_file', metavar='FILE',
        help='path to file with seed numbers to pass to Seq-Gen')
    parser.add_argument(
        '--commands-file',
        type=argparse.FileType('w'),
        dest='commands_file', metavar='FILE',
        help='path to output file with used Seq-Gen commands')
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
    tree_list = dendropy.TreeList.get_from_stream(
        parser.tfile, 'nexus', tree_offset=parser.skip)[:parser.num_records]
    p_dicts = _read_parameter_file(
        parser.pfile.name, parser.skip, parser.num_records)
    if parser.random_seeds_file:
        rng_seeds = parser.random_seeds_file.read().splitlines()
        rng_seeds = [line for line in rng_seeds if line.strip() != '']
        print(rng_seeds)
    else:
        rng_seeds = None
    matrices, seqgen_commands = simulate_multiple_matrices(
        parser.seqgen_path, tree_list, p_dicts, rng_seeds=rng_seeds,
        seq_len=parser.length, gamma_cats=parser.gamma_cats)
    if parser.commands_file:
        parser.commands_file.write('\n'.join(seqgen_commands))
    parser.outfile.write(matrices.as_string(schema='nexus'))


if __name__ == '__main__':  # pragma: no cover
    main()
