#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import tempfile
import os.path

import pandas
import dendropy

from pandas.util.testing import assert_dict_equal

from predsim import (
    _get_skiprows,
    _iterrecords,
    kappa_to_titv,
    get_seqgen_params,
    simulate_matrix,
    simulate_multiple_matrices,
    parse_args,
    main,)

SEQGEN_PATH = 'seq-gen'


seqgen_required = pytest.mark.skipif(
    os.path.isfile(SEQGEN_PATH) is False, reason='Seq-Gen is required')


class TestGetSkiprows():

    def test_noskip(self):
        assert _get_skiprows(0) == [0]

    def test_skip(self):
        assert _get_skiprows(2) == [0, 2, 3]


class TestIterRecords():

    frame = pandas.DataFrame({'a': [1, 2, 3], 'b': [4, 5, 6]})

    def test_iterrecords(self):
        assert isinstance(list(_iterrecords(self.frame))[0], dict)


class TestKappaConversion():

    def test_equal_basefreqs(self):
        assert kappa_to_titv(1, .25, .25, .25, .25) == 0.5

    def test_unequal_basefreqs(self):
        assert kappa_to_titv(1, .2, .2, .3, .3) == 0.48

    def test_low_basefreqs_sum(self):
        assert kappa_to_titv(1, .1, .1, .1, .1) == 0.5

    def test_invalid_basefreqs(self):
        with pytest.raises(ZeroDivisionError):
            kappa_to_titv(1, .5, .0, .5, .0)


class TestGetSeqGenParamaters():

    d1 = {'pi(A)': 0.25, 'pi(C)': 0.25,	'pi(G)': 0.25,	'pi(T)': 0.25}
    d2 = {'state_freqs':  '0.25,0.25,0.25,0.25'}
    d3 = {'pi(A)': 0.25, 'pi(C)': 0.25,	'pi(G)': 0.25}

    def test_equal_basefreqs(self):
        assert get_seqgen_params(self.d1) == self.d2

    def test_missing_basefreq(self):
        with pytest.raises(KeyError):
            get_seqgen_params(self.d3)


@seqgen_required
class TestSingleSimulation():

    tree_string = '((t1:0,t2:0):0,t3:0,t4:0);'
    tree = dendropy.Tree.get_from_string(tree_string, 'newick')

    rates = {'general_rates': '1,1,1,1,1,1'}

    def test_simulate_empty_params(self):
        matrix, command = simulate_matrix(SEQGEN_PATH, {}, self.tree)
        assert len(matrix) == 4
        assert matrix.sequence_size == 1000
        assert len(command.split('\t')) == 2
        assert 'HKY' in command

    def test_gtr(self):
        matrix, command = simulate_matrix(SEQGEN_PATH, self.rates, self.tree)
        assert 'GTR' in command

    def test_ti_tv(self):
        matrix, command = simulate_matrix(
            SEQGEN_PATH, {'ti_tv': 1}, self.tree)
        assert 'HKY' in command
        assert ' -t1 ' in command

    def test_gamma(self):
        matrix, command = simulate_matrix(
            SEQGEN_PATH, {'gamma_shape': 2}, self.tree)
        assert ' -a2 ' in command


@seqgen_required
class TestMultipleSimulations():

    treelist_string = '((t1:0,t2:0):0,t3:0,t4:0);((t1:0,t2:0):0,t3:0,t4:0);'
    treelist = dendropy.TreeList.get_from_string(treelist_string, 'newick')

    p_frame = pandas.DataFrame({
        'pi(A)': [0.25, 0.25],
        'pi(C)': [0.25, 0.25],
        'pi(G)': [0.25, 0.25],
        'pi(T)': [0.25, 0.25]})

    def test_empty_input(self):
        with pytest.raises(ValueError):
            simulate_multiple_matrices(
                SEQGEN_PATH, pandas.DataFrame(), dendropy.TreeList())

    def test_two_trees(self):
        simulate_multiple_matrices(
            SEQGEN_PATH, self.p_frame, self.treelist)

    def test_parameter_mismatch(self):
        with pytest.raises(ValueError):
            simulate_multiple_matrices(
                SEQGEN_PATH, self.p_frame[:1], self.treelist[:2])


@seqgen_required
class TestArgumentParser():

    def test_parser_help(self):
        with pytest.raises(SystemExit):
            parse_args(['-h'])

    def test_parser(self):
        with tempfile.NamedTemporaryFile() as p_file:
            with tempfile.NamedTemporaryFile() as t_file:
                with tempfile.NamedTemporaryFile() as command_file:
                    parse_args([
                        '-l100', '-s1', '-g4', '-c', command_file.name,
                        p_file.name, t_file.name])


@seqgen_required
class TestMain():

    def test_args_help(self):
        with pytest.raises(SystemExit):
            main(['-h'])

    def test_noargs(self):
        with pytest.raises(SystemExit):
            main()
