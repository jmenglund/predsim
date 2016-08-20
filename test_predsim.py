#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import pytest
import shutil
import subprocess
import tempfile

import dendropy

from os import devnull

from predsim import (
    _read_parameter_file,
    kappa_to_titv,
    get_seqgen_params,
    simulate_matrix,
    simulate_multiple_matrices,
    parse_args,
    main,)


SEQGEN_PATH = 'seq-gen'


def seqgen_status(path):
    """
    Return True if Seq-Gen executable is working,
    otherwise return False.
    """
    f = open(devnull, 'w')
    try:
        subprocess.check_call(
            SEQGEN_PATH, stdout=f, stderr=subprocess.STDOUT)
        status = True
    except subprocess.CalledProcessError:
        status = False
    except OSError:
        status = False
    finally:
        f.close()
    return status


seqgen_required = pytest.mark.skipif(
    seqgen_status(SEQGEN_PATH) is False, reason='Seq-Gen is required')


class TestReadParameterFile():

    p_file_string = (
        '[ID: 1234567890]\n'
        'Gen\tLnL\tLnPr\tTL\tkappa\tpi(A)\tpi(C)\tpi(G)\tpi(T)\talpha\n'
        '0\t-5\t50\t0.5\t1\t0.25\t0.25\t0.25\t0.25\t1\n'
        '500\t-5\t50\t0.5\t1\t0.25\t0.25\t0.25\t0.25\t1')

    def test_read_parameter_file(self, tmpdir):
        f = tmpdir.join('p-file.txt')
        f.write(self.p_file_string)
        p_dicts = _read_parameter_file(str(f.dirpath('p-file.txt')))
        assert len(p_dicts) == 2

    def test_read_empty_parameter_file(self, tmpdir):
        f = tmpdir.join('empty-p-file.txt')
        f.write('')
        with pytest.raises(ValueError):
            _read_parameter_file(str(f.dirpath('empty-p-file.txt')))


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

    rates = '1,1,1,1,1,1'

    def test_simulate_empty_params(self):
        matrix, command = simulate_matrix(SEQGEN_PATH, self.tree)
        assert len(matrix) == 4
        assert matrix.sequence_size == 1000
        assert len(command.split('\t')) == 2
        assert 'HKY' in command

    def test_gtr(self):
        matrix, command = simulate_matrix(
            SEQGEN_PATH, self.tree, general_rates=self.rates)
        assert 'GTR' in command

    def test_ti_tv(self):
        matrix, command = simulate_matrix(SEQGEN_PATH, self.tree, ti_tv=1)
        assert 'HKY' in command
        assert ' -t1 ' in command

    def test_ti_tv_and_gtr(self):
        with pytest.raises(ValueError):
            simulate_matrix(
                SEQGEN_PATH, self.tree, ti_tv=1, general_rates=self.rates)

    def test_gamma(self):
        matrix, command = simulate_matrix(
            SEQGEN_PATH, self.tree, gamma_shape=2)
        assert ' -a2 ' in command


@seqgen_required
class TestMultipleSimulations():

    treelist_string = '((t1:0,t2:0):0,t3:0,t4:0);((t1:0,t2:0):0,t3:0,t4:0);'
    treelist = dendropy.TreeList.get_from_string(treelist_string, 'newick')

    p_dicts = [
        {
            'pi(A)': '0.25',
            'pi(C)': '0.25',
            'pi(G)': '0.25',
            'pi(T)': '0.25',
        }, {
            'pi(A)': '0.25',
            'pi(C)': '0.25',
            'pi(G)': '0.25',
            'pi(T)': '0.25',
        }]

    rng_seeds = ['123321', '456654']

    def test_multiple_simulations(self):
        simulate_multiple_matrices(
            SEQGEN_PATH, self.treelist, self.p_dicts, self.rng_seeds)

    def test_empty_input(self):
        with pytest.raises(AssertionError):
            simulate_multiple_matrices(SEQGEN_PATH, dendropy.TreeList(), [])

    def test_parameter_mismatch(self):
        with pytest.raises(AssertionError):
            simulate_multiple_matrices(
                SEQGEN_PATH, self.treelist[:2], self.p_dicts[:1])

    def test_rng_seeds_mismatch(self):
        with pytest.raises(AssertionError):
            simulate_multiple_matrices(
                SEQGEN_PATH, self.treelist, self.p_dicts, rng_seeds=['123321'])


@seqgen_required
class TestArgumentParser():

    def test_parser_help(self):
        with pytest.raises(SystemExit):
            parse_args(['-h'])

    def test_parser(self):
        with tempfile.NamedTemporaryFile() as p_file:
            with tempfile.NamedTemporaryFile() as t_file:
                with tempfile.NamedTemporaryFile() as commands_file:
                    parse_args([
                        '-l100', '-s1', '-g4',
                        '--commands-file', commands_file.name,
                        p_file.name, t_file.name])


@seqgen_required
class TestMain():

    def test_args_help(self):
        with pytest.raises(SystemExit):
            main(['-h'])

    def test_noargs(self):
        with pytest.raises(SystemExit):
            main()
