#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import subprocess
import tempfile

import pytest
import dendropy

from predsim import (
    read_tfile,
    read_pfile,
    kappa_to_titv,
    get_seqgen_params,
    simulate_matrix,
    simulate_multiple_matrices,
    parse_args,
    main,
    is_file,)


SEQGEN_PATH = 'seq-gen'

TESTFILES_DIR = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), 'test_files')


def seqgen_status(path):
    """
    Return True if Seq-Gen executable is working,
    otherwise return False.
    """
    f = open(os.devnull, 'w')
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


class TestReadTreeFile():

    t_file_string = (
        """
        #NEXUS
        [ID: 9409050143]
        [Param: tree]
        begin trees;
            translate
            1 t1,
            2 t2,
            3 t3,
            4 t4;
        tree rep.1 = ((1:0.1,2:0.1):0.1,3:0.1,4:0.1);
        tree rep.2 = ((1:0.1,2:0.1):0.1,3:0.1,4:0.1);
        end;""")

    def test_read_tree_file(self, tmpdir):
        f = tmpdir.join('t-file.txt')
        f.write(self.t_file_string)
        tree_list = read_tfile(str(f.dirpath('t-file.txt')))
        assert len(tree_list) == 2

    def test_read_empty_tree_file(self, tmpdir):
        f = tmpdir.join('empty-t-file.txt')
        f.write('')
        with pytest.raises(ValueError):
            read_pfile(str(f.dirpath('empty-t-file.txt')))


class TestReadParameterFile():

    p_file_string = (
        '[ID: 1234567890]\n'
        'Gen\tLnL\tLnPr\tTL\tkappa\tpi(A)\tpi(C)\tpi(G)\tpi(T)\talpha\n'
        '0\t-5\t50\t0.5\t1\t0.25\t0.25\t0.25\t0.25\t1\n'
        '500\t-5\t50\t0.5\t1\t0.25\t0.25\t0.25\t0.25\t1')

    def test_read_parameter_file(self, tmpdir):
        f = tmpdir.join('p-file.txt')
        f.write(self.p_file_string)
        p_dicts = read_pfile(str(f.dirpath('p-file.txt')))
        assert len(p_dicts) == 2

    def test_read_empty_parameter_file(self, tmpdir):
        f = tmpdir.join('empty-p-file.txt')
        f.write('')
        with pytest.raises(ValueError):
            read_pfile(str(f.dirpath('empty-p-file.txt')))


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

    d1 = {'pi(A)': 0.25, 'pi(C)': 0.25, 'pi(G)': 0.25, 'pi(T)': 0.25}
    d2 = {'state_freqs': '0.25,0.25,0.25,0.25'}
    d3 = {'pi(A)': 0.25, 'pi(C)': 0.25, 'pi(G)': 0.25}

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
        result = simulate_matrix(self.tree, seqgen_path=SEQGEN_PATH)
        assert len(result.char_matrix) == 4
        assert result.char_matrix.sequence_size == 1000
        assert result.tree == self.tree
        assert 'HKY' in result.command_line

    def test_gtr(self):
        result = simulate_matrix(
            self.tree, general_rates=self.rates, seqgen_path=SEQGEN_PATH)
        assert 'GTR' in result.command_line

    def test_ti_tv(self):
        result = simulate_matrix(self.tree, ti_tv=1, seqgen_path=SEQGEN_PATH)
        assert 'HKY' in result.command_line
        assert ' -t1 ' in result.command_line

    def test_ti_tv_and_gtr(self):
        with pytest.raises(ValueError):
            simulate_matrix(
                self.tree, ti_tv=1, general_rates=self.rates,
                seqgen_path=SEQGEN_PATH)

    def test_gamma(self):
        result = simulate_matrix(
            self.tree, gamma_shape=2, seqgen_path=SEQGEN_PATH)
        assert ' -a2 ' in result.command_line

    def test_gamma_shape_and_cats(self):
        result = simulate_matrix(
            self.tree, gamma_shape=2, gamma_cats=5, seqgen_path=SEQGEN_PATH)
        assert ' -a2 ' in result.command_line
        assert ' -g5 ' in result.command_line

    def test_gamma_cats(self):
        with pytest.raises(ValueError):
            simulate_matrix(
                self.tree, gamma_cats=5, seqgen_path=SEQGEN_PATH)


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
            self.treelist, self.p_dicts, self.rng_seeds,
            seqgen_path=SEQGEN_PATH)

    def test_empty_input(self):
        with pytest.raises(AssertionError):
            simulate_multiple_matrices(
                dendropy.TreeList(), [], seqgen_path=SEQGEN_PATH)

    def test_parameter_mismatch(self):
        with pytest.raises(AssertionError):
            simulate_multiple_matrices(
                self.treelist[:2], self.p_dicts[:1], seqgen_path=SEQGEN_PATH)

    def test_rng_seeds_mismatch(self):
        with pytest.raises(AssertionError):
            simulate_multiple_matrices(
                self.treelist, self.p_dicts, rng_seeds=['123321'],
                seqgen_path=SEQGEN_PATH)


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

    def test_is_file(self):
        with tempfile.NamedTemporaryFile() as tmp:
            assert is_file(tmp.name) == tmp.name

    def test_is_file_error(self):
        with pytest.raises(argparse.ArgumentTypeError):
            is_file('')


@seqgen_required
class TestMain():

    def test_args_help(self):
        with pytest.raises(SystemExit):
            main(['-h'])

    def test_noargs(self):
        with pytest.raises(SystemExit):
            main()
