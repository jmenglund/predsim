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
    combine_simulation_input,
    iter_seqgen_results,
    parse_args,
    main,
    is_file,)


SEQGEN_PATH = 'seq-gen'

TESTFILES_DIR = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), 'test_files')


with open(os.path.join(TESTFILES_DIR, 'expected-hky-3.nex'), 'r') as fo:
    EXP_NEX_3 = fo.read()

with open(os.path.join(TESTFILES_DIR, 'expected-hky-1.nex'), 'r') as fo:
    EXP_NEX_1 = fo.read()

with open(os.path.join(TESTFILES_DIR, 'expected-hky-3.phy'), 'r') as fo:
    EXP_PHY_3 = fo.read()

with open(os.path.join(TESTFILES_DIR, 'expected-hky-1.phy'), 'r') as fo:
    EXP_PHY_1 = fo.read()

with open(os.path.join(TESTFILES_DIR, 'expected-jc-1.nex'), 'r') as fo:
    EXP_JC_NEX_1 = fo.read()

with open(os.path.join(TESTFILES_DIR, 'expected-jc-gamma-1.nex'), 'r') as fo:
    EXP_JC_GAMMA_NEX_1 = fo.read()

with open(os.path.join(
        TESTFILES_DIR, 'expected-jc-propinvar-1.nex'), 'r') as fo:
    EXP_JC_PROPINVAR_NEX_1 = fo.read()

with open(os.path.join(
        TESTFILES_DIR, 'expected-jc-invgamma-1.nex'), 'r') as fo:
    EXP_JC_INVGAMMA_NEX_1 = fo.read()

with open(os.path.join(
        TESTFILES_DIR, 'expected-gtr-1.nex'), 'r') as fo:
    EXP_GTR_NEX_1 = fo.read()


def seqgen_status(path):
    """
    Return True if Seq-Gen executable is working,
    otherwise return False.
    """
    fo = open(os.devnull, 'w')
    try:
        subprocess.check_call(
            SEQGEN_PATH, stdout=fo, stderr=subprocess.STDOUT)
        status = True
    except subprocess.CalledProcessError:
        status = False
    except OSError:
        status = False
    finally:
        fo.close()
    return status


seqgen_required = pytest.mark.skipif(
    seqgen_status(SEQGEN_PATH) is False, reason='Seq-Gen is required')


class TestReadTreeFile():

    def test_read_tree_file(self):
        tree_list = read_tfile(os.path.join(TESTFILES_DIR, 'data-hky.t'))
        assert len(tree_list) == 3

    def test_read_empty_tree_file(self, tmpdir):
        fo = tmpdir.join('empty-t-file.txt')
        fo.write('')
        with pytest.raises(ValueError):
            read_pfile(str(fo.dirpath('empty-t-file.txt')))


class TestReadParameterFile():

    def test_read_parameter_file(self):
        p_dicts = read_pfile(os.path.join(TESTFILES_DIR, 'data-hky.p'))
        assert len(p_dicts) == 3

    def test_read_empty_parameter_file(self, tmpdir):
        fo = tmpdir.join('empty-p-file.txt')
        fo.write('')
        with pytest.raises(ValueError):
            read_pfile(str(fo.dirpath('empty-p-file.txt')))


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
        assert get_seqgen_params(self.d3) == self.d2


@seqgen_required
class TestSingleSimulation():

    tree_string = '((t1:0,t2:0):0,t3:0,t4:0);'
    tree = dendropy.Tree.get_from_string(tree_string, 'newick')

    rates = '1,1,1,1,1,1'

    def test_simulate_empty_params(self):
        result = simulate_matrix(self.tree, seqgen_path=SEQGEN_PATH)
        assert len(result.char_matrix) == 4
        assert result.char_matrix.sequence_size == 1000
        assert result.tree == self.tree.as_string('newick')
        assert 'HKY' in result.command

    def test_gtr(self):
        result = simulate_matrix(
            self.tree, general_rates=self.rates, seqgen_path=SEQGEN_PATH)
        assert 'GTR' in result.command

    def test_ti_tv(self):
        result = simulate_matrix(self.tree, ti_tv=1, seqgen_path=SEQGEN_PATH)
        assert 'HKY' in result.command
        assert ' -t1 ' in result.command

    def test_ti_tv_and_gtr(self):
        with pytest.raises(ValueError):
            simulate_matrix(
                self.tree, ti_tv=1, general_rates=self.rates,
                seqgen_path=SEQGEN_PATH)

    def test_gamma(self):
        result = simulate_matrix(
            self.tree, gamma_shape=2, seqgen_path=SEQGEN_PATH)
        assert ' -a2 ' in result.command

    def test_gamma_shape_and_cats(self):
        result = simulate_matrix(
            self.tree, gamma_shape=2, gamma_cats=5, seqgen_path=SEQGEN_PATH)
        assert ' -a2 ' in result.command
        assert ' -g5 ' in result.command

    def test_gamma_cats(self):
        with pytest.raises(ValueError):
            simulate_matrix(
                self.tree, gamma_cats=5, seqgen_path=SEQGEN_PATH)


class TestCombineSimulationInput():

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

    def test_input(self):
        combine_simulation_input(
            self.treelist, self.p_dicts, self.rng_seeds)

    def test_empty_input(self):
        with pytest.raises(AssertionError):
            combine_simulation_input(dendropy.TreeList(), [])

    def test_parameter_mismatch(self):
        with pytest.raises(AssertionError):
            combine_simulation_input(
                self.treelist[:2], self.p_dicts[:1])

    def test_rng_seeds_mismatch(self):
        with pytest.raises(AssertionError):
            combine_simulation_input(
                self.treelist, self.p_dicts, rng_seeds=['123321'])


@seqgen_required
class TestIterSeqgenResults():

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
    simulation_input = zip(treelist, p_dicts, rng_seeds)

    def test_multiple_simulations(self):
        results = iter_seqgen_results(
            self.simulation_input, seqgen_path=SEQGEN_PATH)
        assert len(list(results)) == 2


class TestArgumentParser():

    commands_fo = tempfile.NamedTemporaryFile('w')
    trees_fo = tempfile.NamedTemporaryFile('w')

    def test_parser_help(self):
        with pytest.raises(SystemExit):
            parse_args(['-h'])

    def test_parser(self):
        pfile_path = os.path.join(TESTFILES_DIR, 'data-hky.p')
        tfile_path = os.path.join(TESTFILES_DIR, 'data-hky.t')
        seeds_filepath = os.path.join(TESTFILES_DIR, 'seeds-1.txt')

        parser = parse_args([
            '-l', '2', '-s', '1', '-g', '5', '-n', '1', '-f', 'phylip',
            '-p', 'sg', '--seeds-file', seeds_filepath,
            '--commands-file', self.commands_fo.name,
            '--trees-file', self.trees_fo.name,
            pfile_path, tfile_path])
        assert parser.length == 2
        assert parser.skip == 1
        assert parser.gamma_cats == 5
        assert parser.num_records == 1
        assert parser.out_format == 'phylip'
        assert parser.sg_filepath == 'sg'
        assert parser.commands_filepath == self.commands_fo.name
        assert parser.trees_filepath == self.trees_fo.name
        assert parser.seeds_filepath == seeds_filepath
        assert parser.pfile_path == pfile_path
        assert parser.tfile_path == tfile_path

    def test_is_file(self):
        with tempfile.NamedTemporaryFile() as tmp:
            assert is_file(tmp.name) == tmp.name

    def test_is_file_error(self):
        with pytest.raises(argparse.ArgumentTypeError):
            is_file('')


@seqgen_required
class TestMain():

    outfile = tempfile.NamedTemporaryFile('w')

    def test_args_help(self):
        with pytest.raises(SystemExit):
            main(['-h'])

    def test_noargs(self):
        with pytest.raises(SystemExit):
            main()

    def test_hky(self):
        main([
            '-l', '2',
            os.path.join(TESTFILES_DIR, 'data-hky.p'),
            os.path.join(TESTFILES_DIR, 'data-hky.t')])

    def test_hky_seeds_file(self, capsys):
        main([
            '-l', '2',
            '--seeds-file', os.path.join(TESTFILES_DIR, 'seeds-3.txt'),
            os.path.join(TESTFILES_DIR, 'data-hky.p'),
            os.path.join(TESTFILES_DIR, 'data-hky.t')])
        out, err = capsys.readouterr()
        assert out == EXP_NEX_3
        assert err == ''

    def test_hky_skip(self, capsys):
        main([
            '-l', '2', '-s', '2',
            '--seeds-file', os.path.join(TESTFILES_DIR, 'seeds-1.txt'),
            os.path.join(TESTFILES_DIR, 'data-hky.p'),
            os.path.join(TESTFILES_DIR, 'data-hky.t')])
        out, err = capsys.readouterr()
        assert out == EXP_NEX_1
        assert err == ''

    def test_hky_commands_file(self):
        main([
            '-l', '2', '-s', '2', '--commands-file', self.outfile.name,
            '--seeds-file', os.path.join(TESTFILES_DIR, 'seeds-1.txt'),
            os.path.join(TESTFILES_DIR, 'data-hky.p'),
            os.path.join(TESTFILES_DIR, 'data-hky.t')])
        assert os.path.isfile(self.outfile.name)

    def test_hky_trees_file(self):
        main([
            '-l', '2', '-s', '2', '--trees-file', self.outfile.name,
            '--seeds-file', os.path.join(TESTFILES_DIR, 'seeds-1.txt'),
            os.path.join(TESTFILES_DIR, 'data-hky.p'),
            os.path.join(TESTFILES_DIR, 'data-hky.t')])
        assert os.path.isfile(self.outfile.name)

    def test_hky_phylip(self, capsys):
        main([
            '-l', '2', '-s', '2', '-f', 'phylip',
            '--seeds-file', os.path.join(TESTFILES_DIR, 'seeds-1.txt'),
            os.path.join(TESTFILES_DIR, 'data-hky.p'),
            os.path.join(TESTFILES_DIR, 'data-hky.t')])
        out, err = capsys.readouterr()
        assert out == EXP_PHY_1
        assert err == ''

    @pytest.mark.parametrize(
        'pfile,tfile,expected',
        [
            ('data-jc.p', 'data-jc.t', 'expected-jc-1.nex'),
            ('data-jc-gamma.p', 'data-jc-gamma.t', 'expected-jc-gamma-1.nex'),
            (
                'data-jc-propinvar.p', 'data-jc-propinvar.t',
                'expected-jc-propinvar-1.nex'),
            (
                'data-jc-invgamma.p', 'data-jc-invgamma.t',
                'expected-jc-invgamma-1.nex'),
            ('data-gtr.p', 'data-gtr.t', 'expected-gtr-1.nex'),
        ])
    def test_subst_model(self, capsys, pfile, tfile, expected):
        with open(os.path.join(TESTFILES_DIR, expected), 'r') as fo:
            expected_content = fo.read()
        main([
            '-l', '2', '-n', '1', '--seeds-file',
            os.path.join(TESTFILES_DIR, 'seeds-1.txt'),
            os.path.join(TESTFILES_DIR, pfile),
            os.path.join(TESTFILES_DIR, tfile)])
        out, err = capsys.readouterr()
        assert out == expected_content
        assert err == ''
