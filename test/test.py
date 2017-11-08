import logging
import os
import sys
import unittest

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from module import JsonHelper, ResultAnalyzer

class TestModules(unittest.TestCase):
    
    def setUp(self):
        pass

    def test_getGeneName(self):
        # one gene -> return the gene
        self.assertEqual(ResultAnalyzer.getGeneName(
            "8__wzx__wzx-O157__188 AB602252.1;O antigen flippase;O157"),
            "wzx")
        # no gene -> return ""
        self.assertEqual(ResultAnalyzer.getGeneName(
            "8__wsx__wsx-O157__188 AB602252.1;O antigen flippase;O157"),
            "")
        # more than one gene -> return first in list
        self.assertEqual(ResultAnalyzer.getGeneName(
            "8__wzt__wzm-O157__188 AB602252.1;O antigen flippase;O157"),
            "wzm")
    
    def test_filter_lone_pair(self):
        # null test
        self.assertEqual(ResultAnalyzer.filter_lone_pair({}), {})
        # single gene -> filtered
        self.assertEqual(ResultAnalyzer.filter_lone_pair_helper(
            {'O157':[{'gene':'wzx', 'seq':'ATCG'}]}, 'wzx'), 
            {'O157': []})
        # single gene -> filtered
        self.assertEqual(ResultAnalyzer.filter_lone_pair(
            {'O157':[{'gene':'wzx', 'seq':'ATCG'}]}),
            {'O157': []})
        # pair gene -> not filtered
        self.assertEqual(ResultAnalyzer.filter_lone_pair(
            {'O157': [{'gene': 'wzx', 'seq': 'ATCG'}, {'gene': 'wzy', 'seq': 'ATCG'}]}),
            {'O157': [{'gene': 'wzx', 'seq': 'ATCG'}, {'gene': 'wzy', 'seq': 'ATCG'}]})
        # real data
        data_in = 
    
    def test_getSerotypes(self):
        self.assertEqual(ResultAnalyzer.getAntigen(ResultAnalyzer.getSerotypes('asdH3333fsdfO222H3333')), 'O')
        self.assertEqual(ResultAnalyzer.getAntigen(ResultAnalyzer.getSerotypes('O')), '')

    def test_expand_serotype_dict(self):
        self.assertEqual(ResultAnalyzer.expand_serotype_dict('/home/sam/Projects/MoreSerotype/test/Data/genome_output.json'), 597+1)

    def test_merge_blacklist_entry(self):
        entry1 = {
            'genome name': 'genomeA',
            'allele name': 'alleleA',
            'given O': 1,
            'given H': 2,
            'predicted O': '3',
            'predicted H': '4'
        }
        entry2 = {
            'genome name': 'genomeA',
            'allele name': 'alleleB',
            'given O': 1,
            'given H': 2,
            'predicted O': '5',
            'predicted H': '6'
        }
        merged_entry = {
            'genome name': 'genomeA',
            'allele name': 'alleleA|alleleB',
            'given O': 1,
            'given H': 2,
            'predicted O': '3|5',
            'predicted H': '4|6'
        }
        self.assertEqual(ResultAnalyzer.merge_blacklist_entry(entry1, entry2), merged_entry)

if __name__ == '__main__':
    unittest.main()
