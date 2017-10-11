from module import JsonHelper, ResultAnalyzer
import unittest

class TestModules(unittest.TestCase):
    
    def setUp(self):
        pass

    def test_resultAnalyzer(self):
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
        self.assertEqual(ResultAnalyzer.getAntigen(ResultAnalyzer.getSerotypes('asdH3333fsdfO222H3333')), 'O')
        self.assertEqual(ResultAnalyzer.getAntigen(ResultAnalyzer.getSerotypes('O')), '')

if __name__ == '__main__':
    unittest.main()