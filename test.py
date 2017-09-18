from module import HeaderModifier, JsonHelper, SerotypeHelper
import unittest

STRAINS_FILE = "data/strains.json"
ALLELE_FILE = "data/EcOH.fasta"
class TestModules(unittest.TestCase):
    
    def setUp(self):
        self.strains = JsonHelper.read_from_json(STRAINS_FILE)

    def test_headerModifier(self):
        with self.assertRaises(Warning):
            HeaderModifier.getNewHeader("", self.strains)
            HeaderModifier.getNewHeader("XXXXXX", self.strains)
        self.assertEqual(HeaderModifier.getNewHeader("ESC_BA3688AA_AS", self.strains), "ESC_BA3688AA_AS|O157:H7")

    def test_serotypeHelper(self):
        self.assertEqual(SerotypeHelper.getSerotype(""), "")
        self.assertEqual(SerotypeHelper.getSerotype("1__fliC__fliC-H10__3 AY337482.1;flagellin;H10"), "H10")
        self.assertEqual(SerotypeHelper.getSerotype("OH"), "")
        self.assertEqual(SerotypeHelper.getSerotype("Non-O157"), "")
        self.assertEqual(SerotypeHelper.getSerotype("O12345"), "")
        self.assertEqual(SerotypeHelper.getSerotypes(""), {'O': '', 'H': ''} )
        self.assertEqual(SerotypeHelper.getSerotypes("O123"), {'O': '123', 'H': ''} )
        self.assertEqual(SerotypeHelper.getSerotypes("H234"), {'O': '', 'H': '234'} )
        self.assertEqual(SerotypeHelper.getSerotypes("O123H234"), {'O': '123', 'H': '234'} )
        self.assertEqual(SerotypeHelper.getSerotypes("O123Non-H234"), {'O': '123', 'H': ''} )
        self.assertFalse(SerotypeHelper.isMismatch("", {'O': '', 'H': ''}))
        self.assertFalse(SerotypeHelper.isMismatch("H1", {'O': '1', 'H': '1'}))
        self.assertFalse(SerotypeHelper.isMismatch("H1", {'O': '1', 'H': ''}))
        self.assertTrue(SerotypeHelper.isMismatch("H1", {'O': '1', 'H': '2'}))
        self.assertTrue(SerotypeHelper.isMismatch("O2", {'O': '1', 'H': '2'}))
        self.assertTrue(SerotypeHelper.isSameClass("H1", {'O': '2', 'H': '1'}))
        self.assertFalse(SerotypeHelper.isSameClass("H1", {'O': '2', 'H': ''}))
        self.assertFalse(SerotypeHelper.isSameClass("O1", {'O': '', 'H': '123'}))


if __name__ == '__main__':
    unittest.main()