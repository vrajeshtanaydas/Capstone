import unittest
import Pipeline
import TreeMaker
import ISGDataPuller
import SequenceFinder
import EntropyCalculator
import SequenceDifferentiator
import SequenceToFasta

class TestTreeMakerMethods(unittest.TestCase):

    def setUp(self):
        fo = open('testTreeMaker.txt', 'w')
    
    def test_preprocess_boundary(self):
        fo.write('(a:0.5,b:0.7)')
        self.assertEqual(TreeMaker.TreeTable(fo), ['(', 'a', 'b', ')'])
    
    def test_preprocess_wrong_1(self):
        fo.write('{a:0.5,b:0.7}')
        self.assertRaises(IndexError, TreeMaker.TreeTable(), fo)
    
    def tearDown(self):
        fo.close()