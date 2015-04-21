import unittest
import TreeMaker
import ISGDataPuller
import SequenceFinder
import EntropyCalculator
import SequenceDifferentiator
import SequenceToFasta


# Test the TreeTableMaker class
class TestTreeMakerMethods(unittest.TestCase):

    def setUp(self):
        self.fname = 'testTreeMaker.txt'
        self.fo = open(self.fname, 'w+')


    # Test the normal operation of preprocess()
    def test_preprocess_normal(self):
        self.fo.write('(a:0.5,b:0.7,(c:0.0002,d:0.1,(e:0.99)))')
        # close the file since TreeTable is going to open it
        self.fo.close()
        processed = TreeMaker.TreeTable(self.fname).processed

        self.assertEqual(processed, ['(', 'a', 'b', '(', 'c', 'd', '(', 'e',
            ')', ')', ')'])


    # Boundary condition test
    def test_preprocess_two_genomes(self):
        self.fo.write('(a:0.5,b:0.7)')
        # close the file since TreeTable is going to open it
        self.fo.close()
        processed = TreeMaker.TreeTable(self.fname).processed

        self.assertEqual(processed, ['(', 'a', 'b', ')'])
    

    # Test bad characters in file to preprocess
    def test_preprocess_wrong_1(self):
        self.fo.write('{a:0.5,b:0.7}')
        self.fo.close()

        with self.assertRaises(TreeMaker.NewickFileError) as cm:
            TreeMaker.TreeTable(self.fname)



    def tearDown(self):
        pass
    
if __name__ == '__main__':
    unittest.main()
