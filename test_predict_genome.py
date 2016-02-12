import unittest
import json

from predict_genome import svr_features_from_sequence

class TestPredictGenome(unittest.TestCase):

    def test_svr_features_from_sequence_1mer(self):
        sequence = 'AGC'
        features = svr_features_from_sequence(sequence, [1])
        self.assertEqual(len(features), 12, "Feature matrix size mismatch") # with 4 nucleotides and 3 positions there should be 12 features
        just_matches = filter(lambda x: x['value'] == 1, features)
        self.assertEqual(len(just_matches), 3, "Features of value 1 size mismatch") # For 1mers, there should only be 3 items with value 1
        positions = map(lambda x: x['position'], just_matches)
        just_nucleotides = map(lambda x: x['feature'], just_matches) # Extract the original sequence
        self.assertEqual(''.join(just_nucleotides), sequence, "Reconstructed sequence does not match")

    def test_svr_features_from_sequence_1_2_3mer(self):
        sequence = 'CAGGCTTTGGGAGCCAGCGGGGCGGGAGCGGCGAAG'
        features = svr_features_from_sequence(sequence, [1,2,3])
         # with 4 nucleotides and 1,2,3 positions, total should be 4^1*36 + 4^2 * 35 + 4^3 * 34 = 2880
        self.assertEqual(len(features), 2880, "Feature matrix size mismatch")

        just_matches = filter(lambda x: x['value'] == 1, features)
        print len(just_matches)
        # Should be 36+35+34 = 105 items where value is 1
        self.assertEqual(len(just_matches), 105, "Features of value 1 size mismatch")
        positions = map(lambda x: x['position'], just_matches)

        # Extract the original sequence, just the 1mers
        just_nucleotides = map(lambda x: x['feature'], just_matches[0:36])
        self.assertEqual(''.join(just_nucleotides), sequence, "Reconstructed sequence does not match")

    def test_svr_features_from_sequence_1_2_3mer_exact(self):
        sequence = 'CAGGCTTTGGGAGCCAGCGGGGCGGGAGCGGCGAAG'
        features = svr_features_from_sequence(sequence, [1,2,3])
        with open('test_matrix.json', 'r') as f:
            expected = json.load(f)
        self.assertEqual(features, expected, "Generated feature matrix does not match")


if __name__ == '__main__':
    unittest.main()