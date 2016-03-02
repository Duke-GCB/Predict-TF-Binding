import unittest
import json

from predict_genome import svr_features_from_sequence, generate_matching_sequences, predict, load_model

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

    def _check_generated_matches(self, sequence, core, width, expected_matches):
        count = 0
        for match in generate_matching_sequences(sequence, core, width):
            count+=1
            self.assertIn(match, expected_matches, 'Mismatch : {} not in {}'.format(match, expected_matches))
        self.assertEqual(count, len(expected_matches), 'Unexpected number of matches')

    def test_generates_matching_sequences(self):
        sequence = 'ACCTTAGCCTTGATAT'
        core = 'CCTT'
        expected_matches = [(0,('ACCTTA',)),(6,('GCCTTG',))]
        self._check_generated_matches(sequence, core, 6, expected_matches)

    def test_skips_unknown_bases(self):
        sequence = 'ACCTTAGCCTTNATAT'
        core = 'CCTT'
        expected_matches = [(0,('ACCTTA',))]
        self._check_generated_matches(sequence, core, 6, expected_matches)

    def load_model(self):
        # This model is 178MB in size, so not practical to include in the source repo
        model_dict = load_model('ELK1_100nM_Bound_filtered_normalized_GGAA_1a2a3mer_format.model')
        self.assertEqual(model_dict['size'], 2881)

    def test_generates_matching_reverse_complements(self):
        core = 'GCTG' # must not be palindromic
        sequence = 'ATTCAGCGAA' # Reverse complement in the middle
        expected_matches = [(2,('CGCTGA',))]
        self._check_generated_matches(sequence, core, 6, expected_matches)

    def test_generates_matching_palindromes(self):
        core = 'GGCC' # Must be palindromic
        sequence = 'ATTGGCCGAA' # Core in the middle
        expected_matches = [(2, ('TGGCCG', 'CGGCCA'))]
        # Palindromes yield a position and two sequences
        self._check_generated_matches(sequence, core, 6, expected_matches)


if __name__ == '__main__':
    unittest.main()