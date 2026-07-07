import unittest

from PhenGO.predict.data import PhenotypeLabelEncoder


class PhenotypeLabelEncoderTests(unittest.TestCase):
    def test_lethal_is_positive_class(self):
        encoder = PhenotypeLabelEncoder()

        self.assertEqual(encoder.classes_, ["viable", "lethal"])
        self.assertEqual(
            encoder.fit_transform(["viable", "lethal", "inviable", "essential"]),
            [0, 1, 1, 1],
        )
        self.assertEqual(encoder.inverse_transform([0, 1, 1]), ["viable", "lethal", "lethal"])


if __name__ == "__main__":
    unittest.main()
