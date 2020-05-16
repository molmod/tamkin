import unittest
from get_energy_orca import get_energy_orca


class TestOrcaSPE(unittest.TestCase):

    def test_orca_spe(self):
        self.assertEqual(get_energy_orca("hno_spe.out"), -130.348214495060, "Should be -130.348214495060")

if __name__ == '__main__':
    unittest.main()

