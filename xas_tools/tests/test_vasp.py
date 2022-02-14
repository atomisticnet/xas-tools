import unittest
import tempfile
import json
import os

import pymatgen as mg

from xas_tools.vasp import CHPCalculation

fixtures = os.path.join(os.path.dirname(__file__), 'fixtures')
Li3PS4 = os.path.join(fixtures, 'CONTCAR_Li3PS4')


class VaspTest(unittest.TestCase):

    def test_write_vasp_Li3PS4(self):
        struc = mg.Structure.from_file(Li3PS4)
        chp = CHPCalculation(struc, element="S", n=1, ell=0, z=1.0)
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'Li3PS4_chp_input')
            chp.write_vasp_input(supercell=(1, 2, 1),
                                 band_multiple=2,
                                 path=path)
            with open(os.path.join(path, 'metadata.json')) as fp:
                metadata = json.load(fp)
            self.assertEqual(metadata['multiplicity'][0], 2)

    def test_vasp_warning(self):
        raise NotImplementedError(
            'The VASP test is only a stub and needs to be completed.')


if __name__ == "__main__":
    unittest.main()
