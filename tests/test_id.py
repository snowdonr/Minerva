'''
Created on Nov. 4, 2022

@author: Ryan
'''
import unittest
from gcms import identification
from gcms import mass_spectrum


class IdentificationTest(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_add_spectrum(self):
        test_obj = identification.UserEntry()
        self.assertIsInstance(test_obj.spectrum, mass_spectrum.MassSpectrum_SQL)
        self.assertTrue(len(test_obj.spectrum.values) == 0)
        test_obj.mass1(50)
        self.assertTrue(len(test_obj.spectrum.values) == 1)
        

    def test_pub_chem(self):
        id_type_dict = {"name": ("ethanol", 702, "ethanol"),
                        "inchi": (r"InChI=1S/C2H8O2Si/c1-5(2,3)4/h3-4H,1-2H3", 14014, "dihydroxy(dimethyl)silane"),
                        "inchikey": ("LYBIZMNPXTXVMV-UHFFFAOYSA-N", 12716, "propan-2-yl prop-2-enoate"),
                        "cid": ("1", 1, "3-acetyloxy-4-(trimethylazaniumyl)butanoate")}
        for id_type, values in id_type_dict.items():
            result = identification.UserEntry.scan_pubchem(id_type, values[0])
            self.assertTrue(len(result) > 0, msg="At least one result returned by pubchem")
            self.assertEqual(result[0]["CID"], values[1])
            self.assertEqual(result[0]["IUPACName"], values[2])
            self.assertTrue(len(result[0]["InChI"]) > 0)
            # "InChI", "InChIKey", "MolecularFormula", "MolecularWeight", "MonoisotopicMass", "XLogP", "CanonicalSMILES"


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
