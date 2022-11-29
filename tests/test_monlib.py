"""
Author: "Keitaro Yamashita, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""
from __future__ import absolute_import, division, print_function, generators
import unittest
import os
import glob
from servalcat import utils
import gemmi
import common

monlib_path = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

print("monomer library path= {}".format(monlib_path))
doc = gemmi.cif.read(os.path.join(monlib_path, "list/mon_lib_list.cif"))
comp_list = doc.find_block("comp_list")
link_list = doc.find_block("link_list")
mod_list = doc.find_block("mod_list")
lgroups = {row.str(0):row.str(1) for row in comp_list.find("_chem_comp.", ["id", "group"])}
print("Loading {} monomers".format(len(lgroups)))
monlib = gemmi.read_monomer_lib(monlib_path, list(lgroups))

def check_mpeptide(cc): # this returns True for P-peptide as well!
    mc_atoms = set(['N', "H", 'CA', 'C', 'O', 'OXT']) # "HA", 
    std_bonds = set([tuple(sorted(b)) for b in (("N","H"), ("C", "O"), ("N", "CA"), ("CA", "C"), ("C", "OXT"))]) # ("CA", "HA"),
    
    atoms = set([a.id for a in cc.atoms])
    if not mc_atoms.issubset(atoms):
        return False

    bonds = set([tuple(sorted((b.id1.atom, b.id2.atom))) for b in cc.rt.bonds])
    if not std_bonds.issubset(bonds):
        return False
    
    n_bonds = set([b[0] for b in bonds if b[1]=="N"])
    #if n_bonds.issubset(("H","H2","H3")):
    if n_bonds == set(("H","H2","H3","CA")):
        return False # this is normal peptide
    
    for r in n_bonds - set(("H", "CA")):
        #if r[0] != "H": return True
        if r == "CN": return "M"
        if r == "CD": return "P"
    return False


def test_dna_rna(cc):
    na_atoms = set(["C4'", "C3'", "O3'", "HO3'", "C5'", "O5'", "P", "OP3", "OP2", "OP1"])
    na_tor1 = ["C4'", "C3'", "O3'", "HO3'"]
    na_tor2 = ["C5'", "O5'", "P", "OP3"]

    atoms = set([a.id for a in cc.atoms])
    if not na_atoms.issubset(atoms):
        return False

    # Check 3' atoms and 5' atoms
    flags = [False, False]
    for t in cc.rt.torsions:
        t_atoms = [t.id1.atom, t.id2.atom, t.id3.atom, t.id4.atom]
        if na_tor1 == t_atoms or na_tor1 == t_atoms[::-1]:
            flags[0] = True
        if na_tor2 == t_atoms or na_tor2 == t_atoms[::-1]:
            flags[1] = True
    #print(flags)
    if not all(flags):
        return False

    # Check around P
    # P-OP3 (single), P-OP1, P-OP2 bond
    flags = [False, False, False]
    for b in cc.rt.bonds:
        if sorted([b.id1.atom, b.id2.atom]) == sorted(["P", "OP3"]) and b.type == gemmi.BondType.Single:
            flags[0] = True
        elif sorted([b.id1.atom, b.id2.atom]) == sorted(["P", "OP1"]):
            flags[1] = True
        elif sorted([b.id1.atom, b.id2.atom]) == sorted(["P", "OP2"]):
            flags[2] = True
    #print(flags)
    if not all(flags):
        return False
    return True


class TestMonlib(unittest.TestCase):
    def setUp(self): self.errors = []
    def tearDown(self): self.assertEqual(len(self.errors), 0, msg="\n"+"\n".join(self.errors))

    def test_list_consistency(self):
        cgroups = {}
        for f in glob.glob(os.path.join(monlib_path, "?/*.cif")):
            tmp = gemmi.cif.read(f)
            block = tmp.find_block("comp_list")
            table = block.find("_chem_comp.", ["id", "group"])
            for row in table:
                cgroups[row[0]] = row.str(1)

        only_in_cif = set(cgroups) - set(lgroups)
        try: self.assertFalse(only_in_cif, msg="only in cif files")
        except AssertionError as e: self.errors.append(str(e))
        
        only_in_list = set(lgroups) - set(cgroups)
        try: self.assertFalse(only_in_list, msg="only in list")
        except AssertionError as e: self.errors.append(str(e))

    def test_group_cif(self):
        for mon in monlib.monomers:
            cc = monlib.monomers[mon]
            try: self.assertEqual(lgroups[mon], cc.group, msg="group mismatch: {}".format(mon))
            except AssertionError as e: self.errors.append(str(e))

    def test_group(self):
        lg_set = set(lgroups.values())
        all_set = set(["peptide", "P-peptide", "M-peptide", "DNA", "RNA", "pyranose", "ketopyranose", "furanose", "NON-POLYMER"])
        try: self.assertTrue(lg_set.issubset(all_set), msg="groups: {}".format(str(lg_set)))
        except AssertionError as e: self.errors.append(str(e))

    def test_links(self):
        known_groups = set(lgroups.values())
        known_groups.add("DNA/RNA") # can we really consider these known groups?
        known_groups.add("pept")
        
        ltab = link_list.find("_chem_link.", ["id", "comp_id_1", "mod_id_1", "group_comp_1", "comp_id_2", "mod_id_2", "group_comp_2"])
        mtab = mod_list.find("_chem_mod.", ["id", "comp_id", "group_id"])
        
        link_undef_group = [(row.str(0), row.str(i)) for row in ltab for i in (3,6) if row.str(i) and row.str(i) not in known_groups]
        mod_undef_group =  [(row.str(0), row.str(2)) for row in mtab if row.str(2) and row.str(2) not in known_groups]
        link_undef_comp =  [(row.str(0), row.str(i)) for row in ltab for i in (1,4) if row.str(i) !="" and row.str(i) not in lgroups]
        mod_undef_comp =   [(row.str(0), row.str(1)) for row in mtab if row.str(1) !="" and row.str(1) not in lgroups]

        try: self.assertFalse(link_undef_group, msg="undefined groups referenced in links")
        except AssertionError as e: self.errors.append(str(e))
        try: self.assertFalse(mod_undef_group, msg="undefined groups referenced in mods")
        except AssertionError as e: self.errors.append(str(e))
        try: self.assertFalse(link_undef_comp, msg="undefined comp referenced in links")
        except AssertionError as e: self.errors.append(str(e))
        try: self.assertFalse(mod_undef_comp, msg="undefined comp referenced in mods")
        except AssertionError as e: self.errors.append(str(e))

        for row in ltab:
            if row.str(1):
                self.assertEqual(row.str(3), lgroups[row.str(1)], msg="{} in link {}".format(row.str(1), row.str(0)))
            if row.str(4):
                self.assertEqual(row.str(6), lgroups[row.str(4)], msg="{} in link {}".format(row.str(4), row.str(0)))

        for row in mtab:
            if row.str(1):
                self.assertEqual(row.str(2), lgroups[row.str(1)], msg="{} in mod {}".format(row.str(1), row.str(0)))

    def test_peptide(self): pass
    def test_sugars(self): pass
    def test_m_p_peptide(self):
        for mon in lgroups:
            if lgroups[mon] in ("P-peptide", "M-peptide"):
                check = check_mpeptide(monlib.monomers[mon])
                try: self.assertEqual(check, lgroups[mon][0], msg="{} not qualified as {}".format(mon, lgroups[mon]))
                except AssertionError as e: self.errors.append(str(e))
                
    def test_rna_dna(self):
        for mon in lgroups:
            if lgroups[mon] in ("DNA", "RNA"):
                check = test_dna_rna(monlib.monomers[mon])
                try: self.assertTrue(check, msg="{} not qualified as DNA/RNA".format(mon))
                except AssertionError as e: self.errors.append(str(e))

    def test_sp2sp2(self):
        for mon in monlib.monomers:
            cc = monlib.monomers[mon]
            for tor in cc.rt.torsions:
                if tor.label.startswith("sp2_sp2"):
                    self.assertAlmostEqual(tor.esd, 5, msg="{} {}".format(mon, tor.label))

        for name in monlib.links:
            cc = monlib.links[name]
            for tor in cc.rt.torsions:
                if tor.label.startswith("sp2_sp2"):
                    self.assertAlmostEqual(tor.esd, 5, msg="{} {}".format(name, tor.label))

if __name__ == '__main__':
    #TestServal().test_6kls()
    unittest.main()
test_dna_rna
