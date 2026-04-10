from __future__ import annotations

import unittest

from workflows.smk_scripts._repeatmasker_out import parse_repeat_tags


class TestSmkScriptsConsolidations(unittest.TestCase):
    def test_repeatmasker_tag_parsing_old_and_new_shapes(self) -> None:
        old = "SMPL=KA1|ITER=1|RANK=2|ORIG=CL0001_some"
        new = "SMPL=KA1|ORG=Triticum|GENOMES=AA|ITER=1|RANK=2|ORIG=CL0001_some"

        d_old = parse_repeat_tags(old)
        self.assertEqual(d_old["smpl"], "KA1")
        self.assertEqual(d_old["iter"], "1")
        self.assertEqual(d_old["rank"], "2")
        self.assertEqual(d_old["orig"], "CL0001_some")

        d_new = parse_repeat_tags(new)
        self.assertEqual(d_new["smpl"], "KA1")
        self.assertEqual(d_new["iter"], "1")
        self.assertEqual(d_new["rank"], "2")
        self.assertEqual(d_new["orig"], "CL0001_some")


if __name__ == "__main__":
    unittest.main()

