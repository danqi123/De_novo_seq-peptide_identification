"""Unit test for identifier."""

from project_spectra.identifier import Identifier

class TestIdentifier:
    """Test code for project_spectra.identifier"""

    def test_prot_identifier(self):
        """Test function prot_identifier"""
        peptide_list = ["FEFER", "FEEFR", "FKCCAR", "FKSCSR"]
        ID = Identifier(peptide_list)
        assert isinstance(ID.stout, dict)
        assert "M1FN61_SOYBN" in ID.stout.keys()




