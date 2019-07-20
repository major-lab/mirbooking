import gi
gi.require_version('Mirbooking', '2.2')
from gi.repository import GLib, Mirbooking
import unittest

target_seq = """
GCACACAGAGCAGCATAAAGCCCAGTTGCTTTGGGAAGTGTTTGGGACCAGATGGATTGT
AGGGAGTAGGGTACAATACAGTCTGTTCTCCTCCAGCTCCTTCTTTCTGCAACATGGGGA
AGAACAAACTCCTTCATCCAAGTCTGGTTCTTCTCCTCTTGGTCCTCCTGCCCACAGACG
CCTCAGTCTCTGGAAAACCGCAGTATATGGTTCTGGTCCCCTCCCTGCTCCACACTGAGA
CCACTGAGAAGGGCTGTGTCCTTCTGAGCTACCTGAATGAGACAGTGACTGTAAGTGCTT
CCTTGGAGTCTGTCAGGGGAAACAGGAGCCTCTTCACTGACCTGGAGGCGGAGAATGACG
TACTCCACTGTGTCGCCTTCGCTGTCCCAAAGTCTTCATCCAATGAGGAGGTAATGTTCC
TCACTGTCCAAGTGAAAGGACCAACCCAAGAATTTAAGAAGCGGACCACAGTGATGGTTA
AGAACGAGGACAGTCTGGTCTTTGTCCAGACAGACAAATCAATCTACAAACCAGGGCAGA
CAGTGAAATTTCGTGTTGTCTCCATGGATGAAAACTTTCACCCCCTGAATGAGTTGATTC
CACTAGTATACATTCAGGATCCCAAAGGAAATCGCATCGCACAATGGCAGAGTTTCCAGT
TAGAGGGTGGCCTCAAGCAATTTTCTTTTCCCCTCTCATCAGAGCCCTTCCAGGGCTCCT
ACAAGGTGGTGGTACAGAAGAAATCAGGTGGAAGGACAGAGCACCCTTTCACCGTGGAGG
AATTTGTTCTTCCCAAGTTTGAAGTACAAGTAACAGTGCCAAAGATAATCACCATCTTGG
AAGAAGAGATGAATGTATCAGTGTGTGGCCTATACACATATGGGAAGCCTGTCCCTGGAC
ATGTGACTGTGAGCATTTGCAGAAAGTATAGTGACGCTTCCGACTGCCACGGTGAAGATT
CACAGGCTTTCTGTGAGAAATTCAGTGGACAGCTAAACAGCCATGGCTGCTTCTATCAGC
AAGTAAAAACCAAGGTCTTCCAGCTGAAGAGGAAGGAGTATGAAATGAAACTTCACACTG
AGGCCCAGATCCAAGAAGAAGGAACAGTGGTGGAATTGACTGGAAGGCAGTCCAGTGAAA
TCACAAGAACCATAACCAAACTCTCATTTGTGAAAGTGGACTCACACTTTCGACAGGGAA
TTCCCTTCTTTGGGCAGGTGCGCCTAGTAGATGGGAAAGGCGTCCCTATACCAAATAAAG
TCATATTCATCAGAGGAAATGAAGCAAACTATTACTCCAATGCTACCACGGATGAGCATG
GCCTTGTACAGTTCTCTATCAACACCACCAATGTTATGGGTACCTCTCTTACTGTTAGGG
TCAATTACAAGGATCGTAGTCCCTGTTACGGCTACCAGTGGGTGTCAGAAGAACACGAAG
AGGCACATCACACTGCTTATCTTGTGTTCTCCCCAAGCAAGAGCTTTGTCCACCTTGAGC
CCATGTCTCATGAACTACCCTGTGGCCATACTCAGACAGTCCAGGCACATTATATTCTGA
ATGGAGGCACCCTGCTGGGGCTGAAGAAGCTCTCCTTCTATTATCTGATAATGGCAAAGG
GAGGCATTGTCCGAACTGGGACTCATGGACTGCTTGTGAAGCAGGAAGACATGAAGGGCC
ATTTTTCCATCTCAATCCCTGTGAAGTCAGACATTGCTCCTGTCGCTCGGTTGCTCATCT
ATGCTGTTTTACCTACCGGGGACGTGATTGGGGATTCTGCAAAATATGATGTTGAAAATT
GTCTGGCCAACAAGGTGGATTTGAGCTTCAGCCCATCACAAAGTCTCCCAGCCTCACACG
CCCACCTGCGAGTCACAGCGGCTCCTCAGTCCGTCTGCGCCCTCCGTGCTGTGGACCAAA
GCGTGCTGCTCATGAAGCCTGATGCTGAGCTCTCGGCGTCCTCGGTTTACAACCTGCTAC
CAGAAAAGGACCTCACTGGCTTCCCTGGGCCTTTGAATGACCAGGACGATGAAGACTGCA
TCAATCGTCATAATGTCTATATTAATGGAATCACATATACTCCAGTATCAAGTACAAATG
AAAAGGATATGTACAGCTTCCTAGAGGACATGGGCTTAAAGGCATTCACCAACTCAAAGA
TTCGTAAACCCAAAATGTGTCCACAGCTTCAACAGTATGAAATGCATGGACCTGAAGGTC
TACGTGTAGGTTTTTATGAGTCAGATGTAATGGGAAGAGGCCATGCACGCCTGGTGCATG
TTGAAGAGCCTCACACGGAGACCGTACGAAAGTACTTCCCTGAGACATGGATCTGGGATT
TGGTGGTGGTAAACTCAGCAGGTGTGGCTGAGGTAGGAGTAACAGTCCCTGACACCATCA
CCGAGTGGAAGGCAGGGGCCTTCTGCCTGTCTGAAGATGCTGGACTTGGTATCTCTTCCA
CTGCCTCTCTCCGAGCCTTCCAGCCCTTCTTTGTGGAGCTCACAATGCCTTACTCTGTGA
TTCGTGGAGAGGCCTTCACACTCAAGGCCACGGTCCTAAACTACCTTCCCAAATGCATCC
GGGTCAGTGTGCAGCTGGAAGCCTCTCCCGCCTTCCTAGCTGTCCCAGTGGAGAAGGAAC
AAGCGCCTCACTGCATCTGTGCAAACGGGCGGCAAACTGTGTCCTGGGCAGTAACCCCAA
AGTCATTAGGAAATGTGAATTTCACTGTGAGCGCAGAGGCACTAGAGTCTCAAGAGCTGT
GTGGGACTGAGGTGCCTTCAGTTCCTGAACACGGAAGGAAAGACACAGTCATCAAGCCTC
TGTTGGTTGAACCTGAAGGACTAGAGAAGGAAACAACATTCAACTCCCTACTTTGTCCAT
CAGGTGGTGAGGTTTCTGAAGAATTATCCCTGAAACTGCCACCAAATGTGGTAGAAGAAT
CTGCCCGAGCTTCTGTCTCAGTTTTGGGAGACATATTAGGCTCTGCCATGCAAAACACAC
AAAATCTTCTCCAGATGCCCTATGGCTGTGGAGAGCAGAATATGGTCCTCTTTGCTCCTA
ACATCTATGTACTGGATTATCTAAATGAAACACAGCAGCTTACTCCAGAGATCAAGTCCA
AGGCCATTGGCTATCTCAACACTGGTTACCAGAGACAGTTGAACTACAAACACTATGATG
GCTCCTACAGCACCTTTGGGGAGCGATATGGCAGGAACCAGGGCAACACCTGGCTCACAG
CCTTTGTTCTGAAGACTTTTGCCCAAGCTCGAGCCTACATCTTCATCGATGAAGCACACA
TTACCCAAGCCCTCATATGGCTCTCCCAGAGGCAGAAGGACAATGGCTGTTTCAGGAGCT
CTGGGTCACTGCTCAACAATGCCATAAAGGGAGGAGTAGAAGATGAAGTGACCCTCTCCG
CCTATATCACCATCGCCCTTCTGGAGATTCCTCTCACAGTCACTCACCCTGTTGTCCGCA
ATGCCCTGTTTTGCCTGGAGTCAGCCTGGAAGACAGCACAAGAAGGGGACCATGGCAGCC
ATGTATATACCAAAGCACTGCTGGCCTATGCTTTTGCCCTGGCAGGTAACCAGGACAAGA
GGAAGGAAGTACTCAAGTCACTTAATGAGGAAGCTGTGAAGAAAGACAACTCTGTCCATT
GGGAGCGCCCTCAGAAACCCAAGGCACCAGTGGGGCATTTTTACGAACCCCAGGCTCCCT
CTGCTGAGGTGGAGATGACATCCTATGTGCTCCTCGCTTATCTCACGGCCCAGCCAGCCC
CAACCTCGGAGGACCTGACCTCTGCAACCAACATCGTGAAGTGGATCACGAAGCAGCAGA
ATGCCCAGGGCGGTTTCTCCTCCACCCAGGACACAGTGGTGGCTCTCCATGCTCTGTCCA
AATATGGAGCAGCCACATTTACCAGGACTGGGAAGGCTGCACAGGTGACTATCCAGTCTT
CAGGGACATTTTCCAGCAAATTCCAAGTGGACAACAACAACCGCCTGTTACTGCAGCAGG
TCTCATTGCCAGAGCTGCCTGGGGAATACAGCATGAAAGTGACAGGAGAAGGATGTGTCT
ACCTCCAGACATCCTTGAAATACAATATTCTCCCAGAAAAGGAAGAGTTCCCCTTTGCTT
TAGGAGTGCAGACTCTGCCTCAAACTTGTGATGAACCCAAAGCCCACACCAGCTTCCAAA
TCTCCCTAAGTGTCAGTTACACAGGGAGCCGCTCTGCCTCCAACATGGCGATCGTTGATG
TGAAGATGGTCTCTGGCTTCATTCCCCTGAAGCCAACAGTGAAAATGCTTGAAAGATCTA
ACCATGTGAGCCGGACAGAAGTCAGCAGCAACCATGTCTTGATTTACCTTGATAAGGTGT
CAAATCAGACACTGAGCTTGTTCTTCACGGTTCTGCAAGATGTCCCAGTAAGAGATCTGA
AACCAGCCATAGTGAAAGTCTATGATTACTACGAGACGGATGAGTTTGCAATTGCTGAGT
ACAATGCTCCTTGCAGCAAAGATCTTGGAAATGCTTGAAGACCACAAGGCTGAAAAGTGC
TTTGCTGGAGTCCTGTTCTCAGAGCTCCACAGAAGACACGTGTTTTTGTATCTTTAAAGA
CTTGATGAATAAACACTTTTTCTGGTCAATGTCAAAAAAAAAAAAAAAAAAAAAAAAA
"""

target = Mirbooking.Target(accession='NM_000014.4', sequence=target_seq)
mirna = Mirbooking.Mirna(accession='MIMAT0000001', sequence='UGAGGUAGUAGGUUGUAUAGUU')

def reverse_complement(seq):
    return reversed(seq.translate(bytes.maketrans(b'ACGU', b'TGCA')))

from math import exp
RT = 1.987203611e-3 * 310.15

class SimpleScoreTable(Mirbooking.ScoreTable):
    def do_compute_score(self, mirna, target, position):
        score = Mirbooking.Score()
        score.kf = 1
        score.kr = float('inf')
        score.kcat = 0

        # simple hamming distance
        if position > target.get_sequence_length() - 7:
            return True, score

        d = sum(1 if a == b else 0
                for a, b in zip(reverse_complement(mirna.get_subsequence (1, 7)), target.get_subsequence (position, 7)))

        # at most 2 mismatch
        if d >= 5:
            score.kr = 1e12 * exp(-d/RT)

        return True, score

class MirbookingScoreTableTestCase(unittest.TestCase):
    def test_simple_score_table(self):
        score_table = SimpleScoreTable()
        ret, positions = score_table.compute_positions(mirna, target)
        self.assertTrue(ret)
        for p in positions:
            self.assertTrue(score_table.compute_score(mirna, target, p)[1].kr < float('inf'))

class MirbookingDefaultScoreTableTestCase(unittest.TestCase):
    def test_cutoff_filter(self):
        broker = Mirbooking.Broker(score_table=SimpleScoreTable())
        broker.set_sequence_quantity(target, 5.0)
        broker.set_sequence_quantity(mirna, 5.0)

        score_table = Mirbooking.DefaultScoreTable(seed_scores=GLib.MappedFile('../data/scores-7mer-3mismatch-ending', False).get_bytes())

        ud = Mirbooking.DefaultScoreTableCutoffFilterUserData()
        ud.broker = broker
        ud.cutoff = 0.1;
        ud.relative_cutoff = 0

        score_table.set_filter(Mirbooking.DefaultScoreTable.cutoff_filter, ud)

        ret, positions = score_table.compute_positions(mirna, target)

        self.assertListEqual(positions, [2460, 4078, 4154, 4233])

class MirbookingMcffScoreTableTestCase(unittest.TestCase):
    def test_mcff_score_table(self):
        from shutil import which
        if which('mcff') is None:
            self.skipTest("'mcff' is not found in path.")
        score_table = Mirbooking.McffScoreTable()
        ret, positions = score_table.compute_positions(mirna, target)
        self.assertTrue(ret)
        self.assertTrue(len(positions) > 0)
        for p in positions:
            self.assertLess(score_table.compute_score(mirna, target, p).score.kr, float('inf'))

class MirbookingBrokerTestCase(unittest.TestCase):
    def test_run(self):
        mirbooking = Mirbooking.Broker(score_table=SimpleScoreTable())
        mirbooking.set_sequence_quantity(target, 5.0)
        mirbooking.set_sequence_quantity(mirna, 5.0)
        self.assertIsInstance(mirbooking.get_score_table(), SimpleScoreTable)

        # solve steady-state
        while True:
            ret, etr = mirbooking.evaluate()
            self.assertTrue(ret)
            ret = mirbooking.step(Mirbooking.BrokerStepMode.SOLVE_STEADY_STATE, 1.0)
            self.assertTrue(ret)
            if etr < 1:
                break

        for target_site in mirbooking.get_target_sites():
            for occupant in target_site.occupants:
                self.assertTrue(mirbooking.get_occupant_quantity(occupant) > 0)

    def test_get_target_sites_as_dataframe(self):
        try:
            import pandas
        except ImportError:
            return unittest.skip('Pandas is not available.')
        mirbooking = Mirbooking.Broker(score_table=SimpleScoreTable())
        mirbooking.set_sequence_quantity(target, 5.0)
        mirbooking.set_sequence_quantity(mirna, 5.0)
        ret, norm = mirbooking.evaluate()
        mirbooking.step(Mirbooking.BrokerStepMode.SOLVE_STEADY_STATE, 1.0)
        df = mirbooking.get_target_sites_as_dataframe()
        self.assertEqual(71, len(df))
        self.assertEqual('NM_000014.4', df.index[0][1])

unittest.main()
