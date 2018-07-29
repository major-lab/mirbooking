import gi
gi.require_version('Mirbooking', '1.0')
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

class SimpleScoreTable(Mirbooking.ScoreTable):
    def do_compute_score(self, mirna, target, position):
        # simple hamming distance
        if position > target.get_sequence_length() - 7:
            return float('inf')

        d = sum(1 if a == b else 0
                for a, b in zip(reverse_complement(mirna.get_subsequence (1, 7)), target.get_subsequence (position, 7)))

        # at most 1 mismatch
        if d == 6:
            return -d
        else:
            return float('inf')

class MirbookingScoreTableTestCase(unittest.TestCase):
    def test_simple_score_table(self):
        score_table = SimpleScoreTable()
        scores, positions = score_table.compute_scores(mirna, target)
        self.assertEqual(len(scores), len(positions))
        for i, s in enumerate(scores):
            self.assertEqual(s, score_table.compute_score(mirna, target, positions[i]))

class MirbookingBrokerTestCase(unittest.TestCase):
    pass

def t():
    def test_run(self):
        mirbooking = Mirbooking.Broker(score_table=SimpleScoreTable())
        mirbooking.set_sequence_quantity(target, 5.0)
        mirbooking.set_sequence_quantity(mirna, 5.0)
        norm = 0
        mirbooking.evaluate(0, norm)
        mirbooking.step(0, norm)
        for target_site in mirbooking.get_target_sites():
            for occupant in target_site.occupants:
                print('Target {} receives {} miRNA {} at position {}.'.format(target_site.target.get_accession(), occupant.quantity, occupant.mirna.get_accession(), target_site.position))

    def test_get_target_sites_as_dataframe(self):
        mirbooking = Mirbooking.Broker(score_table=SimpleScoreTable())
        mirbooking.set_sequence_quantity(target, 5.0)
        mirbooking.set_sequence_quantity(mirna, 5.0)
        mirbooking.evaluate(0, norm=None)
        mirbooking.step(0, step_size=1)
        df = mirbooking.get_target_sites_as_dataframe()
        self.assertEqual(5, len(df))
        self.assertEqual('NM_000014.4', df.index[0][0])

unittest.main()
