from ..overrides import override
from ..module import get_introspection_module

Mirbooking = get_introspection_module('Mirbooking')

@override
class Sequence(Mirbooking.Sequence):
    def get_subsequence(self, subsequence_offset, subsequence_len):
        return Mirbooking.Sequence.get_subsequence(self, subsequence_offset, subsequence_len)[:subsequence_len]

@override
class Broker(Mirbooking.Broker):
    def get_target_sites_as_dataframe(self):
        import pandas as pd
        def get_target_sites():
            cur_target = None
            for target_site in self.get_target_sites():
                if target_site.target != cur_target:
                    cur_target = target_site.target
                    target_quantity = self.get_sequence_quantity(target_site.target)
                for occupant in target_site.occupants:
                    yield [target_site.target.get_accession(),
                           target_site.target.get_name(),
                           target_quantity,
                           target_site.position + 1,
                           occupant.mirna.get_accession(),
                           occupant.mirna.get_name(),
                           self.get_sequence_quantity(occupant.mirna) + self.get_bound_mirna_quantity(occupant.mirna),
                           (occupant.score.kr + occupant.score.kcat + self.get_target_site_kother(target_site)) / occupant.score.kf,
                           self.get_occupant_quantity(occupant)]
        columns = ['target_accession',
                   'target_name',
                   'target_quantity',
                   'position',
                   'mirna_accession',
                   'mirna_name',
                   'mirna_quantity',
                   'score',
                   'quantity']
        return pd.DataFrame(get_target_sites(), columns=columns).set_index(['target_accession', 'position', 'mirna_accession'])
