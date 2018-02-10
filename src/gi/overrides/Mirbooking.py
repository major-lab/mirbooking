from ..overrides import override
from ..module import get_introspection_module

Mirbooking = get_introspection_module('Mirbooking')

@override
class Broker(Mirbooking.Broker):
    def get_target_sites_as_dataframe(self):
        import pandas as pd
        def get_target_sites():
            for target_site in self.get_target_sites():
                for occupant in target_site.occupants:
                    yield [target_site.target.get_accession(),
                           target_site.target.get_name(),
                           self.get_sequence_quantity(target_site.target),
                           target_site.position,
                           self.get_target_site_vacancy(target_site),
                           occupant.mirna.get_accession(),
                           occupant.mirna.get_name(),
                           self.get_sequence_quantity(occupant.mirna),
                           self.get_score_table().compute_score(occupant.mirna, target_site.target, target_site.position),
                           occupant.quantity]
        columns = ['target_accession',
                   'target_name',
                   'target_quantity',
                   'position',
                   'vacancy',
                   'mirna_accession',
                   'mirna_name',
                   'mirna_quantity',
                   'probability',
                   'occupancy']
        return pd.DataFrame(get_target_sites(), columns=columns).set_index(['target_accession', 'position', 'mirna_accession'])
