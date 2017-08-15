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
                    yield target_site.target.get_accession(), occupant.mirna.get_accession(), target_site.position, self.get_score_table().compute_score(occupant.mirna, target_site.target, target_site.position), occupant.quantity, self.get_target_site_vacancy(target_site)
        return pd.DataFrame(get_target_sites(), columns=['target', 'mirna', 'position', 'probability', 'occupancy', 'vacancy']).set_index(['target', 'mirna', 'position'])
