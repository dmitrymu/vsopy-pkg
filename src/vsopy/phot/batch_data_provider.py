import numpy as np
from astropy.table import QTable, join, unique

class BatchDataProvider:
    """Common interface to all photometry-lelated data for a single session. 

    The class loads the following tables from the session directory
    - images:
        - image_id 
        - filter 
        - time
        - exposure 
        - airmass 
        - temperature 
        - path
    - batches: 
        - batch_id 
        - time 
        - temperature 
        - airmass 
        - time_range 
        - temperature_range 
        - airmass_range
    - batch_images:
        - batch_id 
        - image_id
    - measured:
        - image_id 
        - auid 
        - M (mag, err)  instrumental magnitude and its uncertainty
        - flux 
        - snr           in dB
        - peak          normalized by ADU maximum, i.e. [0..1]
    - sequence:
        - auid 
        - band 
        - M (mag, err)  catalog magnitude and its uncertainty
    """
    def __init__(self, session):
        self.batches_ = QTable.read(session.batches_file_path)
        self.batch_images_ = QTable.read(session.batch_images_file_path)
        self.images_ = QTable.read(session.images_file_path)
        self.measured_ = QTable.read(session.measured_file_path)
        measured_stars = unique(self.measured_, keys=['auid'])
        sequence = QTable.read(session.sequence_file_path)
        self.sequence_ = sequence[np.isin(sequence['auid'],(measured_stars['auid']))]

    @property
    def target_auid(self):
        return self.sequence_.meta['auid'] if 'auid' in self.sequence_.meta else None

    def instr(self, band):
        return f"instr {band}"

    def targ(self, band):
        return f"targ {band}"

    def sequence_band(self, band):
        """Given a band, return a table with catalog magnitude in that band.
          Columns: 
            - auid, 
            - name of the band (renamed from 'M')
        """
        M = self.sequence_[self.sequence_['band'] == band]['auid', 'M']
        M.rename_column('M', band)
        return M

    def sequence_band_pair(self, band_pair):
        """ Given a pair of band names, return a table with catalog magnitudes in those bands.
          Columns: 
            - auid, 
            - name of band[0] (renamed from 'M'), 
            - name of band[1] (renamed from 'M')
        """
        bandA, bandB = band_pair
        return join(self.sequence_band(bandA), self.sequence_band(bandB), 'auid')

    def check_band_pair(self, band_pair, auid):
        bandA, bandB = band_pair
        s = join(self.sequence_band(bandA), self.sequence_band(bandB), 'auid')
        return s[s['auid'] == auid][bandA, bandB]

    def batch_band(self, band):
        fltr = self.images_[self.images_['filter'] == band]['image_id',]
        image_M = join(fltr, self.measured_['image_id', 'auid', 'M'])
        batch_M = join(image_M, self.batch_images_, 'image_id')['batch_id', 'auid', 'M']
        batch_M.rename_column('M', self.instr(band))
        return batch_M

    def batch_band_pair(self, band_pair):
        bandA, bandB = band_pair
        return join(self.batch_band(bandA), self.batch_band(bandB), ['batch_id', 'auid'])

    def batch_and_sequence_band_pair(self, band_pair):
        return join(self.batch_band_pair(band_pair), self.sequence_band_pair(band_pair), 'auid')

    def batch_comp_star(self, band_pair, auid):
        batches = self.batch_and_sequence_band_pair(band_pair)
        result = batches[batches['auid'] == auid]['batch_id', band_pair[0], band_pair[1],
                                                 self.instr(band_pair[0]), self.instr(band_pair[1])]
        result.sort('batch_id')
        return result

    def batch_target_star(self, band_pair, auid):
        batches = self.batch_band_pair(band_pair)
        result = batches[batches['auid'] == auid]['batch_id', self.instr(band_pair[0]), self.instr(band_pair[1])]
        result.rename_columns([self.instr(band_pair[0]), self.instr(band_pair[1])],
                              [self.targ(band_pair[0]), self.targ(band_pair[1])])
        result.sort('batch_id')
        return result