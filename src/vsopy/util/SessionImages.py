import astropy.units as u
import ccdproc as ccdp
import numpy as np

from astropy.io.typing import PathLike
from astropy.table import QTable
from astropy.time import Time
from vsopy.util import TargetLayout


def session_image_list(image_layout:TargetLayout) -> QTable:
    """Traverse image folder structure and return a table of images.

    :param image_layout: Layout object that provides lights_dir property
                         for accessing the image directories.
    :type image_layout: :py:class:`vsopy.util.TargetLayout`
    :return: Table containing image metadata.
    :rtype: :py:class:`~astropy.table.QTable` with columns:
        - image_id: Unique identifier for the image (sequential int).
        - filter: Filter used for the image.
        - time: Observation time of the image.
        - exposure: Exposure time of the image in seconds.
        - airmass: Air mass at the time of observation.
        - temperature: CCD temperature during the exposure.
        - path: Full path to the image file.
    """
    KEYS = set(['file', 'exptime', 'ccd-temp', 'filter', 'airmass', 'date-obs'])

    files = QTable([dict({key: row[key]
                          for key in row.colnames
                          if key in KEYS},
                         dir=d)
                    for d in image_layout.lights_dir.iterdir()
                    if d.is_dir()
                    for row in ccdp.ImageFileCollection(d).summary])
    files['image_id'] = [n+1 for n in range(len(files))]
    files['time'] = Time(files['date-obs'])
    files['exposure'] = files['exptime'] * u.second
    files['temperature'] = files['ccd-temp'] * u.deg_C
    files.add_column([str(row['dir'] / row['file']) for row in files], name='path')

    images = files['image_id', 'filter', 'time', 'exposure', 'airmass', 'temperature', 'path']
    images.meta = {'start': np.min(images['time']), 'finish': np.max(images['time'])}

    return images

def batch_session_images(image_list_path:PathLike) -> tuple[QTable, QTable]:
    """Create batches of images based on filter and time.

    This function assume that images are taken with N diffrent filters and
    filter are used in the same order repeatedly.  E.g. if there are 3 filters
    V, R, I, then images are taken in the order V, R, I, V, R, I, etc.
    In this case each batch (V, R, I) may be considered as a single observation
    for the purpose of diffrential photometry.

    This function groups sequential images into batches.  It is assumed
    that the first batch contains images in the correct order of filters.
    During further processing, if batch of N images has incorrect order,
    first is skipped and the next N images are tried.

    :param image_list_path: path to the serialized image list table created
                            by :py:func:`session_image_list`.
    :type image_list_path: PathLike
    :return: tuple of two QTables:

        - batches: table with batch metadata, fields:
            - batch_id: unique identifier for the batch
            - time: mean observation time for the batch
            - time_range: total duration of the batch
            - temperature: mean CCD temperature for the batch
            - temperature_range: temperature range of the batch
            - airmass: mean airmass for the batch
            - airmass_range: airmass range of the batch
        - batch_images: table with image IDs in each batch, fields:
            - batch_id: identifier of the batch
            - image_id: identifier of the image in the batch

    :rtype: tuple[QTable, QTable]
    """
    images = QTable.read(image_list_path)
    images.sort('time')
    bands = set(images['filter'])
    batch_size = len(bands)
    order = list(images[0:batch_size]['filter'])

    images.add_column(0, name='batch_id')
    images.add_column(0*u.second, name='time_range')
    images.add_column(0*images['temperature'].unit, name='temperature_range')
    images.add_column(0., name='airmass_range')
    batches = images['batch_id', 'time', 'temperature', 'airmass',
                     'time_range', 'temperature_range', 'airmass_range'][:0].copy()
    batch_images = images['batch_id', 'image_id'][:0].copy()

    next = 0
    id = 0
    while next < len(images) - batch_size + 1:
        batch = images[next:next+batch_size]
        if np.all(batch['filter'] == order):
            id += 1
            batches.add_row({
                'batch_id': id,
                'time': np.mean(batch['time']),
                'time_range': np.max(batch['time']) - np.min(batch['time']),
                'temperature': np.mean(batch['temperature']),
                'temperature_range': np.max(batch['temperature']) - np.min(batch['temperature']),
                'airmass': np.mean(batch['airmass']),
                'airmass_range': np.max(batch['airmass']) - np.min(batch['airmass'])
            })
            for image in batch:
               batch_images.add_row(dict(
                   batch_id=id,
                   image_id=image['image_id']
               ))
            next += batch_size
        else:
            print(f"Skipped file #{next+1} '{images['file']}'")
            next += 1

    batches.meta = {'start': np.min(images['time']), 'finish': np.max(images['time']),
                    'bands': [str(x) for x in order]}

    return batches, batch_images
