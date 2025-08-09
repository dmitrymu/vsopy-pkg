from os import PathLike
from pathlib import Path
from typing_extensions import deprecated

class LayoutBase:
    """ Generic layout
        Supports root directory and enforces directory creation.
    """
    def __init__(self, root:PathLike, create:bool=True):
        """ Set up layout at root directory.

            :param root: path to root directory.
            :type root: path-like
            :param create: whether to enforce directory creation. Defaults to True.
            :type create: bool, optional
        """
        self.create_ = create
        self.root_ = Path(root)

    def _enforce(dir):
        """ Decorator enforcing directory creation.

            dir - a method creating directory
        """
        def enforcer(self, *args, **kwargs):
            result = Path(dir(self, *args, **kwargs))
            if self.create_ and not result.exists():
                result.mkdir(parents=True)
            return result
        return enforcer

    @property
    @_enforce
    def root_dir(self) -> Path:
        """ Returns path to the root directory
        """
        return self.root_


class Session:
    """ Session object representing a measurement session.
    """
    def __init__(self, tag:str=None, name:str=None):
        """Create a session object.

        :param tag: session tag, defaults to None
        :type tag: str, optional
        :param name: session object name, defaults to None
        :type name: str, optional
        """
        self.tag_ = tag
        self.name_ = name.replace(' ', '_')

    @property
    def rel_path(self) -> Path:
        """ Returns relative path to the session directory

        :return: relative path to the session directory
        :rtype: Path
        """
        return Path(self.tag_) / Path(self.name_)

    @property
    def name(self):
        return self.name_.replace('_', '')


class TargetLayout(LayoutBase):
    """ Layout for target images in a session.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @property
    @LayoutBase._enforce
    def lights_dir(self) -> Path:
        return self.root_dir / 'Light'


class ImageLayout(LayoutBase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_images(self, session:Session) -> TargetLayout:
        return TargetLayout(self.root_dir / session.rel_path,
                           create=self.create_)


class ChartsLayout(LayoutBase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @property
    def charts_file_path(self):
        return self.root_dir / 'charts.ecsv'

    @property
    def std_fields_file_path(self):
        return self.root_dir / 'std_fields.ecsv'

    @property
    def targets_file_path(self):
        return self.root_dir / 'targets.ecsv'

    def get_centroid_file_path(self, name):
        return self.root_dir / f'{name}_c.ecsv'

    def get_sequence_file_path(self, name):
        return self.root_dir / f'{name}_s.ecsv'


class SessionLayout(LayoutBase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @property
    @LayoutBase._enforce
    def solved_dir(self):
        return self.root_dir / 'solved'

    @property
    def blacklist_file_path(self):
        return self.root_dir / 'blacklist.json'

    @property
    def batches_file_path(self):
        return self.root_dir / 'batches.ecsv'

    @property
    def batch_images_file_path(self):
        return self.root_dir / 'batch_images.ecsv'

    @property
    @deprecated("Use centroid_file_path and sequence_file_path instead")
    def chart_file_path(self):
        return self.root_dir / 'chart.ecsv'

    @property
    def centroid_file_path(self):
        return self.root_dir / 'centroids.ecsv'

    @property
    def sequence_file_path(self):
        return self.root_dir / 'sequence.ecsv'

    @property
    def images_file_path(self):
        return self.root_dir / 'images.ecsv'

    @property
    def settings_file_path(self):
        return self.root_dir / 'settings.json'

    @property
    def measured_file_path(self):
        return self.root_dir / 'measured.ecsv'

    @property
    @deprecated("Use measured_file_path instead")
    def photometry_file_path(self):
        return self.root_dir / 'photometry.ecsv'


class WorkLayout(LayoutBase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @property
    @LayoutBase._enforce
    def tmp_dir(self):
        return self.root_dir / 'tmp'

    @property
    @LayoutBase._enforce
    @deprecated("Use charts instead")
    def charts_dir(self):
        return self.root_dir / 'charts'

    @property
    def charts(self):
        return ChartsLayout(self.root_dir / 'charts', create=self.create_)

    @property
    @LayoutBase._enforce
    def calibr_dir(self):
        return self.root_dir / 'calibr'

    def get_session(self, session) -> SessionLayout:
        return SessionLayout(self.root_dir / Path('session') / session.rel_path,
                             create=self.create_)

