import numpy as np
from pathlib import Path
from os import PathLike
from astropy.table import QTable, Row
from collections.abc import Callable
from typing import Any

class PersistentTable:
    """ A wrapper around QTable serialized to disk.
    """
    def __init__(self, path:PathLike, initializer:Callable[...,QTable]=QTable) -> None:
        """Create the wrapper

        :param path: path to file with the table data.
        :type path: PathLike
        :param initializer: callable with no arguments returning
                            :py:class:`~astropy.table.QTable`,
                            defaults to QTable constructor
        :type initializer: callable, optional
        """
        self.path_ = Path(path)
        self.table_ = None
        self.init_ = initializer
        self.format_ = 'ascii.ecsv'

    def flush(self) -> None:
        """ Flush the table to the disk
        """
        if self.table_ is not None:
            self.table_.write(self.path_, format=self.format_, overwrite=True)

    def get(self) -> QTable:
        """Get underlying table object

        :return: Table read from disk or initialized by the initializer.
        :rtype: :py:class:`~astropy.table.QTable`
        """
        if not self.table_:
            if self.path_.exists():
                self.table_ = QTable.read(self.path_, format=self.format_)
            else:
                self.table_ = self.init_()
                self.flush()
        return self.table_

    @staticmethod
    def init_from_template(template:QTable) -> QTable:
        """Convenience method to create an empty table from template

        :param template: template table to be used for initialization.
        :type template: :py:class:`~astropy.table.QTable`
        :return: Empty table replicating the structure of the template.
        :rtype: :py:class:`~astropy.table.QTable`
        """
        return QTable(template)[[]]

    def append(self, row:Any):
        """Append a row to the table and flush

        :param row: data row to be added
        :type row: dictionary-like which is acceptable by
                   :py:meth:`~astropy.table.QTable.add_row`
        :return: updated table
        :rtype: :py:class:`~astropy.table.QTable`
        """
        self.get().add_row(row)
        self.flush()
        return self.table_

    def row_by_key(self, field:str, key:Any) -> Row | None:
        """Select row by the value of the field

        :param field: the name of the table column
        :type field: str
        :param key: the value, type must match the data type of the column
        :type key: Any
        :return: table row if table[field]==key;
                 None if key does not match any value in 'field';
        :rtype: Row | None
        :raises: KeyError if 'key' matches multiple rows.
        """
        if self.table_ is None:
            return None
        rows = self.get()[self.table_[field] == key]
        if len(rows) > 1:
            raise KeyError(f"{len(rows)} rows found for {field}={key}")
        return None if len(rows) == 0 else rows[0]

    def row_by_keys(self, keys:dict[str, Any]) -> Row | None:
        """Select row matching values of multiple fields

        :param keys: dictionary with field names as keys and values to match
        :type keys: dict[str, Any]
        :return: table row matching all key-value pairs;
                 None if there is no match;
        :rtype: Row | None
        :raises: KeyError if multiple rows match the criteria.
        """
        if self.table_ is None:
            return None
        filter = np.all([self.get()[field] == val for field, val in keys.items()], axis=0)
        rows = self.get()[filter]
        if len(rows) > 1:
            raise KeyError(f"{len(rows)} rows found for {keys}")
        return None if len(rows) == 0 else rows[0]
