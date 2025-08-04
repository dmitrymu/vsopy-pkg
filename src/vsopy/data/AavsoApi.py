import functools
import io
import astropy.units as u
from astropy.utils.data import download_file


def name_for_api(star_name):
    return star_name.replace(' ', '+')

def q_value(q):
    return q.value if hasattr(q, 'value') else q

class AavsoApi:
    """ Encapsulate HTTP APIs to AAVSO database.

        See https://www.aavso.org/apis-aavso-resources

        All download operations by default are cached using
        :py:func:`astropy.utils.data.download_file`.
    """

    def __init__(self, cache_web_content=True) -> None:
        """Create API client for AAVSO database

        :param cache_web_content: specify whether to cache downloaded content, defaults to True
        :type cache_web_content: bool, optional

        """
        self.cache_web_content = cache_web_content
        pass

    def fetch(self, uri: str):
        """Fetch content from the given URI.

        This method downloads the content from the specified URI and returns
        a file-like object. If the download fails, it returns an
        :py:class:`~io.StringIO` object containing an error message.

        :param uri: URI to fetch content from.
        :type uri: str
        :return: file-like object containing the downloaded content or an error message.
        :rtype: file-like
        :raises Exception: If the download fails, an error message is returned
                           as the content of :py:class:`~io.StringIO`
        """
        try:
            return open(download_file(uri, cache=self.cache_web_content))
        except Exception as ex:
            msg = f"Downloading {uri} failed: {ex}"
            return io.StringIO(msg)

    def fetch_content(self, uri: str) -> str:
        """Fetch content from the given URI and return it as a string.

        This method uses the `fetch` method to download the content from the
        specified URI and reads it into a string. If the download fails, it
        returns an error message.

        :param uri: URI to fetch content from.
        :type uri: str
        :return: Content of the URI as a string or an error message.
        :rtype: str
        """
        with self.fetch(uri) as stream:
            return stream.read()

    def _fetch_content(uri) -> str:
        """Decorator to fetch content from a given URI.
        This decorator wraps a method to automatically fetch content from the
        specified URI when the method is called. The URI can be constructed
        using the method's arguments.

        Usage:
            @AavsoApi._fetch_content
            def my_method(self, arg1, arg2):
                return f"http://example.com/api?arg1={arg1}&arg2={arg2}"

        :param uri: callable that returns a URI string.
        :type uri: callable
        :return: A wrapper function that fetches content from the URI
        :rtype: function
        """
        @functools.wraps(uri)
        def fetcher(self, *args, **kwargs):
            return self.fetch_content(uri(self, *args, **kwargs))
        return fetcher

    @_fetch_content
    def get_vsx_votable(self, star_name:str) -> str:
        """Get star information in VOTable format from AAVSO VSX database.

        This queries AAVSO VSX database for a specific star's information
        using HTTP API (https://www.aavso.org/direct-web-query-vsxvsp).
        The result is VOTable XML (https://www.ivoa.net/documents/VOTable/20191021/index.html).

        :param star_name: Name of the star to query.
        :type star_name: str
        :return: XML string containing star information in VOTable format.
        :rtype: str
        """
        return ("http://www.aavso.org/vsx/index.php?view=query.votable"
                f"&ident={name_for_api(star_name)}")

    @_fetch_content
    def get_star_chart(self, star_name: str, fov=60*u.arcmin, maglimit=16*u.mag) -> str:
        """Get comparison stars for a given target star frrom AAVSO VSP tool

        :param star_name: Name of the target star.
        :type star_name: str
        :param fov: Field of view, defaults to 60 :py:class:`~astropy.units.arcmin`
        :type fov: :py:class:`~astropy.units.Quantity`, optional
        :param maglimit: Limit the list of stars by magnitude, defaults to 16 :py:class:`~astropy.units.mag`
        :type maglimit: :py:class:`~astropy.units.Quantity`, optional
        :return: JSON string containing the list of comparison stars.
        :rtype: str
        """
        return ("https://www.aavso.org/apps/vsp/api/chart/?format=json"
                f"&star={name_for_api(star_name)}&fov={q_value(fov)}&maglimit={q_value(maglimit)}")

    @_fetch_content
    def get_std_field_chart(self, ra=0, dec=0, fov=60, maglimit=16) -> str:
        """Get standard field chart from AAVSO VSP tool.

        This queries AAVSO VSX database for photometry data of stars in
        a standard field. The chart is centered at the specified RA and Dec
        coordinates, with a given field of view and magnitude limit.

        :param ra: RA of the field center, deg, defaults to 0
        :type ra: int, optional
        :param dec: Dec  of the field center, deg, defaults to 0
        :type dec: int, optional
        :param fov: Field of view, arcmin, defaults to 60
        :type fov: int, optional
        :param maglimit: Magnitude limit, defaults to 16
        :type maglimit: int, optional
        :return: JSON string containing the photometry data of the standard field.
        :rtype: str
        """
        return ("https://www.aavso.org/apps/vsp/api/chart/?format=json&special=std_field"
                f"&ra={q_value(ra)}&dec={q_value(dec)}"
                f"&fov={q_value(fov)}&maglimit={q_value(maglimit)}")

    @_fetch_content
    def get_chart_by_id(self, chart_id: str) -> str:
        """Get chart by its ID from AAVSO VSP tool.

        AAVSO VSP caches chart data using char ID which is derived from
        the star name, field of view, and magnitude limit. Reetrieving
        by ID may reduce the load on AAVSO servers.

        :param chart_id: Chart ID to retrieve.
        :type chart_id: str
        :return: JSON string containing the chart data.
        :rtype: str
        """
        return (f"https://apps.aavso.org/vsp/api/chart/{chart_id}/?format=json")

    @_fetch_content
    def get_std_fields(self) -> str:
        """Get the list of standard calibration fields from AAVSO VSX database.

        This queries AAVSO VSX database for the list of standard star fields
        in JSON format.

        :return: JSON string containing the list of standard star fields.
        :rtype: str
        """
        return "https://www.aavso.org/vsx/index.php?view=api.std_fields&format=json"
