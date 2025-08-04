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
        """Create API client for AAVSO database.

        Args:
            cache_web_content (bool, optional): specify whether to cache
            downloaded content. Defaults to True.
        """
        self.cache_web_content = cache_web_content
        pass

    def fetch(self, uri: str):
        """Fetch content from the given URI.

        This method downloads the content from the specified URI and returns
        a file-like object. If the download fails, it returns an
        :py:class:`io.StringIO` object containing an error message.

        Args:
            uri (str): The URI to fetch content from.

        Returns:
            file-like: A file-like object containing the downloaded content or an error message.

        Raises:
            Exception: If the download fails, an error message is returned
            as the content of :py:class:`io.StringIO` object.
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

        Args:
            uri (str): The URI to fetch content from.

        Returns:
            _type_: _description_
        """
        with self.fetch(uri) as stream:
            return stream.read()

    def _fetch_content(uri) -> str:
        """Decorator to fetch content from a given URI.

        This decorator wraps a method to automatically fetch content from the
        specified URI when the method is called. The URI can be constructed
        using the method's arguments.

        Args:
            uri (str): The URI to fetch content from, which is a method that
                       takes arguments and returns a formatted URI.

        Returns:
            function: A wrapper function that fetches content from the URI
                      constructed by calling the provided method with the given arguments.

        Usage:
            @AavsoApi._fetch_content
            def my_method(self, arg1, arg2):
                return f"http://example.com/api?arg1={arg1}&arg2={arg2}"

        The decorated method should return a URI string that can be used to
        fetch content. The decorator will handle the fetching and return the
        content as a string.
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

        Args:
            star_name (str): Name of the star to query.

        Returns:
            str: VOTable XML string containing star information.
        """
        return ("http://www.aavso.org/vsx/index.php?view=query.votable"
                f"&ident={name_for_api(star_name)}")

    @_fetch_content
    def get_star_chart(self, star_name: str, fov=60*u.arcmin, maglimit=16*u.mag) -> str:
        """Get comparison stars for a given target star frrom AAVSO VSP tool.

        Args:
            star_name (str): Name of the target star.
            fov (Quantity, optional): Field of view. Defaults to 60*u.arcmin.
            maglimit (Quantity, optional): Limit the list of stars by magnitude. Defaults to 16*u.mag.

        Returns:
            str: JSON string containing the list of comparison stars.
        """
        return ("https://www.aavso.org/apps/vsp/api/chart/?format=json"
                f"&star={name_for_api(star_name)}&fov={q_value(fov)}&maglimit={q_value(maglimit)}")

    @_fetch_content
    def get_std_field_chart(self, ra=0, dec=0, fov=60, maglimit=16) -> str:
        """Get standard field chart from AAVSO VSP tool.

        This queries AAVSO VSX database for photometry data of stars in
        a standard field. The chart is centered at the specified RA and Dec
        coordinates, with a given field of view and magnitude limit.

        Args:
            ra (float, optional): RA of the field center. Defaults to 0.
            dec (float, optional): Dec  of the field center. Defaults to 0.
            fov (float, optional): Field of view Defaults to 60.
            maglimit (float, optional): Magnitude limit. Defaults to 16.

        Returns:
            str: Photometry data in JSON format.
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

        Args:
            chart_id (str): Chart ID to retrieve.

        Returns:
            str: JSON string containing the chart data.
        """
        return (f"https://apps.aavso.org/vsp/api/chart/{chart_id}/?format=json")

    @_fetch_content
    def get_std_fields(self) -> str:
        """Get the list of standard calibration fields from AAVSO VSX database.

        This queries AAVSO VSX database for the list of standard star fields.
        in JSON format.

        Returns:
            str: JSON string containing the list of standard star fields.
        """
        return "https://www.aavso.org/vsx/index.php?view=api.std_fields&format=json"
