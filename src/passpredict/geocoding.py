from __future__ import annotations
import typing

import httpx

from .caches import MemoryCache
from .locations import Location

if typing.TYPE_CHECKING:
    from .caches import BaseCache


class NominatimGeocoder:
    url = "https://nominatim.openstreetmap.org/search"

    def __init__(self, cache: BaseCache = None):
        if cache:
            self.cache = cache
        else:
            self.cache = MemoryCache()

    def query(self, q: str) -> Location:
        """
        Get latitude and longitude for search query
        Can't use Nominatim geocoder in production!
        """
        # check cache for geocoding response
        key = f"location:{q}"
        with self.cache as cache:
            res = cache.get(key)
            if res:
                location = Location(
                    name=res['name'],
                    latitude_deg=res['lat'],
                    longitude_deg=res['lon'],
                    elevation_m=res['h'],
                )
            else:
                location = self._query_nominatim(q)
                # cache response for 30 days
                cache.set(key, location.dict(), ttl=86400 * 30)
        return location

    def _query_nominatim(self, q: str):
        params = {
            'q': q,
            'format': 'json',
            'limit': 1,
        }
        try:
            response = httpx.get(self.url, params=params)
        except httpx.RequestError:
            return None
        location = self._serialize_response(response.json()[0])
        return location

    def _serialize_response(self, res: dict) -> Location:
        """
        Serialize response dictionary from Nominatim
        """
        location = Location(
            latitude_deg=float(res['lat']),
            longitude_deg=float(res['lon']),
            elevation_m=0,
            name=res['display_name'],
        )
        return location
