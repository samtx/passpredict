from __future__ import annotations
import typing

import httpx

from passpredict.caches import JsonCache

from .locations import Location


class NominatimGeocoder:
    url = "https://nominatim.openstreetmap.org/search"
    cache = JsonCache('locations.json')

    @classmethod
    def query(cls, q: str) -> Location:
        """
        Get latitude and longitude for search query
        Can't use Nominatim geocoder in production!
        """
        # check cache for geocoding response
        key = f"location:{q}"
        with cls.cache as cache:
            res = cache.get(key)
            if res:
                location = Location(
                    name=res['name'],
                    latitude_deg=res['lat'],
                    longitude_deg=res['lon'],
                    elevation_m=res['h'],
                )
            else:
                location = cls._query_nominatim(q)
                # cache response for 30 days
                cache.set(key, location.dict(), ttl=86400 * 30) 
        return location

    @classmethod
    def _query_nominatim(cls, q: str):
        params = {
            'q': q,
            'format': 'json',
            'limit': 1,
        }
        try:
            response = httpx.get(cls.url, params=params)
        except httpx.RequestError:
            return None
        location = cls._serialize_response(response.json()[0])
        return location

    def _serialize_response(res: dict) -> Location:
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
