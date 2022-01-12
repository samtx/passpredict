import httpx

nominatim_url = "https://nominatim.openstreetmap.org/search"

def geocoder(q, source=nominatim_url):
    """
    Get latitude and longitude for search query

    Can't use Nominatim geocoder in production!
    """
    params = {
        'q': q,
        'format': 'json',
        'limit': 1,
    }
    try:
        response = httpx.get(source, params=params)
        data = response.json()[0]
    except httpx.RequestError:
        data = None
    return data