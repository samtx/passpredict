import requests

nominatim_url = "https://nominatim.openstreetmap.org/search"

def geocoder(q, source=nominatim_url):
    """
    Get latitude and longitude for search query

    Can't use Nominatim geocoder in production!
    """
    params = {
        'q': q,
        'format': 'json',
    }
    response = requests.get(source, params=params)
    data = response.json()
    return data[0]