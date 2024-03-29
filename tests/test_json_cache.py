import json
import time

import pytest

from passpredict import caches

@pytest.fixture(scope='function')
def tle_json_cache():
    data = {
        'sat:12345': {
            'ttl': time.time() + 86400,
            'data': {
                'an aribtrary key': 'value',
                'ttl': 4566987,
            }
        },
        'sat:48852': {
            'ttl': time.time() + 5,
            'data': {
                'tle1': '1 48852U 21053A   21197.14114025  .00007723  00000-0  87850-4 0  9997',
                'tle2': '2 48852  41.4696 106.4609 0006818 163.3920 306.0561 15.62530664 12228',
                'name': 'SHENZHOU-12',
            }
        }

    }
    return data

@pytest.fixture(scope="function")
def cache(tmp_path, tle_json_cache):
    fname = tmp_path / 'tle.json'
    with open(fname, 'w') as f:
        json.dump(tle_json_cache, f)
    with caches.JsonCache(fname) as cache:
        yield cache


def test_jsoncache_opens(tmp_path, tle_json_cache):
    fname = tmp_path / 'tle.json'
    with open(fname, 'w') as f:
        json.dump(tle_json_cache, f)
    with caches.JsonCache(fname) as cache:
        res = cache.get('sat:12345')
        assert res['an aribtrary key'] == 'value'
        assert res['ttl'] == 4566987


def test_jsoncache_filenotfound(tmp_path):
    """ If the file doesn't exist, create empty cache """
    fname = tmp_path / 'fake.json'
    cache = caches.JsonCache(fname)
    cache.load()
    assert not cache._cache
    assert len(cache._cache.items()) == 0


def test_set_jsoncache_ttl(cache):
    key = 'sat:123'
    item = {
        'tle1': 'a long tle',
        'tle2': 'another long tle'
    }
    cache.set(key, item, ttl=-10)
    res = cache.get(key)
    assert res is None
    assert key not in cache


def test_get_jsoncache_item(cache):
    key = 'sat:12345'
    res = cache.get(key)
    assert res
    assert key in cache


def test_pop_jsoncache_item(cache):
    key = 'sat:12345'
    assert key in cache
    res = cache.pop(key)
    assert res
    assert key not in cache


def test_pop_default_jsoncache_item(cache):
    key = 'sat:12345'
    assert key in cache
    res = cache.pop(key)
    assert res
    assert key not in cache
    res = cache.pop(key)
    assert res is None


def test_set_jsoncache_ttl_to_none(cache):
    key = 'sat:123'
    item = {
        'tle1': 'a long tle',
        'tle2': 'another long tle'
    }
    cache.set(key, item, ttl=None)
    res = cache.get(key)
    assert res['tle1'] == 'a long tle'
    assert res['tle2'] == 'another long tle'
    assert key in cache
