from importlib.resources import path
import sys
import pathlib

from examples.celestrak_source import celestrak_source


# Set example directory on python path
example_dir = pathlib.Path(__file__).parent.parent / 'examples'
sys.path.append(example_dir)

from examples import *


def test_brute_force_observer():
    """
    Make sure example doesn't error
    """
    brute_force_observer()


def test_vallado_predict_11_6():
    """
    Make sure example doesn't error
    """
    vallado_predict_11_6()


def test_standard_observer():
    """
    Make sure example doesn't error
    """
    standard_observer()


def test_celestrak_source():
    """
    Make sure example doesn't error
    """
    celestrak_source()