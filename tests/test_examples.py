from importlib.resources import path
import sys
import pathlib


# Set example directory on python path
example_dir = pathlib.Path(__file__).parent.parent / 'examples'
sys.path.append(example_dir)

from examples import *

"""
Make sure examples don't error
"""

def test_brute_force_observer():
    brute_force_observer()


def test_vallado_predict_11_6():
    vallado_predict_11_6()


def test_standard_observer():
    standard_observer()


def test_celestrak_source():
    celestrak_source()

def test_all_visual_satellites():
    all_visual_satellites()