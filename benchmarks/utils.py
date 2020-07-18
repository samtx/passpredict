# benchmark utilities

import pickle
from pathlib import Path

DIRPATH = Path(__file__).parent / 'data'
    

def load_pickle(name):
    with open(DIRPATH / '{}.pkl'.format(name), 'rb') as f:
        data = pickle.load(f)
    return data


def save_pickle(data, name):
    with open(DIRPATH / '{}.pkl'.format(name), 'wb') as f:
        pickle.dump(data, f)