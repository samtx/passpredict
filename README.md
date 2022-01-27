# passpredict

[![Tests](https://github.com/samtx/passpredict/actions/workflows/main.yml/badge.svg)](https://github.com/samtx/passpredict/actions/workflows/main.yml)

Predict upcoming satellite overpasses over a point on Earth.

This library exposes a command-line interface and a backend API to generate overpass predictions.


## Install

The package and command line tool can be installed from PyPI using `pip`.

    pip install passpredict

If you intend to use the package as a command line tool, it is a good idea to install the package with its dependencies using `pipx`.

    pipx install passpredict


## Command Line Usage

The command line output is generated using Rich tables.

## Example

Predict upcoming visible overpasses of the International Space Station. The location is entered using decimal coordinates, with positive values East and North.

- International Space Station (ID 25544)
- Location: 30.2711&deg; N, 97.1234&deg; W
- Visible passes only

```
$ passpredict -lat 30.2711 -lon -97.1234 -s 25544

Satellite ID 25544 ISS (ZARYA) overpasses 
Lat=30.2711°, Lon=-97.1234°, Timezone America/Chicago
Using TLE with epoch 2022-01-23T21:18:30.062880+00:00
1 25544U 98067A   22023.88784795  .00005671  00000-0  10872-3 0  9995
2 25544  51.6443 331.1875 0006753  55.2214  46.3522 15.49604727322805

┏━━━━━━━━━┳━━━━━━━━━━┳━━━━━━┳━━━━━━┳━━━━━━━━━━┳━━━━━━┳━━━━━━┳━━━━━━━━━━┳━━━━━━┳━━━━━━┳━━━━━━━━━┓
┃         ┃  Start   ┃ Sta… ┃ Sta… ┃    Max   ┃ Max  ┃ Max  ┃    End   ┃ End  ┃ End  ┃         ┃
┃  Date   ┃   Time   ┃   El ┃   Az ┃   Time   ┃  El  ┃  Az  ┃   Time   ┃  El  ┃  Az  ┃  Type   ┃
┡━━━━━━━━━╇━━━━━━━━━━╇━━━━━━╇━━━━━━╇━━━━━━━━━━╇━━━━━━╇━━━━━━╇━━━━━━━━━━╇━━━━━━╇━━━━━━╇━━━━━━━━━┩
│ 2/02/22 │ 19:56:17 │  10° │  NNW │ 19:59:25 │  38° │   NE │ 20:02:33 │  10° │  ESE │ visible │
└─────────┴──────────┴──────┴──────┴──────────┴──────┴──────┴──────────┴──────┴──────┴─────────┘

