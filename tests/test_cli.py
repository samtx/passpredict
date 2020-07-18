
import numpy as np
from click.testing import CliRunner

import passpredict.cli


def test_cli():
    runner = CliRunner()
    result = runner.invoke(passpredict.cli.main,'-s 25544 -l "austin, texas" -u -5 -d 2')
    assert result.exit_code == 0