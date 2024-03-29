""" BioSimulators-compliant command-line interface to the `MySimulator <https://url.for.my.simulator>`_ simulation program.

:Author: Author name <email@organization>
:Date: YYYY-MM-DD
:Copyright: YYYY, Owner
:License: <License, e.g., MIT>
"""

from ._version import __version__
from .core import get_simulator_version, exec_sedml_docs_in_combine_archive
from biosimulators_utils.simulator.cli import build_cli

App = build_cli('biosimulators-mpbn', __version__,
                'mpbn', get_simulator_version(), 'https://mpbn.readthedocs.io/',
                exec_sedml_docs_in_combine_archive)


def main():
    with App() as app:
        app.run()
