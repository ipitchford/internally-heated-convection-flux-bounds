"""Project paths, relative to this file, so that the bundle runs from wherever it is unpacked."""
import os
CODE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(CODE)
DNS = os.path.join(ROOT, 'data', 'dns')
STEADY = os.path.join(ROOT, 'data', 'steady')
PAPER = os.path.join(ROOT, 'paper')
FIGS = os.path.join(ROOT, 'figs')
LOGS = os.path.join(ROOT, 'logs')
