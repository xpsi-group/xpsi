import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../xpsi")))

from xpsi import __version__

print(".".join(__version__.split(".")[:2]))
