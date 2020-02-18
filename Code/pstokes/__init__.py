# -*- coding: utf-8 -*-
"""
pstokes
=======

A small module using FEniCS to solve the p-Stokes equations on
periodic domains. In particular, this could be suitable for
glacier-flow modelling on ISMIP-HOM-type domains.

"""

# star import everything
from .pstokes import *
from .domain_utils import *
from .helpers import *

__package__ = 'pstokes'
__title__ = 'p-Stokes'
__version__ = '0.1'
__author__ = 'Christian Helanow'
__license__ = 'GPL2'
__copyright__ = 'Copyright 2020 Christian Helanow'
__docformat__ = 'restructuredtext'
