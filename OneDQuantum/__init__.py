#!/usr/bin/env python
# -*- coding:utf-8 -*-
from . import OneDSchrodinger, OneDThermal, OneDMaxwell
from .OneDSchrodinger import *
from .OneDThermal import *
from .OneDMaxwell import *
__all__ = (OneDSchrodinger.__all__ 
           + OneDThermal.__all__ 
           + OneDMaxwell.__all__  
           + ['OneDSchrodinger', 'OneDThermal', 'OneDMaxwell'])

try: 
    OneDSchrodinger.bindOpenMP(True)
except OSError:
    OneDSchrodinger.bindOpenMP(False)

# vim: ts=4 sw=4 sts=4 expandtab
