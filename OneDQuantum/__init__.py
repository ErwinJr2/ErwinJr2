#!/usr/bin/env python
# -*- coding:utf-8 -*-
from . import OneDSchrodinger
from .OneDSchrodinger import *
__all__ = (OneDSchrodinger.__all__
           + ['OneDSchrodinger'])

try:
    OneDSchrodinger.bindOpenMP(True)
except OSError:
    OneDSchrodinger.bindOpenMP(False)

# vim: ts=4 sw=4 sts=4 expandtab
