#!/usr/bin/env python
# -*- coding:utf-8 -*-
from PyQt5.QtWidgets import QWidget, QComboBox, QFrame
class mtrlComboBox(QComboBox):
    """A QComboBox that pops up on the right"""
    def __init__(self, parent=None):
        super(mtrlComboBox, self).__init__(parent)

    def showPopup(self):
        super(mtrlComboBox, self).showPopup()
        popup = self.findChild(QFrame)
        popup.move(popup.x()+self.size().width(), popup.y())

# vim: ts=4 sw=4 sts=4 expandtab
