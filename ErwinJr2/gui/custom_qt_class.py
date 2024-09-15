"""This file defined customized Qt widgets for ErwinJr2"""

from PyQt5.QtWidgets import QComboBox, QFrame


class MtrlComboBox(QComboBox):
    """A QComboBox that pops up on the right"""

    def showPopup(self):  # override QComboBox.showPopup pylint: disable=invalid-name
        super().showPopup()
        popup = self.findChild(QFrame)
        popup.move(popup.x() + self.size().width(), popup.y())
