from ErwinJr2.OneDQuantum import OneDSchrodinger

__all__ = (OneDSchrodinger.__all__
           + ['OneDSchrodinger'])

try:
    OneDSchrodinger.bindOpenMP(True)
except OSError:
    OneDSchrodinger.bindOpenMP(False)
