from abc import ABCMeta, abstractmethod

import openmm as mm

__all__ = ['Restraint']


class Restraint(metaclass=ABCMeta):
    def __init__(self):
        pass

    @abstractmethod
    def apply(self, system: mm.System, *args, **kwargs):
        pass
