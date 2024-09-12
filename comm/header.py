from typing import List

__all__ = ['Header']


class Header:
    """
    PDB header for CASP15 submission.
    """

    def __init__(self, pfrmat, target, author, remark, method):
        """
        PFRMAT: str = TS, QA, LG
        TARGET: str = Txxx
        AUTHOR: str = XXXX-XXXX-XXXX
        METHOD: str =  Description of methods used
        """
        kwargs = dict(PFRMAT=pfrmat, TARGET=target, AUTHOR=author, REMARK=remark, METHOD=method)
        self.__dict__.update(kwargs)

    def get_header(self) -> List[str]:
        header = ['PFRMAT {}'.format(self.PFRMAT),
                  'TARGET {}'.format(self.TARGET),
                  'AUTHOR {}'.format(self.AUTHOR),
                  'REMARK {}'.format(self.REMARK),
                  'METHOD {}'.format(self.METHOD)]
        return header
