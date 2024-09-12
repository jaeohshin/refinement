from typing import Optional, Union
import os
import shutil


class Renamer:
    def __init__(self, filename: Union[os.PathLike[str], str]):
        self.filename = filename
        self.candidate = self.filename + '.bak'

    def rename(self):
        if os.path.exists(self.filename):
            shutil.move(self.filename, self.candidate)
        self.candidate = None


def backup(filename: Union[os.PathLike[str], str]) -> Optional[str]:
    r = Renamer(filename)
    r.rename()
    return r.candidate
