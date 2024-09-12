import contextlib, os, tempfile

__all__ = ['change_working_directory', 'use_temporary_directory']


@contextlib.contextmanager
def change_working_directory(path):
    cwd = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(cwd)


@contextlib.contextmanager
def use_temporary_directory():
    tmp_dir = tempfile.TemporaryDirectory()
    with change_working_directory(tmp_dir.name):
        try:
            yield
        finally:
            pass
