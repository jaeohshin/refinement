import subprocess
from typing import AnyStr, Callable, List

__all__ = ['run']


def identity(x: AnyStr) -> str:
    return str(x)


def run(args: List[AnyStr],
        parser: Callable[[AnyStr], AnyStr] = identity,
        input_text: AnyStr = '') -> str:
    result = subprocess.run(args, input=input_text,
                            stdout=subprocess.PIPE,
                            text=True)
    return parser(result.stdout)
