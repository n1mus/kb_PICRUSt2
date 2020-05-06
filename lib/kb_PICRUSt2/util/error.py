"""
Exception library to reduce hardcoding and categorize exceptions
Useful for testing
Minimize namespace so can import *
"""


class NonZeroReturnException(Exception):
    pass


class NoWorkspaceReferenceException(Exception):
    pass


class ArgumentException(Exception):
    pass








