import sys
import os


def check_pershombox_availability():
    try:
        import pershombox

    except ImportError:
        sys.path.append(os.path.dirname(os.getcwd()))

        try:
            import pershombox

        except ImportError as ex:
            raise ImportError(
"""
Could not import pershombox. Running your python interpreter in the 'tutorials' sub folder could resolve this issue. 
"""
            ) from ex

