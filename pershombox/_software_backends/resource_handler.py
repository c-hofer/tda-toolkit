import os
import warnings
from configparser import ConfigParser
from subprocess import call, DEVNULL
from enum import Enum


__CFG_FILE_NAME = 'software_backends.cfg'


__cfg_path = os.path.join(os.path.dirname(__file__), __CFG_FILE_NAME)


class Backends(Enum):
    dipha = 'dipha'
    perseus = 'perseus'
    hera_wasserstein_dist = 'hera_wasserstein_dist'


__fall_backs = {
    Backends.dipha: 'dipha',
    Backends.hera_wasserstein_dist: 'hera',
    Backends.perseus: 'perseus'
}


__paths_or_errors = {
   Backends.dipha: None,
   Backends.hera_wasserstein_dist: None,
   Backends.perseus: None
}


parser = ConfigParser()
parser.read(__cfg_path)


class SoftwareBackendError(Exception):
    pass


def init_backend(backend: Backends):
    path = parser.get('paths', backend.value)

    if path == '':
        path = __fall_backs[backend]

    try:
        call([path], stdout=DEVNULL, stderr=DEVNULL)

    except Exception as ex:
        __paths_or_errors[backend] = ex

    else:
        __paths_or_errors[backend] = path


def get_path(backend: Backends)->str:
    path_or_error = __paths_or_errors[backend]

    if isinstance(path_or_error, Exception):
        ex_text = "{} backend software is not available.".format(backend.value)
        new_ex = SoftwareBackendError(ex_text)
        raise new_ex from path_or_error

    else:
        return path_or_error


def get_backend_cfg_errors():
    return [(b.value, e) for b, e in __paths_or_errors.items() if isinstance(e, Exception)]


def init_software_backends():

    for software_backend in Backends:
        init_backend(software_backend)

    backend_errors = get_backend_cfg_errors()
    if len(backend_errors) > 0:

        error_text = ""
        error_text += "The following backends are not properly configured\n"

        for b, _ in backend_errors:
            error_text += (b) + '\n'

        error_text += "Using stuff dependent on those backends will cause runtime errors.\n"
        error_text += "You can get all errors by calling pershombox.get_backend_cfg_errors().\n"

        warnings.warn(error_text, UserWarning)