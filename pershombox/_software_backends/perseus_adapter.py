import os
import numpy
from subprocess import call
from tempfile import TemporaryDirectory
from subprocess import DEVNULL
from .resource_handler import get_path, Backends


__stdout = DEVNULL
__stderr = DEVNULL


def _get_perseus_path():
    return get_path(Backends.perseus)


def _call_perseus(complex_type, complex_file_string):
    def get_dim_from_dgm_file(name):
        x = name.split('.txt')[0]
        x = x.split('_')[1]
        return int(x)

    with TemporaryDirectory() as tmp_dir:
        comp_file_path = os.path.join(tmp_dir, 'complex.txt')
        perseus_path = _get_perseus_path()

        with open(comp_file_path, 'w') as comp_file:
            comp_file.write(complex_file_string)

        call([perseus_path, complex_type, comp_file_path, tmp_dir + '/'], stdout=__stdout, stderr=__stderr)

        # dgm file names are assumed to be like output_0.txt
        diagram_files = [name for name in os.listdir(tmp_dir) if name.startswith('_') and name != '_betti.txt']

        dgms = {}
        for name in diagram_files:
            dim = get_dim_from_dgm_file(name)
            dgm_file_path = os.path.join(tmp_dir, name)

            if os.stat(dgm_file_path).st_size == 0:
                dgms[dim] = []

            else:
                dgm = numpy.loadtxt(dgm_file_path)

                if dgm.ndim == 2:
                    dgms[dim] = dgm.tolist()
                elif dgm.ndim == 1:
                    dgms[dim] = [dgm.tolist()]
                else:
                    raise ValueError('Oddly shaped array read from dgm_file_path.')

        return dgms


class PerseusAdapterException(Exception):
    pass
