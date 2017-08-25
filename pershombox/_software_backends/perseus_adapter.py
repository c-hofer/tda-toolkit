from subprocess import call
from tempfile import TemporaryDirectory
from subprocess import DEVNULL
import os
import numpy
import warnings


__ext_lib_path = os.path.join(os.path.dirname(__file__), 'ext_lib')

__perseus_path = os.path.join(__ext_lib_path, 'perseusLin')

if not os.path.exists(__perseus_path):

    text  = """Found no perseusLin executable in  {}. Get it from
            http://people.maths.ox.ac.uk/nanda/perseus/index.html and copy
            it into {} or set the perseus_path variable
            of this module after import manually.
            """.format(__perseus_path, __ext_lib_path)

    warnings.warn(text, ImportWarning)


def _call_perseus(complex_type, complex_file_string):
    def get_dim_from_dgm_file(name):
        x = name.split('.txt')[0]
        x = x.split('_')[1]
        return int(x)

    with TemporaryDirectory() as tmp_dir:
        comp_file_path = os.path.join(tmp_dir, 'complex.txt')
        with open(comp_file_path, 'w') as comp_file:
            comp_file.write(complex_file_string)

        call([__perseus_path, complex_type, comp_file_path, tmp_dir + '/'], stdout=DEVNULL)

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
