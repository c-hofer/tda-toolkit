import os
import importlib.util


ext_lib_path = os.path.join(os.path.dirname(__file__), 'ext_lib')

try:
    spec = importlib.util.spec_from_file_location("_dionysus", os.path.join(ext_lib_path, '_dionysus.so'))
    _dionysus = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(_dionysus)

except ImportError as ex:
    import os
    new_ex = ImportError("It seems that there is no _dionysus.so in {}.".format(os.path.dirname(__file__)) +
                         "You can find instructions to solve this issue in this file: {}.".format(__file__))
    raise new_ex from ex

"""
Dionysus installation:
Documentation: http://www.mrzv.org/software/dionysus/
Repository: http://hg.mrzv.org/Dionysus/

Hints for compilation
{
    install dependencies:
        libboost-all-dev

    ccmake:
        use_cgal = OFF
        Boost_PYTHON_LIBRARY_RELEASE = /usr/lib/x86_64-linux-gnu/libboost_python-pyXX.so   (XX =  version of python)
        Boost_PYTHON_LIBRARY_DEBUG = /usr/lib/x86_64-linux-gnu/libboost_python-pyXX.so     (XX =  version of python)
        PYTHON_INCLUDE_DIR = /usr/include/pythonXX  (XX =  version of python)
        PYTHON_LIBRARY =  /usr/lib/x86_64-linux-gnu/libpythonXX.so  (XX =  version of python)
}

Copy the file _dionysus.so from your_build_dir/bindings/python/dionysus to the directory of THIS file.

Do NOT remove the build, since _dionysus.so links to resources there!
"""

def PersistenceDiagram(dimension: int, points: list):
    if not isinstance(dimension, int):
        raise ValueError('dimension expected to be int but got {}'.format(type(dimension)))

    try:
        points = [tuple(p) for p in points]

    except Exception as ex:
        raise ValueError('points expected to be list of tuples.') from ex

    return _dionysus.PersistenceDiagram(dimension, points)


def _check_for_persistence_dgm(param_name: str, object):
    if not isinstance(object, _dionysus.PersistenceDiagram):
        try:
            old_object = object
            object = PersistenceDiagram(1, object)

        except Exception as ex:
            raise ValueError('{} exptected to be PersistenceDiagram but got {}'.format(param_name, type(old_object))) from ex

    return object


def wasserstein_distance(dgm_1: list, dgm_2: list, p: int=2):
    if not isinstance(p, int):
        raise ValueError('dimension expected to be int but got {}'.format(type(p)))

    dgm_1 = _check_for_persistence_dgm('dgm_1', dgm_1)
    dgm_2 = _check_for_persistence_dgm('dgm_2', dgm_2)

    return _dionysus.wasserstein_distance(dgm_1, dgm_2, p)


def bottleneck_distance(dgm_1, dgm_2):
    dgm_1 = _check_for_persistence_dgm('dgm_1', dgm_1)
    dgm_2 = _check_for_persistence_dgm('dgm_2', dgm_2)

    return _dionysus.bottleneck_distance(dgm_1, dgm_2)
