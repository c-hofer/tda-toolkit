"""
This module implements high level access to the functionality of the DIPHA project. https://github.com/pkuwwt/dipha

Hints for compilation of DIPHA
{
install openmpi-bin    (not lib-openmpi-dev)
}
"""
import io, os, struct
import numpy
import functools
import subprocess
import warnings

from subprocess import DEVNULL
from tempfile import TemporaryDirectory


# region module resources


__ext_lib_path = os.path.join(os.path.dirname(__file__), 'ext_lib')

__dipha_path = os.path.join(__ext_lib_path, 'dipha')


if not os.path.isfile(__dipha_path):
    try:
        subprocess.call('dipha', stdout=DEVNULL, stderr=DEVNULL)

    except FileNotFoundError:

        text = """Found no dipha executable in  {}. Get it from
                http://people.maths.ox.ac.uk/nanda/perseus/index.html and copy
                it into {} or set the __dipha_path variable
                of this module manually after import.
                """.format(__dipha_path, __ext_lib_path)

        warnings.warn(text, ImportWarning)


__tmp_dir_fact = TemporaryDirectory


# endregion

# region Helpers


class DiphaAdapterException(Exception):
    pass


def _pack(value, type: type):
    if type is float:
        return struct.pack('d', value)

    if type is int:
        return struct.pack('<q', value)

    raise ValueError('Packing is not implemented for type {}.'.format(type(value)))


def _unpack(bytes_to_unpack, type: type):
    if type is float:
        return struct.unpack('d', bytes_to_unpack)

    if type is int:
        return struct.unpack('<q', bytes_to_unpack)

    raise ValueError('Unpacking is not implemented for type {}.'.format(type))


def _unpack_from_file(f: io.BufferedReader, type):
    read_bytes = f.read(8)
    return _unpack(read_bytes, type)[0]


def _run_dipha(input_file, output_file, limit_dimensions: int=None, dual: bool=False, benchmark: bool=False):
    try:
        args = []

        if limit_dimensions is not None:
            args.append("--upper_dim")
            args.append(str(int(limit_dimensions)))

        if dual:
            args.append("--dual")

        if benchmark:
            args.append("--benchmark")

        args += [input_file, output_file]

        p = subprocess.Popen(["dipha",  *args], stdout=DEVNULL, stderr=DEVNULL)
        p.wait()
    except FileNotFoundError:
        raise DiphaAdapterException("The program dipha seems not to be installed. "
                                    "Get it from https://github.com/DIPHA/dipha")

# endregion

# region file classes


class _DIPHAFile:
    """
    Abstract base class for all DIPHA file types.
    """
    _magic_number = 8067171840

    def __init__(self):
        pass

    def write_to_binary_file(self, f):
        f.write(_pack(self._magic_number, int))

    @staticmethod
    def load_from_binary_file(f):
        read_magic_number = _unpack_from_file(f, int)

        if read_magic_number != _DIPHAFile._magic_number:
            raise ValueError("Argument is not a valid DIPHA file.")


class _ImageDataFile(_DIPHAFile):
    """

    """
    _file_type_number = 1

    def __init__(self, img: numpy.array):
        """

        :param img: The array specifies an image in n = img.ndim dimensions. Where
        img[x_1, ... ,x_n] specifies the value at (x_1, ... ,x_n) in canonical coordinates.
        (This means the image is assumed to be defined in the positive quadrant of the coordinate system)
        """
        super(_ImageDataFile, self).__init__()
        img = numpy.array(img)

        self._img = img

    def write_to_binary_file(self, f):
        # Write header
        super(_ImageDataFile, self).write_to_binary_file(f)   #Base class part

        number_of_data_values = functools.reduce(lambda x, y: x*y, self._img.shape, 1)
        dimension = self._img.ndim
        lattice_resolutions = self._img.shape

        header = [_ImageDataFile._file_type_number, number_of_data_values, dimension, *lattice_resolutions]
        for value in header:
            f.write(_pack(value, int))

        # Write data
        values = []
        if self._img.ndim == 2:
            for column_number in range(self._img.shape[1]):
                values += list(self._img[:, column_number])

        elif self._img.ndim == 3:
            for z in range(self._img.shape[2]):
                z_slice = self._img[:, :, z]
                for column_number in range(z_slice.shape[1]):
                    values += list(self._img[:, column_number, z])

        else:
            raise DiphaAdapterException('Only image data of 2 or 3 dimensions possible!')

        for value in values:
            f.write(_pack(value, float))

    @staticmethod
    def load_from_binary_file(f):
        raise NotImplementedError()


class _DistanceMatrixFile(_DIPHAFile):
    """

    """
    _file_type_number = 7

    def __init__(self, distance_matrix: numpy.array):
        super().__init__()
        self._distance_matrix = numpy.array(distance_matrix)

    def write_to_binary_file(self, f):
        # Write header
        super(_DistanceMatrixFile, self).write_to_binary_file(f)  # Base class part

        f.write(_pack(self._file_type_number, int))

        number_of_points = self._distance_matrix.shape[0]
        f.write(_pack(number_of_points, int))

        for i in range(number_of_points):
            for j in range(number_of_points):
                f.write(_pack(self._distance_matrix[i, j], float))

    @staticmethod
    def load_from_binary_file(f):
        raise NotImplementedError()


class _PersistenceDiagramFile(_DIPHAFile):
    """

    """
    _file_type_number = 2

    def __init__(self, points):
        super(_PersistenceDiagramFile, self).__init__()
        self.points = list(points)

    def write_to_binary_file(self, f):
        raise NotImplementedError()

    @staticmethod
    def load_from_binary_file(f):
        super(_PersistenceDiagramFile, _PersistenceDiagramFile).load_from_binary_file(f)

        file_type_number = _unpack_from_file(f, int)
        if file_type_number != _PersistenceDiagramFile._file_type_number:
            raise ValueError("Argument is not a Persistence Diagram DIPHA file.")

        number_of_points = _unpack_from_file(f, int)
        points = []
        for i in range(number_of_points):
            dimension = _unpack_from_file(f, int)
            birth_time = _unpack_from_file(f, float)
            death_time = _unpack_from_file(f, float)
            points.append((dimension, birth_time, death_time))

        return _PersistenceDiagramFile(points)


# endregion


# region public functional interface


def persistence_diagrams_of_filtrated_cubical_complex(filtrated_cubical_complex: numpy.array,
                                                      limit_dimensions: int=None,
                                                      dual: bool=False,
                                                      benchmark: bool=False,
                                                      set_inf_to_max_filt_val=False)->[[tuple]]:
    """
    Calculates the persistence diagram for a cubical complex.

    :param filtrated_cubical_complex: An n dimensional array, such that

        filtrated_cubical_complex[x_1, ... , x_d] = f(x_1, ... , x_d)

    where f is the filtration and x_i are the coordinates of the vertex with respect to the canonical basis in the
    positive quadrant on the unit spaced grid.

    :return:
    List with the points of the persistence diagram of dimension k at position k.
    """
    filtrated_cubical_complex = numpy.array(filtrated_cubical_complex)
    dimension = filtrated_cubical_complex.ndim

    with __tmp_dir_fact() as tmp_dir:

        image_data_file_path = os.path.join(tmp_dir, "image_data")
        persistence_diagram_file_path = os.path.join(tmp_dir, "persistence_diagram")

        with open(image_data_file_path, "bw") as f:
            _ImageDataFile(filtrated_cubical_complex).write_to_binary_file(f)

        _run_dipha(image_data_file_path,
                   persistence_diagram_file_path,
                   limit_dimensions,
                   dual,
                   benchmark)

        with open(persistence_diagram_file_path, "rb") as f:
            diagram = _PersistenceDiagramFile.load_from_binary_file(f)

        tmp = {}
        for i in range(dimension):
            tmp[i] = []

        for dimension, birth, death in diagram.points:
            if dimension < 0:
                dimension = -dimension - 1
                if not set_inf_to_max_filt_val:
                    death = float('inf')

            tmp[dimension].append((birth, death))

        return [tmp[key] for key in tmp.keys()]


def persistence_diagrams_of_VR_complex_from_distance_matrix(distance_matrix: numpy.array,
                                                            upper_dimension: int,
                                                            dual: bool = False,
                                                            benchmark: bool = False) -> [[tuple]]:
    distance_matrix = numpy.array(distance_matrix)

    with __tmp_dir_fact() as tmp_dir:
        distance_matrix_file_path = os.path.join(tmp_dir, "distance_matrix")
        persistence_diagram_file_path = os.path.join(tmp_dir, "persistence_diagram")

        with open(distance_matrix_file_path, "bw") as f:
            _DistanceMatrixFile(distance_matrix).write_to_binary_file(f)

        _run_dipha(distance_matrix_file_path,
                   persistence_diagram_file_path,
                   upper_dimension,
                   dual,
                   benchmark)

        with open(persistence_diagram_file_path, "rb") as f:
            diagram = _PersistenceDiagramFile.load_from_binary_file(f)

        tmp = {}
        for i in range(upper_dimension):
            tmp[i] = []

        for dimension, birth, death in diagram.points:
            if dimension < 0:
                dimension = -dimension - 1
                death = float('inf')

            tmp[dimension].append((birth, death))

        return [tmp[key] for key in tmp.keys()]


# endregion



