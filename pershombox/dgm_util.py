def de_essentialize(persistence_diagram: [[]], maximal_filtration_value: float)->[[]]:
    """
    Replaces all essental classes with points dying at  maximal_filtration_value.
    Example:
        de_essentialize([(1, 2), (0, inf)], 10) = [(1, 2), (0, 10)]

    :param persistence_diagram:
    :param maximal_filtration_value: value which replaces inf
    :return:
    """
    if len(persistence_diagram) == 0:
        return persistence_diagram

    else:
        dgm = [tuple(p) for p in persistence_diagram]

        de_essentialized_dgm = [p for p in dgm if p[1] != float('inf')] + \
                               [(p[0], maximal_filtration_value) for p in dgm if p[1] == float('inf')]

        return de_essentialized_dgm
