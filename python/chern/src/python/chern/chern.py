import numpy as np
import inspect


def triv(kx, ky):
    return np.array([[0]])


def check(hkfunc):
    """
    Check the input if callable is
        * square matrix on random input hkfunc(kx, ky)
        * hermitian
        * real
    :param hkfunc: callable to accept two positional args.
    :return: 0 if not symmetric; 1 if symmetric.
    """
    xx = np.random.rand(32) * 2 * np.pi - np.pi
    yy = np.random.rand(32) * 2 * np.pi - np.pi
    dim = np.array(hkfunc(0, 0)).shape[0]
    for x in xx:
        for y in yy:
            ham = np.array(hkfunc(x, y))
            assert ham.shape == (dim, dim), "Not square matrix"
            assert np.allclose(ham.conjugate().T, ham), "Not hermitian matrix"
            if not np.allclose(ham.T, ham):
                return 0
    return 1


class Hamiltonian(object):
    """Hamiltonian container.
    """

    def __init__(self, hk=triv, name="myhk"):
        assert callable(hk)
        r = check(hk)
        if r == 0:
            self._real = False
        elif r == 1:
            self._real = True
        self._name = name
        self._hk = hk
        self._dim = len(hk(0, 0))

    @property
    def name(self):
        return self._name

    @property
    def dim(self):
        return self._dim

    @property
    def real(self):
        return self._real

    def __repr__(self):
        return f"Hamiltonian({self.name}:{self.dim}-bands)"

    def __call__(self, *args, **kwargs):
        return self._hk(*args, **kwargs)


class Chern(object):
    """Chern container.
    """

    def __init__(self, hk=None, n_discretized=16):
        self._hk = hk
        self._discrete = n_discretized
        self._qq = self.discretize()

    def __repr__(self):
        return f"{type(self).__name__}({self.hk.name}:{self.n_discretized}x{self.n_discretized})"

    @property
    def hk(self):
        return self._hk

    @property
    def n_discretized(self):
        return self._discrete

    @property
    def qq(self):
        return self._qq

    def discretize(self):
        nq = self.n_discretized
        dq = 2. * np.pi / nq
        qq = np.arange(-np.pi, np.pi + dq, dq)
        return qq

    def diagn(self):
        """
        diagonalize self.hk Hamiltonian across all discretized brillioun zone.
        :return: corresponding eigen-energy and eigen-states on each (kx, ky),
        sorted ascendant.
        """
        dim = len(self.qq)
        # hErg = np.empty((dim, dim, self.hk.dim))
        # hVec = np.empty((dim, dim, self.hk.dim, self.hk.dim))
        hErg = []
        hVec = []
        for iqy in range(dim):
            qy = self.qq[iqy]
            for iqx in range(dim):
                qx = self.qq[iqx]
                ham = self.hk(qx, qy)
                eigsys = diagn_sort(ham)
                hErg.append(np.array([erg for erg, _ in eigsys]))
                hVec.append(np.array([vec for _, vec in eigsys]).T)
        hErg = np.array(hErg).reshape(dim, dim, -1)
        hVec = np.array(hVec).reshape(dim, dim, self.hk.dim, -1)
        return hErg, hVec

    def calc_chern(self):
        """
        Calculate the 1st Chern number of given Hamiltonian. (self.hk).
        Using Japanese algorithm (2005).
        :return: list of Chern number of all bands. Should summed to 0.
        """
        if self.hk is None:
            print("Hamiltonian not given")
            return None
        if not isinstance(self.hk, Hamiltonian):
            self._hk = Hamiltonian(self.hk)
        if self.hk.real:
            return np.zeros(self.hk.dim)
        if not hasattr(self, "_hVec"):
            _, self._hVec = self.diagn()
        return calc_chern(self._hVec)

    @property
    def chern(self):
        self._chern = self.calc_chern()
        return self._chern


def diagn_sort(mat):
    """
    Diagonalize matrix.
    :param mat: matrix to be diagonalized.
    :return: ((val, vec), ...), sorted ascendant according to val.
    """
    vals, vecs = np.linalg.eig(mat)
    idx = vals.argsort()
    # idx = vals.argsort()[::-1]
    eigsys = tuple((vals[i], vecs[:, i]) for i in idx)
    return eigsys


def calc_chern(hVec):
    """
    Calculate 1st Chern number given eigen-vectors. Using Japanese algorithm (2005).
    :param hVec: eigen-vectors of all band (assuming gapless).
    :return: list of Chern number of all bands. Should summed to 0.
    """
    hVec = np.array(hVec)
    dimy, dimx, nlevels, _ = hVec.shape
    cnlist = np.zeros(nlevels)

    for iy in range(dimy-1):
        for ix in range(dimx-1):
            u12 = hVec[iy, ix + 1].conjugate().T @ hVec[iy, ix]
            u23 = hVec[iy + 1, ix + 1].conjugate().T @ hVec[iy, ix + 1]
            u34 = hVec[iy + 1, ix].conjugate().T @ hVec[iy + 1, ix + 1]
            u41 = hVec[iy, ix].conjugate().T @ hVec[iy + 1, ix]
            t12 = np.diag(u12.diagonal())
            t23 = np.diag(u23.diagonal())
            t34 = np.diag(u34.diagonal())
            t41 = np.diag(u41.diagonal())
            tplaquet = t41 @ t34 @ t23 @ t12
            cnlist += np.angle(tplaquet.diagonal())

    cnlist /= 2 * np.pi
    cnlist = chop(cnlist)
    return cnlist


def chop(array, tol=1e-7):
    """
    realize Mathematica Chop[].
    :param array: 1D array.
    :param tol: tolerance to be chopped. default to 1e-7
    :return: chopped array. (original array alse modified.)
    """
    for i in range(len(array)):
        a = array[i]
        if np.abs(a-round(a)) < tol:
            array[i] = round(a)
    return array
