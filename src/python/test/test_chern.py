import unittest
import numpy as np
from chern import Hamiltonian, Chern, check


class TestHamiltonian(unittest.TestCase):
    def test_create_hamiltonian(self):
        hk = Hamiltonian()
        kx = np.random.rand() * 2 * np.pi - np.pi
        ky = np.random.rand() * 2 * np.pi - np.pi
        self.assertTrue(hk(kx, ky) == 0.)
        self.assertTrue(hk.dim == 1)
        self.assertTrue(Chern(hk).chern == [0.])


class TestChern(unittest.TestCase):

    # Todo: test check function
    def test_check(self):
        pass

    def test_discretize(self):
        cn = Chern()
        qq = cn.discretize()
        self.assertEqual(qq[0], -np.pi)
        self.assertEqual(qq[-1], np.pi)
        dq = qq[1] - qq[0]
        self.assertEqual(dq, 2 * np.pi / cn.n_discretized)
        self.assertTrue(np.allclose(np.diff(qq), [dq] * cn.n_discretized))

    # Todo: test calc_chern function
    def test_calc_chern(self):
        pass


class TestHaldaneHoneycomb(unittest.TestCase):
    """test Haldane honeycomb. """

    def haldane_honeycomb(self, kx, ky, m=0.5, phi=np.pi/2):

        k = np.array([kx / np.sqrt(3.), ky * 2. / 3.])

        t1 = t2 = 1.

        a1 = np.array([np.sqrt(3) * 0.5, 0.5])
        a2 = np.array([0, -1])
        a3 = np.array([-np.sqrt(3) * 0.5, 0.5])

        b1 = a2 - a3
        b2 = a3 - a1
        b3 = a1 - a2

        pauli0 = np.eye(2)
        pauli1 = np.array([[0, 1], [1, 0]])
        pauli2 = np.array([[0, -1j], [1j, 0]])
        pauli3 = np.array([[1, 0], [0, -1]])

        hk = 2 * t2 * np.cos(phi) * (
                np.cos(k @ b1) + np.cos(k @ b2) + np.cos(k @ b3)
        ) * pauli0 + t1 * (
                (np.cos(k @ a1) + np.cos(k @ a2) + np.cos(k @ a3)) * pauli1 +
                (np.sin(k @ a1) + np.sin(k @ a2) + np.sin(k @ a3)) * pauli2
        ) + (m - 2 * t2 * np.sin(phi) * (
                np.sin(k @ b1) + np.sin(k @ b2) + np.sin(k @ b3)
        )) * pauli3

        return hk

    def test_create_haldane(self):
        hk = Hamiltonian(self.haldane_honeycomb, "haldane")
        self.assertEqual(check(hk), 0)
        kx = np.random.rand() * 2 * np.pi - np.pi
        ky = np.random.rand() * 2 * np.pi - np.pi
        self.assertTrue(np.allclose(self.haldane_honeycomb(kx, ky), hk(kx, ky)))

    def test_haldane_chern(self):
        hk = Hamiltonian(self.haldane_honeycomb, "haldane")
        self.assertTrue(np.allclose(Chern(hk).chern, [-1, 1]))


class TestReal(unittest.TestCase):
    """test real Hamiltonian. """

    def realhk(self, kx, ky):
        mat = np.array([
            [3 - np.cos(kx) - np.cos(ky), np.sin(kx), np.sin(ky)],
            [np.sin(kx), np.cos(kx) * np.cos(ky), np.sin(kx) * np.sin(ky)],
            [np.sin(ky), np.sin(kx) * np.sin(ky), 1]
        ])
        return mat

    def test_create_real(self):
        hk = Hamiltonian(self.realhk, name="realhk")
        self.assertTrue(check(hk), 1)
        kx = np.random.rand() * 2 * np.pi - np.pi
        ky = np.random.rand() * 2 * np.pi - np.pi
        self.assertTrue(np.allclose(self.realhk(kx, ky), hk(kx, ky)))

    def test_real_chern(self):
        hk = Hamiltonian(self.realhk, name="realhk")
        self.assertTrue(np.allclose(Chern(hk).chern, np.zeros(3)))
