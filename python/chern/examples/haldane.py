import numpy as np


def haldane_honeycomb(kx, ky, m=0.5, phi=np.pi/2):

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


if __name__ == "__main__":
    from chern import Hamiltonian, Chern
    from functools import partial
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--m', type=float, help='haldane mass', default=0.5)
    parser.add_argument('--phi', type=float, help='phi flux', default=np.pi/2)
    args = parser.parse_args()
    hk = Hamiltonian(partial(haldane_honeycomb, m=args.m, phi=args.phi), "haldane")
    cn = Chern(hk)
    print(f"\nThe Chern number for haldane(m={args.m:5.2f},phi={args.phi:5.2f}) is:  {cn.chern}\n")
