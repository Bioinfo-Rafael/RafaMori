"""Geometry primitives for the configured target coordinates."""
from __future__ import annotations

import numpy as np


ATOMIC_MASSES = {
    1: 1.00784, 5: 10.81, 6: 12.011, 7: 14.007, 8: 15.999, 9: 18.998403,
    14: 28.085, 15: 30.973762, 16: 32.06, 17: 35.45, 35: 79.904, 53: 126.90447,
}


def distance(point_xyz: dict[int, np.ndarray], atom_a: int, atom_b: int) -> float:
    return float(np.linalg.norm(point_xyz[atom_a] - point_xyz[atom_b]))


def angle(point_xyz: dict[int, np.ndarray], atom_a: int, vertex: int, atom_c: int) -> float:
    vector_a = point_xyz[atom_a] - point_xyz[vertex]
    vector_c = point_xyz[atom_c] - point_xyz[vertex]
    denominator = np.linalg.norm(vector_a) * np.linalg.norm(vector_c)
    if denominator == 0:
        return np.nan
    cosine = np.dot(vector_a, vector_c) / denominator
    return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))


def center_of_mass(point_xyz: dict[int, np.ndarray], atomic_numbers: dict[int, int], atom_numbers: list[int]) -> np.ndarray:
    masses = np.asarray([ATOMIC_MASSES.get(atomic_numbers[number], float(atomic_numbers[number])) for number in atom_numbers])
    coordinates = np.asarray([point_xyz[number] for number in atom_numbers])
    return np.average(coordinates, axis=0, weights=masses)


def com_distance(point_xyz: dict[int, np.ndarray], atomic_numbers: dict[int, int], fragment_a: list[int], fragment_b: list[int]) -> float:
    return float(np.linalg.norm(center_of_mass(point_xyz, atomic_numbers, fragment_a) - center_of_mass(point_xyz, atomic_numbers, fragment_b)))
