"""Independent graph pairing oracle retained from the exact validation suite."""

from __future__ import annotations

from functools import cache

Pairing = tuple[int, ...]

def components(adjacency: list[list[int]]) -> list[set[int]]:
    """Connected components by graph traversal, without a disjoint-set library."""
    unseen = set(range(len(adjacency)))
    result = []
    while unseen:
        start = unseen.pop()
        found = {start}
        stack = [start]
        while stack:
            for neighbor in adjacency[stack.pop()]:
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    found.add(neighbor)
                    stack.append(neighbor)
        result.append(found)
    return result

def connect(adjacency: list[list[int]], left: int, right: int) -> None:
    adjacency[left].append(right)
    adjacency[right].append(left)

def validate_pairing(diagram: Pairing) -> int:
    """Require a perfect matching; top and bottom are each indexed left to right."""
    if len(diagram) % 2 or not diagram:
        raise ValueError("a pairing must have a positive even number of endpoints")
    for index, partner in enumerate(diagram):
        if (type(partner) is not int or not 0 <= partner < len(diagram)
                or partner == index or diagram[partner] != index):
            raise ValueError("pairing must be a fixed-point-free involution")
    return len(diagram) // 2

@cache
def identity(strands: int) -> Pairing:
    return tuple(range(strands, 2 * strands)) + tuple(range(strands))

@cache
def cup_cap(strands: int, index: int) -> Pairing:
    """The generator e_index, with index starting at one."""
    if not 1 <= index < strands:
        raise ValueError("generator index must be between one and strands-1")
    result = list(identity(strands))
    left, right = index - 1, index
    result[left], result[right] = right, left
    result[strands + left], result[strands + right] = strands + right, strands + left
    return tuple(result)

@cache
def compose(upper: Pairing, lower: Pairing) -> tuple[Pairing, int]:
    """Stack two matchings; return the external matching and internal loop count."""
    strands = validate_pairing(upper)
    if validate_pairing(lower) != strands:
        raise ValueError("cannot compose diagrams of different widths")
    # The shared middle boundary has indices n..2n-1; external bottom is 2n..3n-1.
    adjacency = [[] for _ in range(3 * strands)]
    for offset, diagram in ((0, upper), (strands, lower)):
        for endpoint, partner in enumerate(diagram):
            if endpoint < partner:
                connect(adjacency, endpoint + offset, partner + offset)
    result = [-1] * (2 * strands)
    loops = 0
    for connected in components(adjacency):
        boundary = sorted(
            vertex if vertex < strands else vertex - strands
            for vertex in connected if vertex < strands or vertex >= 2 * strands
        )
        if not boundary:
            loops += 1
        elif len(boundary) == 2:
            left, right = boundary
            result[left], result[right] = right, left
        else:
            raise AssertionError(f"composition component has {len(boundary)} endpoints")
    matching = tuple(result)
    validate_pairing(matching)
    return matching, loops

@cache
def trace_loops(diagram: Pairing) -> int:
    """Join each top endpoint to the bottom endpoint in the same position."""
    strands = validate_pairing(diagram)
    adjacency = [[] for _ in range(2 * strands)]
    for endpoint, partner in enumerate(diagram):
        if endpoint < partner:
            connect(adjacency, endpoint, partner)
    for endpoint in range(strands):
        connect(adjacency, endpoint, endpoint + strands)
    return len(components(adjacency))
