# language: python
from fractions import Fraction
import math

ELEMENTS = [
    'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr',
    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
    'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
    'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
    'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
    'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
    'Rg', 'Cn', 'Fl', 'Lv', 'Ts', 'Og'
]

def parse_molecule(molecule):
    stack = []
    element_count = {}
    i = 0
    while i < len(molecule):
        if molecule[i] == "(":
            stack.append(element_count)
            element_count = {}
            i += 1
        elif molecule[i] == ")":
            i += 1
            num = ""
            while i < len(molecule) and molecule[i].isdigit():
                num += molecule[i]
                i += 1
            factor = int(num) if num else 1
            prev = stack.pop()
            for k, v in element_count.items():
                prev[k] = prev.get(k, 0) + v * factor
            element_count = prev
        else:
            if i+1 < len(molecule) and molecule[i+1].islower():
                elem = molecule[i:i+2]
                i += 2
            else:
                elem = molecule[i]
                i += 1
            if elem not in ELEMENTS:
                raise ValueError(f"Unknown element: {elem}")
            num = ""
            while i < len(molecule) and molecule[i].isdigit():
                num += molecule[i]
                i += 1
            element_count[elem] = element_count.get(elem, 0) + (int(num) if num else 1)
    return dict(element_count)

def parse_side(side):
    molecules = [m.strip() for m in side.split("+")]
    return molecules, [parse_molecule(m) for m in molecules]

# -- small rational linear algebra to get nullspace of integer matrix --
def rref(matrix):
    # matrix: list of rows (each row is list of Fractions)
    m = [row[:] for row in matrix]
    rows = len(m)
    cols = len(m[0]) if rows > 0 else 0
    r = 0
    pivots = []
    for c in range(cols):
        # find pivot
        pivot = None
        for i in range(r, rows):
            if m[i][c] != 0:
                pivot = i
                break
        if pivot is None:
            continue
        m[r], m[pivot] = m[pivot], m[r]
        pv = m[r][c]
        m[r] = [val / pv for val in m[r]]
        for i in range(rows):
            if i != r and m[i][c] != 0:
                factor = m[i][c]
                m[i] = [m[i][j] - factor * m[r][j] for j in range(cols)]
        pivots.append(c)
        r += 1
        if r == rows:
            break
    return m, pivots

def nullspace_basis_int(matrix):
    # matrix: list of rows of integers
    # returns one integer basis vector for nullspace (smallest integer vector), or [] if trivial
    if not matrix:
        return []
    rows = len(matrix)
    cols = len(matrix[0])
    # convert to Fractions
    A = [[Fraction(x) for x in row] for row in matrix]
    rref_m, pivots = rref(A)
    pivot_set = set(pivots)
    free_cols = [c for c in range(cols) if c not in pivot_set]
    if not free_cols:
        return []
    # construct one basis vector by setting first free var = 1, others = 0
    vec = [Fraction(0) for _ in range(cols)]
    free0 = free_cols[0]
    vec[free0] = Fraction(1)
    for row_idx, pcol in enumerate(pivots):
        s = Fraction(0)
        for fc in free_cols:
            s += rref_m[row_idx][fc] * vec[fc]
        vec[pcol] = -s
    # scale to integers
    den_lcm = 1
    for v in vec:
        den_lcm = den_lcm * v.denominator // math.gcd(den_lcm, v.denominator)
    int_vec = [int(v * den_lcm) for v in vec]
    # normalize by gcd and make first nonzero positive
    g = 0
    for x in int_vec:
        g = math.gcd(g, abs(x))
    if g == 0:
        return []
    int_vec = [x // g for x in int_vec]
    # ensure all coefficients non-negative (or at least first nonzero positive)
    for x in int_vec:
        if x != 0:
            if x < 0:
                int_vec = [-y for y in int_vec]
            break
    return int_vec

def balance_equation(equation):
    if "->" not in equation:
        raise ValueError("Equation must contain '->'")
    left_str, right_str = equation.split("->")
    left_mols, left_dicts = parse_side(left_str)
    right_mols, right_dicts = parse_side(right_str)

    elements = set()
    for m in left_dicts + right_dicts:
        elements.update(m.keys())
    elements = sorted(elements)

    # build integer matrix (rows = elements, cols = molecules)
    matrix = []
    for el in elements:
        row = []
        for m in left_dicts:
            row.append(m.get(el, 0))
        for m in right_dicts:
            row.append(-m.get(el, 0))
        matrix.append(row)

    coeffs = nullspace_basis_int(matrix)
    if not coeffs:
        raise ValueError("Cannot balance equation")

    left_coeffs = coeffs[:len(left_mols)]
    right_coeffs = coeffs[len(left_mols):]

    def fmt(coeff, mol):
        coeff_str = f"{coeff}" if coeff != 1 else ""
        return f"{coeff_str}{mol.strip()}"

    left_out = " + ".join(fmt(c, m) for c, m in zip(left_coeffs, left_mols))
    right_out = " + ".join(fmt(c, m) for c, m in zip(right_coeffs, right_mols))
    return f"{left_out} -> {right_out}"

# main
print("Enter a chemical equation using '->'. Example: H2 + O2 -> H2O")
equation = input("Equation: ")

try:
    balanced = balance_equation(equation)
    print("Balanced equation: ")
    print(balanced)
except Exception as e:
    print("Error:", e)
