from collections import defaultdict
from sympy import Matrix, lcm
import re

# ----------------------------
# CHEMICAL ELEMENT LIST
# ----------------------------
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

# ----------------------------
# PARSE MOLECULE
# ----------------------------
def parse_molecule(molecule):
    """
    Parse a molecule string like "Ca(OH)2" into a dict of element counts.
    """
    def multiply_dict(d, factor):
        return {k: v * factor for k, v in d.items()}

    stack = []
    elem_count = defaultdict(int)
    i = 0
    while i < len(molecule):
        if molecule[i] == "(":
            stack.append(elem_count)
            elem_count = defaultdict(int)
            i += 1
        elif molecule[i] == ")":
            i += 1
            num = ""
            while i < len(molecule) and molecule[i].isdigit():
                num += molecule[i]
                i += 1
            factor = int(num) if num else 1
            prev = stack.pop()
            for k, v in elem_count.items():
                prev[k] += v * factor
            elem_count = prev
        else:
            # parse element
            if i+1 < len(molecule) and molecule[i+1].islower():
                elem = molecule[i:i+2]
                i += 2
            else:
                elem = molecule[i]
                i += 1
            if elem not in ELEMENTS:
                raise ValueError(f"Unknown element: {elem}")
            # parse number
            num = ""
            while i < len(molecule) and molecule[i].isdigit():
                num += molecule[i]
                i += 1
            elem_count[elem] += int(num) if num else 1
    return dict(elem_count)

# ----------------------------
# PARSE SIDE
# ----------------------------
def parse_side(side):
    """
    Split a side like "H2 + O2" into molecule dicts
    """
    molecules = side.split("+")
    return [parse_molecule(m.strip()) for m in molecules]

# ----------------------------
# BALANCE EQUATION
# ----------------------------
def balance_equation(equation):
    left, right = equation.split("->")
    left_mols = parse_side(left)
    right_mols = parse_side(right)

    # list all elements
    elements = set()
    for m in left_mols + right_mols:
        elements.update(m.keys())
    elements = sorted(elements)

    # build matrix
    matrix = []
    for el in elements:
        row = []
        for m in left_mols:
            row.append(m.get(el, 0))
        for m in right_mols:
            row.append(-m.get(el, 0))
        matrix.append(row)

    M = Matrix(matrix)
    sol = M.nullspace()
    if not sol:
        raise ValueError("Cannot balance equation")

    coeffs = sol[0]
    lcm_val = lcm([r.q for r in coeffs])
    coeffs = [int(r*lcm_val) for r in coeffs]

    left_coeffs = coeffs[:len(left_mols)]
    right_coeffs = coeffs[len(left_mols):]

    # build balanced equation
    left_str = " + ".join(f"{c if c>1 else ''}{m}" for c, m in zip(left_coeffs, left.split("+")))
    right_str = " + ".join(f"{c if c>1 else ''}{m}" for c, m in zip(right_coeffs, right.split("+")))

    return f"{left_str} -> {right_str}"

# ----------------------------
# MAIN LOOP
# ----------------------------
print("Enter a chemical equation using '->'. Example: H2 + O2 -> H2O")
equation = input("Equation: ")

try:
    balanced = balance_equation(equation)
    print("Balanced equation:")
    print(balanced)
except Exception as e:
    print("Error:", e)