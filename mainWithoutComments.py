from collections import defaultdict
from sympy import Matrix, lcm

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
    element_count = defaultdict(int)
    i = 0
    while i < len(molecule):
        if molecule[i] == "(":
            stack.append(element_count)
            element_count = defaultdict(int)
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
                prev[k] += v * factor
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
            element_count[elem] += int(num) if num else 1
    return dict(element_count)

def parse_side(side):
    molecules = side.split("+")
    return [parse_molecule(m.strip()) for m in molecules]

def balance_equation(equation):
    left, right = equation.split("->")
    left_moles = parse_side(left)
    right_moles = parse_side(right)

    elements = set()
    for m in left_moles + right_moles:
        elements.update(m.keys())
    elements = sorted(elements)


    matrix = []
    for i in elements:
        row = []
        for m in left_moles:
            row.append(m.get(i, 0))
        for m in right_moles:
            row.append(-m.get(i, 0))
        matrix.append(row)

    matrix = Matrix(matrix)
    solution = matrix.nullspace()
    if not solution:
        raise ValueError("Cannot balance equation")

    coefficients = solution[0]
    lcm_val = lcm([r.q for r in coefficients])
    coeffs = [int(r * lcm_val) for r in coefficients]

    def _gcd(a, b):
        a, b = abs(a), abs(b)
        while b:
            a, b = b, a % b
        return a

    g = 0
    for x in coeffs:
        g = _gcd(g, x)
    if g > 1:
        coeffs = [x // g for x in coeffs]

    n_left = len(left_moles)
    left_coeffs = coeffs[:n_left]
    right_coeffs = coeffs[n_left:]

    def all_positive(lst):
        return all(x > 0 for x in lst)

    if not (all_positive(left_coeffs) and all_positive(right_coeffs)):
        neg = [-x for x in coeffs]
        left_n = neg[:n_left]
        right_n = neg[n_left:]
        if all_positive(left_n) and all_positive(right_n):
            coeffs = neg
            left_coeffs = left_n
            right_coeffs = right_n
        else:
            raise ValueError("Cannot balance equation without moving species across the arrow")

    left_str = " + ".join(f"{c if c>1 else ''}{m}" for c, m in zip(left_coeffs, left.split("+")))
    right_str = " + ".join(f"{c if c>1 else ''}{m}" for c, m in zip(right_coeffs, right.split("+")))

    return f"{left_str} -> {right_str}"


print("Enter a chemical equation using '->'. Example: H2 + O2 -> H2O")
equation = input("Equation: ")

try:
    balanced = balance_equation(equation)
    print("Balanced equation: ")
    print(balanced)
except Exception as e:
    print("Error:", e)
