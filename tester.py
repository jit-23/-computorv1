#!/usr/bin/env python3
"""
computor v1 tester

usage:
    python3 tester.py ./computor
    python3 tester.py python3 computor.py
    python3 tester.py --bonus python3 computor.py     # + free form, cubics, bad input
    python3 tester.py -v python3 computor.py          # show output of every test

What it checks per test:
    - no crash (traceback / non-zero exit / timeout)
    - "Reduced form: ... = 0" is mathematically correct (formatting-agnostic)
    - "Polynomial degree: N" is correct
    - the discriminant message (positive / zero / negative) or special message
      (any real number / no solution / can't solve)
    - the solutions, numerically (order-free, tolerance 1e-4 relative).
      Accepts floats, fractions (-1/5), and complex forms like
      "-0.2 + 0.4i", "-1/5 + 2i/5", "(-1+3j)", "x1 = 2.0"
    - warns if "-0" / "-0.0" is printed
"""
import ast
import cmath
import re
import subprocess
import sys

TIMEOUT = 5
TOL = 1e-4

USE_COLOR = sys.stdout.isatty()
def c(code, s): return f"\033[{code}m{s}\033[0m" if USE_COLOR else s
OK, KO, WARN, BOLD = (lambda s: c("32", s)), (lambda s: c("31", s)), (lambda s: c("33", s)), (lambda s: c("1", s))


# ───────────────────────── helpers ─────────────────────────

def quad(a, b, c_):
    """Reference roots of a*x^2 + b*x + c (tester only)."""
    d = cmath.sqrt(b * b - 4 * a * c_)
    return [(-b + d) / (2 * a), (-b - d) / (2 * a)]


def safe_eval(expr):
    """Evaluate a numeric expression (+ - * / and complex literals) without eval()."""
    def ev(n):
        if isinstance(n, ast.Constant) and isinstance(n.value, (int, float, complex)) \
                and not isinstance(n.value, bool):
            return n.value
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, (ast.UAdd, ast.USub)):
            v = ev(n.operand)
            return v if isinstance(n.op, ast.UAdd) else -v
        if isinstance(n, ast.BinOp) and isinstance(n.op, (ast.Add, ast.Sub, ast.Mult, ast.Div)):
            l, r = ev(n.left), ev(n.right)
            if isinstance(n.op, ast.Add): return l + r
            if isinstance(n.op, ast.Sub): return l - r
            if isinstance(n.op, ast.Mult): return l * r
            return l / r
        raise ValueError("not a number")
    return ev(ast.parse(expr, mode="eval").body)


PREFIX = re.compile(r'^\s*(?:x|s|sol|solution|root|r)\s*_?\d*\s*[:=]\s*', re.I)

def parse_number(line):
    """Return a complex if the line is just a (possibly complex) number, else None."""
    low = line.lower()
    if "discriminant" in low or "delta" in low or "degree" in low or "reduced" in low:
        return None
    s = PREFIX.sub("", line).strip()
    if not s or "=" in s:
        return None
    # letters other than i/j (allow scientific notation like 1e-05)
    if re.search(r'[a-hk-z]', re.sub(r'\d[eE][+-]?\d', '', s), re.I):
        return None
    s = re.sub(r'\*\s*[ij]\b', '*1j', s)         # 2*i   -> 2*1j
    s = re.sub(r'(\d)\s*i', r'\1j', s)            # 2i/5  -> 2j/5
    s = re.sub(r'(?<![\w.])[ij](?!\w)', '1j', s)  # + i   -> + 1j
    try:
        return complex(safe_eval(s))
    except Exception:
        return None


def parse_poly(expr):
    """'4 * X^0 + 4 * X^1 - 9.3 * X^2' or '5 + 4 * X' -> {power: coef}. None if unparseable."""
    s = expr.replace(" ", "").lower()
    if not s:
        return None
    coeffs = {}
    for term in re.findall(r'[+-]?(?:[^+-]|(?<=\de)[+-])+', s):   # keep 1e+80 together
        m = re.fullmatch(r'([+-]?)(\d+(?:\.\d+)?(?:e[+-]?\d+)?)?(\*?x(?:\^(\d+))?)?', term)
        if not m or (m.group(2) is None and m.group(3) is None):
            return None
        coef = float(m.group(2)) if m.group(2) else 1.0
        if m.group(1) == "-":
            coef = -coef
        power = int(m.group(4)) if m.group(4) else (1 if m.group(3) else 0)
        coeffs[power] = coeffs.get(power, 0.0) + coef
    return coeffs


def close(a, b):
    return abs(a - b) <= TOL * max(1.0, abs(b))


def covers(got, expected):
    """Every printed root matches an expected one and vice versa (order-free)."""
    return all(any(close(g, e) for e in expected) for g in got) and \
           all(any(close(g, e) for g in got) for e in expected)


def fmt(z):
    z = complex(z)
    return f"{z.real:g}" if abs(z.imag) < 1e-12 else f"{z.real:g}{z.imag:+g}i"


# ───────────────────────── tests ─────────────────────────
# kind: pos | zero | neg | linear | all | none | high | error
# coeffs: expected reduced form {power: coef} (zero terms may be omitted)

BIG = "1" + "0" * 80  # 1e80, catches sqrt functions that don't converge

MANDATORY = [
    # --- subject examples
    dict(name="subject: positive discriminant", eq="5 * X^0 + 4 * X^1 - 9.3 * X^2 = 1 * X^0",
         coeffs={0: 4, 1: 4, 2: -9.3}, degree=2, kind="pos", roots=quad(-9.3, 4, 4)),
    dict(name="subject: degree 1", eq="5 * X^0 + 4 * X^1 = 4 * X^0",
         coeffs={0: 1, 1: 4}, degree=1, kind="linear", roots=[-0.25]),
    dict(name="subject: degree 3", eq="8 * X^0 - 6 * X^1 + 0 * X^2 - 5.6 * X^3 = 3 * X^0",
         coeffs={0: 5, 1: -6, 3: -5.6}, degree=3, kind="high",
         roots=[0.61559813, -0.30779906 + 1.16432357j, -0.30779906 - 1.16432357j]),
    dict(name="subject: any real number", eq="6 * X^0 = 6 * X^0",
         coeffs={}, degree=0, kind="all"),
    dict(name="subject: no solution", eq="10 * X^0 = 15 * X^0",
         coeffs={0: -5}, degree=0, kind="none"),
    dict(name="subject: complex solutions", eq="1 * X^0 + 2 * X^1 + 5 * X^2 = 0",
         coeffs={0: 1, 1: 2, 2: 5}, degree=2, kind="neg", roots=[-0.2 + 0.4j, -0.2 - 0.4j]),
    dict(name="subject: 42 = 42", eq="42 * X^0 = 42 * X^0",
         coeffs={}, degree=0, kind="all"),

    # --- discriminant cases
    dict(name="discriminant zero", eq="1 * X^0 - 2 * X^1 + 1 * X^2 = 0",
         coeffs={0: 1, 1: -2, 2: 1}, degree=2, kind="zero", roots=[1]),
    dict(name="discriminant zero at x=0", eq="0 * X^0 + 0 * X^1 + 3 * X^2 = 0",
         coeffs={2: 3}, degree=2, kind="zero", roots=[0]),
    dict(name="negative a, real roots", eq="4 * X^0 + 0 * X^1 - 1 * X^2 = 0",
         coeffs={0: 4, 2: -1}, degree=2, kind="pos", roots=[2, -2]),
    dict(name="b = 0, complex roots", eq="1 * X^0 + 0 * X^1 + 1 * X^2 = 0",
         coeffs={0: 1, 2: 1}, degree=2, kind="neg", roots=[1j, -1j]),
    dict(name="negative a, complex roots", eq="1 * X^0 + 2 * X^1 - 1 * X^2 = 4 * X^0",
         coeffs={0: -3, 1: 2, 2: -1}, degree=2, kind="neg", roots=quad(-1, 2, -3)),
    dict(name="decimal coefficients", eq="0.5 * X^0 - 1.5 * X^1 + 1 * X^2 = 0",
         coeffs={0: 0.5, 1: -1.5, 2: 1}, degree=2, kind="pos", roots=[1, 0.5]),
    dict(name="terms mostly on right side", eq="2 * X^0 = 1 * X^0 + 3 * X^1 - 1 * X^2",
         coeffs={0: 1, 1: -3, 2: 1}, degree=2, kind="pos", roots=quad(1, -3, 1)),
    dict(name="huge coefficient (sqrt precision)", eq=f"1 * X^2 + {BIG} * X^1 = 0",
         coeffs={1: float(BIG), 2: 1}, degree=2, kind="pos", roots=[0, -float(BIG)]),

    # --- degree reduction
    dict(name="X^2 cancels -> degree 1", eq="1 * X^0 + 2 * X^1 + 3 * X^2 = 3 * X^2",
         coeffs={0: 1, 1: 2}, degree=1, kind="linear", roots=[-0.5]),
    dict(name="X^3 cancels -> degree 2", eq="1 * X^3 + 1 * X^2 - 4 * X^0 = 1 * X^3",
         coeffs={0: -4, 2: 1}, degree=2, kind="pos", roots=[2, -2]),
    dict(name="X^3 cancels -> degree 1", eq="1 * X^3 + 1 * X^1 = 1 * X^3",
         coeffs={1: 1}, degree=1, kind="linear", roots=[0]),
    dict(name="everything cancels", eq="3 * X^2 = 3 * X^2",
         coeffs={}, degree=0, kind="all"),
    dict(name="0 * X^1 = 5", eq="0 * X^1 = 5 * X^0",
         coeffs={0: -5}, degree=0, kind="none"),
    dict(name="constant = 0", eq="5 * X^0 = 0",
         coeffs={0: 5}, degree=0, kind="none"),

    # --- linear
    dict(name="linear, solution 0", eq="4 * X^1 = 0",
         coeffs={1: 4}, degree=1, kind="linear", roots=[0]),
    dict(name="repeated power summed", eq="1 * X^1 + 1 * X^1 = 4 * X^0",
         coeffs={0: -4, 1: 2}, degree=1, kind="linear", roots=[2]),

    # --- too high
    dict(name="degree 4", eq="1 * X^4 = 1 * X^0",
         coeffs={0: -1, 4: 1}, degree=4, kind="high", roots=None),

    # --- stdin (no argument)
    dict(name="STDIN when no argument", eq="5 * X^0 + 4 * X^1 = 4 * X^0", stdin=True,
         coeffs={0: 1, 1: 4}, degree=1, kind="linear", roots=[-0.25]),
]

BONUS = [
    # --- cubics (accepted either as "can't solve" or with correct roots)
    dict(name="cubic: 3 real roots", eq="1 * X^3 - 6 * X^2 + 11 * X^1 - 6 * X^0 = 0",
         coeffs={0: -6, 1: 11, 2: -6, 3: 1}, degree=3, kind="high", roots=[1, 2, 3]),
    dict(name="cubic: 1 real + 2 complex", eq="1 * X^3 + 6 * X^1 - 20 * X^0 = 0",
         coeffs={0: -20, 1: 6, 3: 1}, degree=3, kind="high", roots=[2, -1 + 3j, -1 - 3j]),
    dict(name="cubic: double root", eq="1 * X^3 - 3 * X^1 + 2 * X^0 = 0",
         coeffs={0: 2, 1: -3, 3: 1}, degree=3, kind="high", roots=[-2, 1]),
    dict(name="cubic: a != 1, b != 0", eq="2 * X^3 - 4 * X^2 - 22 * X^1 + 24 * X^0 = 0",
         coeffs={0: 24, 1: -22, 2: -4, 3: 2}, degree=3, kind="high", roots=[1, 4, -3]),

    # --- bad input: must not crash and must not claim a solution
    dict(name="bad input: empty", eq="", kind="error"),
    dict(name="bad input: garbage", eq="hello world", kind="error"),
    dict(name="bad input: missing right side", eq="5 * X^0 =", kind="error"),
    dict(name="bad input: no '='", eq="5 * X^0 + 2 * X^1", kind="error"),
    dict(name="bad input: two '='", eq="5 * X^0 = 3 * X^0 = 2 * X^0", kind="error"),
    dict(name="bad input: double operator", eq="5 * X^0 + + 2 * X^1 = 0", kind="error"),
    dict(name="bad input: negative exponent", eq="5 * X^-1 = 0", kind="error"),

    # --- strict syntax: free-form notation must be rejected, not auto-parsed
    dict(name="bad input: free form (implicit X, missing X^0)", eq="5 + 4 * X + X^2= X^2", kind="error"),
    dict(name="bad input: free form (implicit coefficient/exponent)", eq="X^2 - 4 = 0", kind="error"),
    dict(name="bad input: free form (implicit X^1, no '*')", eq="-X = 3", kind="error"),
    dict(name="bad input: free form (implicit X^0/X^1)", eq="X^2 = X", kind="error"),
]

KEYWORDS = {
    "pos":  r"positive",
    "zero": r"\bzero\b|\bnull\b|\bnul\b|equal to 0|\bis 0\b",
    "neg":  r"negative",
    "all":  r"any real|all real|every real|infinit|all numbers|each real",
    "none": r"no solution|no real solution",
    "high": r"can'?t solve|cannot solve|greater than [23]|unable",
}
ROOT_COUNT = {"pos": (2, 2), "neg": (2, 2), "zero": (1, 2), "linear": (1, 1), "high": (1, 3)}


# ───────────────────────── runner ─────────────────────────

def run(cmd, t):
    args = list(cmd) if t.get("stdin") else list(cmd) + [t["eq"]]
    stdin = t["eq"] + "\n" if t.get("stdin") else ""
    try:
        p = subprocess.run(args, input=stdin, capture_output=True, text=True, timeout=TIMEOUT)
        return p.stdout, p.stderr, p.returncode, None
    except subprocess.TimeoutExpired:
        return "", "", None, f"timeout after {TIMEOUT}s"
    except FileNotFoundError:
        sys.exit(f"cannot run: {' '.join(cmd)}")


def check(t, out, err, rc):
    fails, warns = [], []
    crashed = "Traceback" in err or "Segmentation" in err or (rc is not None and rc < 0)

    if t["kind"] == "error":
        if crashed:
            fails.append("crashed on bad input")
        if not (out.strip() or err.strip()):
            fails.append("printed nothing (should report the error)")
        if re.search(r"solution", out, re.I) and not re.search(r"no solution", out, re.I):
            fails.append("claims a solution for invalid input")
        return fails, warns

    if crashed or rc not in (0, None):
        fails.append(f"crashed (exit {rc})" + (": " + err.strip().splitlines()[-1] if err.strip() else ""))
        return fails, warns

    lines = out.splitlines()
    body = "\n".join(l for l in lines if not l.lower().startswith("reduced"))

    # reduced form
    red = next((l for l in lines if l.lower().startswith("reduced form")), None)
    if red is None:
        fails.append("no 'Reduced form:' line")
    else:
        expr = red.split(":", 1)[1]
        if expr.count("=") != 1:
            fails.append(f"reduced form not 'expr = 0': {expr.strip()}")
        else:
            left, right = expr.split("=")
            got, rhs = parse_poly(left), parse_poly(right)
            if got is None or rhs is None:
                fails.append(f"can't parse reduced form: {expr.strip()}")
            elif any(abs(v) > 1e-12 for v in rhs.values()):
                fails.append(f"right side of reduced form is not 0: {right.strip()}")
            else:
                exp = {k: float(v) for k, v in t["coeffs"].items()}
                bad = [p for p in set(got) | set(exp) if not close(got.get(p, 0.0), exp.get(p, 0.0))]
                if bad:
                    fails.append(f"reduced form wrong for X^{sorted(bad)}: {expr.strip()}")

    # degree
    m = re.search(r"degree\s*:?\s*(\d+)", out, re.I)
    if t["degree"] == 0:
        if m and int(m.group(1)) != 0:
            fails.append(f"degree is {m.group(1)}, expected 0")
    elif not m:
        fails.append("no 'Polynomial degree:' line")
    elif int(m.group(1)) != t["degree"]:
        fails.append(f"degree is {m.group(1)}, expected {t['degree']}")

    # message keyword
    kind = t["kind"]
    has_kw = kind in KEYWORDS and re.search(KEYWORDS[kind], body, re.I)

    # roots
    got = [z for z in (parse_number(l) for l in body.splitlines()) if z is not None]
    expected = t.get("roots")

    if kind == "high":
        solved = expected is not None and got and covers(got, expected)
        if not has_kw and not solved:
            if got and expected is not None:
                fails.append(f"wrong cubic roots: got {[fmt(z) for z in got]}, "
                             f"expected {[fmt(z) for z in expected]}")
            else:
                fails.append("expected a 'can't solve' message (degree > 2)")
    else:
        if kind in KEYWORDS and not has_kw:
            fails.append(f"missing message matching /{KEYWORDS[kind]}/")
        if kind in ROOT_COUNT:
            lo, hi = ROOT_COUNT[kind]
            if not (lo <= len(got) <= hi):
                fails.append(f"expected {lo if lo == hi else f'{lo}-{hi}'} solution(s), "
                             f"found {len(got)}: {[fmt(z) for z in got]}")
            elif not covers(got, expected):
                fails.append(f"wrong solutions: got {[fmt(z) for z in got]}, "
                             f"expected {[fmt(z) for z in expected]}")
        if kind in ("all", "none") and got:
            warns.append(f"prints numbers as solutions: {[fmt(z) for z in got]}")

    if re.search(r"(?<![\w.])-0(?:\.0+)?(?![\d.eE])", body):
        warns.append("prints '-0' (use x + 0.0 or normalise)")

    return fails, warns


def main():
    args = sys.argv[1:]
    verbose = "-v" in args
    bonus = "--bonus" in args
    cmd = [a for a in args if a not in ("-v", "--bonus")] or ["./main"]

    suites = [("MANDATORY", MANDATORY)] + ([("BONUS", BONUS)] if bonus else [])
    total = passed = 0
    for title, tests in suites:
        print(BOLD(f"\n── {title} ──"))
        for t in tests:
            total += 1
            out, err, rc, timeout = run(cmd, t)
            fails, warns = ([timeout], []) if timeout else check(t, out, err, rc)
            label = f'{t["name"]:<38} {repr(t["eq"]) if not t.get("stdin") else "(stdin) " + repr(t["eq"])}'
            if fails:
                print(f"{KO('FAIL')} {label}")
                for f in fails:
                    print(f"       ✗ {f}")
            else:
                passed += 1
                print(f"{OK('PASS')} {label}")
            for w in warns:
                print(f"       {WARN('! ' + w)}")
            if (fails or verbose) and (out or err):
                shown = (out + (("\n[stderr]\n" + err) if err.strip() else "")).rstrip().splitlines()
                for l in shown[:12]:
                    print(f"       │ {l}")
                if len(shown) > 12:
                    print(f"       │ ... ({len(shown) - 12} more lines)")

    color = OK if passed == total else KO
    print(BOLD(f"\n{color(f'{passed}/{total}')} tests passed"))
    sys.exit(0 if passed == total else 1)


if __name__ == "__main__":
    main()