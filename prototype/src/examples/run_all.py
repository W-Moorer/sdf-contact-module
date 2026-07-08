import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.examples import (case01_fixed_body, case02_revolute_pendulum,
                          case03_driven_revolute, case04_four_bar,
                          case06_falling_block, case07_ring_torsion,
                          case08_rmd_import_be, case09_validation)


def main():
    cases = [
        ("Case 01: Fixed body (static)", case01_fixed_body.run, 5),
        ("Case 02: Revolute pendulum", case02_revolute_pendulum.run, 30),
        ("Case 03: Driven revolute", case03_driven_revolute.run, 30),
        ("Case 04: Double pendulum", case04_four_bar.run, 60),
        ("Case 06: Falling block contact", case06_falling_block.run, 30),
        ("Case 07: Ring-on-cube torsion", case07_ring_torsion.run, 90),
        ("Case 08: RMD import + implicit BE", case08_rmd_import_be.run, 60),
        ("Case 09: Validation suite", case09_validation.run, 120),
    ]

    passed, total = 0, len(cases)
    for name, func, timeout in cases:
        print(f"{'='*60}")
        print(f"  {name}")
        print(f"{'='*60}")
        t0 = time.time()
        try:
            func()
            passed += 1
            print(f"  Time: {time.time()-t0:.2f}s")
        except Exception as e:
            print(f"  FAILED after {time.time()-t0:.1f}s: {e}")
            import traceback; traceback.print_exc()
        print()

    print(f"{'='*60}")
    print(f"  Result: {passed}/{total} cases passed")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
