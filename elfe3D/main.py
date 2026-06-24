import numpy as np
import pyelfe3d

def main():

    print("Launching Fortran elfe3d solver...")

    pyelfe3d.elfe3d.solve() # This executes your Fortran 'subroutine solve'

    print("Solver finished successfully!")

    # Grab the arrays straight out of the module object
    Jvec  = pyelfe3d.elfe3d.Jvec
    JTvec = pyelfe3d.elfe3d.JTvec
    Jrows = pyelfe3d.elfe3d.Jrows
    Jcols = pyelfe3d.elfe3d.Jcols
    forward_data = pyelfe3d.elfe3d.forward_data

    print(f"Jvec shape: {Jvec.shape}")
    print(f"JTvec shape: {JTvec.shape}")
    print(f"Jrows shape: {Jrows.shape}")
    print(f"Jcols shape: {Jcols.shape}")
    print(f"forward_data shape: {forward_data.shape}")

    # Check each condition individually
    jvec_ok = np.allclose(Jvec, 666.0)
    jtvec_ok = np.allclose(JTvec, 444.0)
    jrows_ok = np.all(Jrows == 6)
    jcols_ok = np.all(Jcols == 4)
    forward_data_ok = np.allclose(forward_data, 999.0)

    print(f"Jvec all 666.0: {jvec_ok}")
    print(f"JTvec all 444.0: {jtvec_ok}")
    print(f"Jrows all 6: {jrows_ok}")
    print(f"Jcols all 4: {jcols_ok}")
    print(f"forward_data all 999.0: {forward_data_ok}")

    if (jvec_ok and jtvec_ok and jrows_ok and jcols_ok and forward_data_ok):
        print("All tests passed!")
    else:
        print("Some tests failed!")

if __name__ == "__main__":
    main()
