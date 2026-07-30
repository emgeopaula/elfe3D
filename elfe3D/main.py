import numpy as np
import pygimli
import pyelfe3d

def test_jacobian(Jvec, JTvec, Jrows, Jcols) -> None:
    assert Jvec.shape == (6,), f"Jvec shape: expected (6,), got {Jvec.shape}"
    assert JTvec.shape == (40,), f"JTvec shape: expected (40,), got {JTvec.shape}"
    assert Jrows.shape == (6,), f"Jrows shape: expected (6,), got {Jrows.shape}"
    assert Jcols.shape == (40,), f"Jcols shape: expected (40,), got {Jcols.shape}"

    # Value checks for Jrows (1..6) and Jcols (1..40)
    assert np.array_equal(Jrows, np.arange(1, 7)), f"Jrows values mismatch: got {Jrows}"
    assert np.array_equal(Jcols, np.arange(1, 41)), f"Jcols values mismatch: got {Jcols}"

    Jvec_expected = np.array([
        -9.6063503026389437E-009,
        -1.0023161150621066E-008,
        -1.0459021573926215E-008,
        -6.2032370033720292E-011,
        -6.4146007720678829E-011,
        -6.6452234674365855E-011
    ])
    assert np.allclose(Jvec, Jvec_expected)

    JTvec_expected = np.array([
        -1.0287583090042686E-010, -2.2051249959716073E-012, -3.9104659123266598E-011,
        -1.2146399434867539E-011, -8.9838761919397479E-012, -8.6492132482740330E-012,
        -3.9297915655294007E-011, -7.4127375808313418E-012, -3.0816741696238783E-011,
        -1.3142557566667370E-010, -6.9864590914746287E-012, -6.6382125543174224E-011,
        -5.9294384272321313E-012, -1.1855369834650099E-010, -5.8724322221009524E-012,
        -5.9965269807425682E-011, -4.3311237247384389E-011, -6.3002352971932053E-011,
        -1.2092527696776448E-010, -7.5048437528987134E-012, -1.3913681974172897E-011,
        -1.1304033935342531E-011, -3.1142732227444033E-011, -7.4760133283074636E-012,
        -1.5264106874454601E-011, -2.2755008520981667E-011, -3.4062929031111071E-012,
        -2.3809115503825035E-011, -6.2889517271978067E-012, -3.5476162936231666E-012,
        -3.1467298430821801E-012, -9.9893511107386989E-012, -6.1578575761608676E-012,
        -4.2133754432473104E-012, -3.2144189603364880E-011, -3.5026061750504144E-011,
        -2.2892161400620746E-011, -7.5016700256356401E-011, -7.0978984772833998E-011,
        -2.6904212349915077E-011
    ])
    assert np.allclose(JTvec, JTvec_expected)

def test_forward_data(forward_data) -> None:
    assert forward_data.shape == (6,), f"forward_data shape: expected (6,), got {forward_data.shape}"
    expected = np.array([
        4.44562188E-05,
        2.80089204E-05,
        1.97479819E-05,
        3.22833081E-09,
        2.76825052E-09,
        2.51720481E-09
    ])
    assert np.allclose(forward_data, expected)

def test_inv_model(inv_model) -> None:
    assert inv_model.shape == (40,), f"inv_model shape: expected (40,), got {inv_model.shape}"
    assert np.array_equal(inv_model, np.ones(40)), f"inv_model values mismatch: got {inv_model}"

def test_elfe3d_output() -> None:
    # Check that elfe3d produces the same output on consecutive calls
    # This checks for correctness on the python side
    for i in range(3):
        pyelfe3d.elfe3d.solve()
        Jvec  = pyelfe3d.elfe3d.Jvec
        JTvec = pyelfe3d.elfe3d.JTvec
        Jrows = pyelfe3d.elfe3d.Jrows
        Jcols = pyelfe3d.elfe3d.Jcols
        forward_data = pyelfe3d.elfe3d.forward_data
        inv_model = pyelfe3d.elfe3d.inv_model
        test_jacobian(Jvec, JTvec, Jrows, Jcols)
        test_forward_data(forward_data)
        test_inv_model(inv_model)


def main():

    print("PyGimli version:", pygimli.__version__)

    # Check that consecutive calls to elfe3d produce the same (expected) results
    test_elfe3d_output()

    # Call elfe3d
    print("Launching Fortran elfe3d solver...")
    #pyelfe3d.elfe3d.solve() # This executes your Fortran 'subroutine solve'

    # Modify the values of inv_model
    #inv_model = pyelfe3d.elfe3d.inv_model

    #print("Solver finished successfully!")

if __name__ == "__main__":
    main()
