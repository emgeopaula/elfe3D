import pygimli as pg
import pygimli.meshtools as mt
import pyelfe3d
import numpy as np

import os
import sys
import ctypes
from contextlib import contextmanager

def _flush_c_streams():
    """Flush all C / Fortran stdio buffers to prevent leaked output."""
    try:
        # Load C standard library across Unix/Linux systems
        libc = ctypes.CDLL(None)
        libc.fflush(None)  # NULL flushes all open C/Fortran output streams
    except Exception:
        pass

@contextmanager
def silence_fortran():
    # Flush Python and C/Fortran buffers before changing descriptors
    sys.stdout.flush()
    sys.stderr.flush()
    _flush_c_streams()

    # Get OS file descriptors for stdout (1) and stderr (2)
    stdout_fd = sys.stdout.fileno()
    stderr_fd = sys.stderr.fileno()

    # Save copies of original descriptors
    saved_stdout_fd = os.dup(stdout_fd)
    saved_stderr_fd = os.dup(stderr_fd)

    # Open null device
    null_fd = os.open(os.devnull, os.O_WRONLY)

    try:
        # Redirect both stdout and stderr at the OS level
        os.dup2(null_fd, stdout_fd)
        os.dup2(null_fd, stderr_fd)
        yield
    finally:
        # Flush null streams before swapping back
        _flush_c_streams()

        # Restore original descriptors
        os.dup2(saved_stdout_fd, stdout_fd)
        os.dup2(saved_stderr_fd, stderr_fd)

        # Close temporary descriptors
        os.close(null_fd)
        os.close(saved_stdout_fd)
        os.close(saved_stderr_fd)


class elfe3dJacobian(pg.Matrix):
    def __init__(self, Jrows, Jcols, Jvec, JTvec, cached_model, silence_fortran = True):
        super().__init__()
        self._Jrows = Jrows
        self._Jcols = Jcols
        self._Jvec = Jvec
        self._JTvec = JTvec
        self._cached_model = cached_model
        self._silence_fortran = silence_fortran


    def rows(self):
        return len(self._Jrows)


    def cols(self):
        return len(self._Jcols)


    def mult(self, x):
        """Multiply the Jacobian with a vector, Jm."""

        # Set Fortran model to x 
        x_arr = np.array(x)
        pyelfe3d.elfe3d.inv_model = x_arr
        
        # Call solver
        if self._silence_fortran:
            with silence_fortran():
                pyelfe3d.elfe3d.solve()
        else:
            pyelfe3d.elfe3d.solve()
        
        # Reset Fortran model to value from previous iteration
        pyelfe3d.elfe3d.inv_model = self._cached_model
        
        return pyelfe3d.elfe3d.Jvec


    def transMult(self, x):
        """Multiply  Jacobian transposed with a vector, Jᵀd = (dJᵀ)ᵀ."""

        # Set Fortran test_data_vec to x
        x_arr = np.array(x)
        pyelfe3d.elfe3d.test_data_vec = x_arr
        
        # Call solver
        if self._silence_fortran:
            with silence_fortran():
                pyelfe3d.elfe3d.solve()
        else:
            pyelfe3d.elfe3d.solve()
        
        # Reset Fortran model to value from previous iteration
        pyelfe3d.elfe3d.inv_model = self._cached_model

        # Reset test_data_vec to default value
        # (see compute_pseudo_fwd in mod_sensitivities.f90)
        pyelfe3d.elfe3d.test_data_vec[:] = 77.0
        
        return pyelfe3d.elfe3d.JTvec



class elfe3DModelling(pg.Modelling):
    """Use elfe3D solver as a forward operator within pyGIMLI inversion framework."""

    def __init__(self, verbose=True, silence_fortran = True):
       super().__init__(verbose=verbose)

       # Cache variables to avoid duplicate calls during inversion
       self._cached_model = None
       self._cached_response = None #forward_data
       self._cached_jvec = None
       self._cached_jtvec = None
       self._cached_jrows = None
       self._cached_jcols = None
       self._observed_data = None
       self._errors = None
       self._silence_fortran = silence_fortran


    def _run_elfe3d_solver(self, model) -> None:
        """
        Executes the elfe3d solver and updates self._cached_* variables afterwards
        
        Parameters:
            - model: The current model proposed by the inversion framework 
        """

        if model is not None:
            m_arr = np.array(model)
            pyelfe3d.elfe3d.inv_model = m_arr
        
        if self._silence_fortran:
            with silence_fortran():
                pyelfe3d.elfe3d.solve()
        else:
            pyelfe3d.elfe3d.solve()

        # Update cached results
        self._cached_jvec  = pyelfe3d.elfe3d.Jvec
        self._cached_jtvec = pyelfe3d.elfe3d.JTvec
        self._cached_jrows = pyelfe3d.elfe3d.Jrows
        self._cached_jcols = pyelfe3d.elfe3d.Jcols
        self._cached_response = pg.Vector(pyelfe3d.elfe3d.forward_data)
        self._cached_model = pg.Vector(pyelfe3d.elfe3d.inv_model)
        self._observed_data = pg.Vector(pyelfe3d.elfe3d.observed_data)
        self._errors = pg.Vector(pyelfe3d.elfe3d.errors)


    def createStartModel(self, dataVals):
        """
        Define an initial model for the inversion process.
        Triggers initial run of elfe3D forward solver.
        elfe3D will read the start model from the file specified in 'elfe3D_input.txt'
        
        Parameters:
            - dataVals: observed data
        
        pyGIMLi requires dataVals to be an argument to the createStartModel function,
        but elfe3D doesn't need it as it reads the initial model from a file, so dataVals is ignored.

        """

        self._run_elfe3d_solver(model=None)
        return self._cached_model


    def response(self, model):
        """
        Calculates forward response for a given model.

        Parameters:
            - model: The current model proposed by the inversion framework

        Returns:
            - response: forward_data
        """

        self._run_elfe3d_solver(model)
        return self._cached_response 


    def createJacobian(self, model):
        """Instantiates and sets matrix-free Jacobian operator for pyGIMLi."""

        self._jacobian = elfe3dJacobian(
            self._cached_jrows,
            self._cached_jcols,
            self._cached_jvec,
            self._cached_jtvec,
            self._cached_model,
            self._silence_fortran)
        self.setJacobian(self._jacobian)



def main():

    # Define observed data:
    observed_data = pg.Vector([
        3.01376403E-06,
        2.22861987E-06,
        1.62550486E-06,
        5.16305677E-10,
        3.94458122E-10,
        3.06961268E-10
    ]) 

    # Define errors:
    errors = pg.Vector([
        3.01376403E-07,
        2.22861987E-07,
        1.62550486E-07,
        5.16305677E-11,
        3.94458122E-11,
        3.06961268E-11
    ])
    relative_errors = errors / observed_data
    
    # Define fast forward operator
    fop = elfe3DModelling(silence_fortran=True)

    # Initialize inversion with fop
    inv = pg.Inversion(fop=fop, verbose=True)

    # Run inversion
    inv_model = inv.run(dataVals=observed_data, errorVals=relative_errors)

    # When done, print model history
    print(['model-history'])
    for i, model in enumerate(inv.modelHistory):
        print(f"Iteration {i}: Model = {model}")
        
    # Write models into .vtk files
    # reads tetgen-generated mesh with region markers
    mesh = mt.readTetgen("./in/mesh") 

    # region 3 is the inversion region that is updated during iterations
    # we have to verify if pygimli counts the same way as elfe3D
    region3 = [c.id() for c in mesh.cells() if c.marker() == 3]

    # for each model in the inversion history
    for i, model in enumerate(inv.modelHistory):
    
        # check that number of cells in region 3 is equal the number of model entries
        print(len(region3))
        print(len(model))
        assert len(model) == len(region3)
    
        # fill entire model cells with nan
        vtkdata = np.full(mesh.cellCount(), np.nan)
        vtkdata[region3] = model

        # use current iteration model
        mesh["FEM"] = vtkdata

        # write it to .vtk file (one file per iteration)
        mesh.exportVTK(f"iter_{i:03d}.vtk")

if __name__ == "__main__":
    main()
