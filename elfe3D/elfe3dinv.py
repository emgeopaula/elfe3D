import pygimli as pg
import pyelfe3d
import numpy as np
#import scipy.sparse as sp

class elfe3dJacobian(pg.Matrix):
    def __init__(self, Jrows, Jcols, Jvec, JTvec, cached_model):
        super().__init__()
        self._Jrows = Jrows
        self._Jcols = Jcols
        self._Jvec = Jvec
        self._JTvec = JTvec
        self._cached_model = cached_model
    
    def rows(self):
        return len(self._Jrows)

    def cols(self):
        return len(self._Jcols)

    def mult(self, x):
        """Multiply the Jacobian with a vector, Jm."""
        print('[JACOB-mult]')
        print('[JACOB-mult, x]', x)

        # Set Fortran model to x 
        x_arr = np.array(x)
        pyelfe3d.elfe3d.inv_model = x_arr
        
        # Call solver
        pyelfe3d.elfe3d.solve()
        
        # Reset Fortran model to value from previous iteration
        pyelfe3d.elfe3d.inv_model = self._cached_model
        
        return pyelfe3d.elfe3d.Jvec

    def transMult(self, x):
        """Multiply  Jacobian transposed with a vector, Jᵀd = (dJᵀ)ᵀ."""
        print('[JACOB-transMult]')
        print('[JACOB-transMult, x]', x)

        # Set Fortran test_data_vec to x
        x_arr = np.array(x)
        pyelfe3d.elfe3d.test_data_vec = x_arr
        
        # Call solver
        pyelfe3d.elfe3d.solve()
        
        # Reset Fortran model to value from previous iteration
        pyelfe3d.elfe3d.inv_model = self._cached_model

        # Reset test_data_vec to default value
        # (see compute_pseudo_fwd in mod_sensitivities.f90)
        pyelfe3d.elfe3d.test_data_vec[:] = 77.0
        
        return pyelfe3d.elfe3d.JTvec


class elfe3DModelling(pg.Modelling):
    """Use elfe3D solver as a forward operator within pyGIMLI inversion framework."""

    def __init__(self, verbose=True):
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


    def _run_elfe3d_solver(self, model) -> None:
        """
        Executes the elfe3d solver and updates self._cached_* variables afterwards
        
        Parameters:
            - model: The current model proposed by the inversion framework 
        """

        if model is not None:
            m_arr = np.array(model)
            pyelfe3d.elfe3d.inv_model = m_arr
            print('[changed-model]')
            pyelfe3d.elfe3d.solve()
        else:
            # First iteration: let elfe3d guess the initial model
            print('[none-model]')
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
        print('[start-model]')
        self._run_elfe3d_solver(model=None)
        print('[observed_data]', self._observed_data)
        print('[errors]', self._errors)
        return self._cached_model


    def response(self, model):
        """
        Calculates forward response for a given model.

        Parameters:
            - model: The current model proposed by the inversion framework

        Returns:
            - response: forward_data
        """
        print('[response]')
        self._run_elfe3d_solver(model)
        return self._cached_response 


    def createJacobian(self, model):
        """Instantiates and sets matrix-free Jacobian operator for pyGIMLi."""
        print('[JACOB]')
        self._jacobian = elfe3dJacobian(
            self._cached_jrows,
            self._cached_jcols,
            self._cached_jvec,
            self._cached_jtvec,
            self._cached_model)
        self.setJacobian(self._jacobian)


def main():

    # Define errors:
    errors = pg.Vector([
        3.01376403E-07,
        2.22861987E-07,
        1.62550486E-07,
        5.16305677E-11,
        3.94458122E-11,
        3.06961268E-11
    ])

    # Define observed data:
    observed_data = pg.Vector([
        3.01376403E-06,
        2.22861987E-06,
        1.62550486E-06,
        5.16305677E-10,
        3.94458122E-10,
        3.06961268E-10
    ]) 
    
    # Define fast forward operator
    fop = elfe3DModelling()

    # Use the fop (elfe3d) to define the start model
    #start_model = fop.createStartModel()

    inv = pg.Inversion(fop=fop)
    
    # Run inversion
    inv_model = inv.run(dataVals=observed_data, errorVals=errors, maxIter=3)

    # When done, print model history
    print(['model-history'])
    for i, model in enumerate(inv.modelHistory):
        print(f"Iteration {i}: Model = {model}")



if __name__ == "__main__":
    main()
