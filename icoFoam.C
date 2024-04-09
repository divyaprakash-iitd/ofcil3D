/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    icoFoam

Group
    grpIncompressibleSolvers

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

    \heading Solver details
    The solver uses the PISO algorithm to solve the continuity equation:

        \f[
            \div \vec{U} = 0
        \f]

    and momentum equation:

        \f[
            \ddt{\vec{U}}
          + \div \left( \vec{U} \vec{U} \right)
          - \div \left(\nu \grad \vec{U} \right)
          = - \grad p
        \f]

    Where:
    \vartable
        \vec{U} | Velocity
        p       | Pressure
    \endvartable

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"
// My inclusions
#include "sphereToCell.H"
#include <set>
#include "ListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "diracdelta.H"
// Declare the external subroutine
extern "C" {
    void generateellipse(int *noelpts);
    // void sayhello();
    void getpositions(double *pposx, double *pposy, double *pposz, int *noelpts);
    // subroutine getpositions(XC,YC,ZC) bind(C)
    void calculateforces(double *pfx, double *pfy, double *pfz, int *noelpts);
    // subroutine calculateforces(FXC,FYC,FZC) bind(C)
    void updatepositions(double *pvx, double *pvy, double *pvz, double *dt, int *noelpts);
    // subroutine updatepositions(U,V,W,dt) bind(C)
    void arraycheck(double *pxyz, int* n);
    // void arraycheck(int *pxyz, int* n);
}

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, laminar flow"
        " of Newtonian fluids."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    pisoControl piso(mesh);

    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    int noelpts;
    generateellipse(&noelpts);
    Info << "No. of ellipse points: " << noelpts << endl;

    // Initial position of the point source
    scalar pzmid = (mesh.C()[1][0] - mesh.C()[0][0])/2.0;
    Info << "ZZ: " << pzmid << endl;
    vector rr(0.0,0.0,pzmid); // Define the point

    // // Create vectors for storing the particle's Positions, Forces and Velocity
    std::vector<double> pposx(noelpts,0), pposy(noelpts,0), pposz(noelpts,0);   //Position
    std::vector<double> pfx(noelpts,0), pfy(noelpts,0), pfz(noelpts,0);         // Force
    std::vector<double> pvx(noelpts,0), pvy(noelpts,0), pvz(noelpts,0);         // Velocity

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        // Foam::sphereToCell sph(mesh,rr,0.005,0);

        // Info << "SPH: " << sph.typeName() << endl;

	    //# define omega 0.05
	    //const dimensionedVector mySource("mySource", dimensionSet(0,1,-2,0,0,0,0), 1000*Foam::sin(runTime.value()*omega)*vector(0,1,0));

        ///////////////////////////////////////////////////////////////////////////////////////////        
        // Get positions of the ellipse
        getpositions(pposx.data(),pposy.data(),pposz.data(),&noelpts); 
        // Calculate the forces in the particle
        calculateforces(pfx.data(),pfy.data(),pfz.data(),&noelpts);

        // Create a list of lists to store the neighbours for 
        // every node of the particle
        Foam::DynamicList<Foam::labelList> neighborsList;
        for (int inoelpts = 0; inoelpts < noelpts; ++inoelpts) {
            // Get the node's location
            rr[0] = pposx[inoelpts]; rr[1] = pposy[inoelpts]; rr[2] = pposz[inoelpts];
            #include "getNeighbours.H"
            neighborsList.append(appended);
            appended.clear();
        }


        // Initialize the force source term to zero
        F = F*0;
        // Go through each of the nodes and then their neighbours
        for (int inoelpts = 0; inoelpts < noelpts; ++inoelpts) {
            // Get the node's location
            rr[0] = pposx[inoelpts]; rr[1] = pposy[inoelpts]; rr[2] = pposz[inoelpts];
            // Get the node's force
            vector pf(0.0,0.0,0.0);
            pf[0] = pfx[inoelpts]; pf[1] = pfy[inoelpts]; pf[2] = pfz[inoelpts];
            // Get the neighbors list for that node alone
            labelList singleNodeNeighbors = neighborsList[inoelpts];

            // Iterate through each of the neighbouring cells of the selected node 
            forAll(singleNodeNeighbors,idx) {
                int icell = singleNodeNeighbors[idx];   
                #include "interpolateForces.H"
            }
        }
        ///////////////////////////////////////////////////////////////////////////////////////////

        // Momentum predictor
        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U) 
          - F // IBM Source Term 
        );

        if (piso.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));
        }

        // --- PISO loop
        while (piso.correct())
        {
            volScalarField rAU(1.0/UEqn.A());
            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                fvc::flux(HbyA)
              + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
            );

            adjustPhi(phiHbyA, U, p);

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(p, U, phiHbyA, rAU);

            // Non-orthogonal pressure corrector loop
            while (piso.correctNonOrthogonal())
            {
                // Pressure corrector

                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);

                pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

                if (piso.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
        }

        ///////////////////////////////////////////////////////////////////////////////////////////
        // Go through each of the nodes and then their neighbours
        for (int inoelpts = 0; inoelpts < noelpts; ++inoelpts) {
            // Get the node's location
            rr[0] = pposx[inoelpts];
            rr[1] = pposy[inoelpts];
            rr[2] = pposz[inoelpts];
            // Get the neighbors list for that node alone
            labelList singleNodeNeighbors = neighborsList[inoelpts];

            // Define a variable to store the point source velocity
            vector pu(0,0,0);
            // Iterate through each of the neighbouring cells of the selected node 
            forAll(singleNodeNeighbors,idx) {
                int icell = singleNodeNeighbors[idx];   
                #include "interpolateVelocity.H"
            }
            // Store the interpolated velocity at each node
            pvx[inoelpts] = pu[0]; 
            pvy[inoelpts] = pu[1];
            pvz[inoelpts] = pu[2];
        }
        // Clear the list of lists of neighbours
        neighborsList.clear();
        double simdt = runTime.deltaTValue();

        // Update the position of all the nodes
        updatepositions(pvx.data(),pvy.data(),pvz.data(),&simdt,&noelpts);
        ///////////////////////////////////////////////////////////////////////////////////////////

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
