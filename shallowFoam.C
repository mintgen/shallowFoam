/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    shallowFoam

Description
    Transient solver for depth averaged shallow water equations.

Author
    KM-Turbulenz GmbH, 2009
    Florian Mintgen, 2012

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"

#   include "readTimeControls.H"
#   include "CourantNo.H"
#   include "setInitialDeltaT.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

///////////////////////////////////////////////////////////////////////////////
///     Transport equations for H und HU
	fvScalarMatrix HEqn
        (
	    fvm::ddt(H)
	    + fvm::div(phi,H)
        );

        HEqn.solve();

        // Hclip = max (H, Hdry)
        Hclip = (H-Hdry)*pos(H-Hdry) + Hdry;

	// Bottom friction based on Strickler-equation
        alpha = (g*dim_s*dim_s/dim_m)*mag(U)/pow(kst,2.0)/pow(Hclip/dim_m,1.0/3.0)/Hclip;

        // nu_t = C_nu * u_tau * H
        nut = Cnu*sqrt(alpha*mag(HU))*H;
        nut = nut*pos(nutmax-nut) + nutmax*pos(nut-nutmax);
        nut.correctBoundaryConditions();

	// Zero-Gradient for faces at wet-dry interface
	#include "zeroGradientWetDry.H"

         fvVectorMatrix HUEqn
         (
             fvm::ddt(HU)          // d(HU_i)/dt                    Local derivative
	     + fvm::div(phi, HU)   // + d/dx_j ( U_j * HU_i )       Convection
	     + H*gradgSpHf          // + H * d/dx_i ( gSpH )        Downhill-slope force and acceleration due to change in flowdepth
	     + fvm::Sp(alpha,HU)   // + alpha * HU                  Bottom friction
	     - fvm::laplacian(nut,HU) // - d^2/dx_j^2 ( nut * HU )  Turbulent stresses
         );  

	HUEqn.solve();

///////////////////////////////////////////////////////////////////////////////
///     Calculation of U and phi

        // if Hclip < Hdry2, then U = 0
        HU = HU*pos(Hclip - Hdry2); 
        U = HU/Hclip;
       
        phi = (fvc::interpolate(U) & mesh.Sf());

///////////////////////////////////////////////////////////////////////////////

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
