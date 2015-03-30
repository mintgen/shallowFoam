/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2006-2008 OpenCFD Ltd.
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

Class
    Foam::dischargeFlowDepthOutletFvPatchVectorField

Description

    When H is specified at the outlet, this BC sets HU such that the Strickler equation is satisfied.

    Example of the boundary condition specification:
    @verbatim
    inlet
    {
        type            dischargeFlowDepthOutlet;
        value           uniform (0 0 0); // placeholder
    }
    @endverbatim

    Note: Since this BC makes use of the Strickler equation, the gradient of the bottom level S should be set properly.

Author
    KM-Turbulenz GmbH, 2009
    dervied from flowRateInletVelocityFvPatchVectorField.H

SourceFiles
    dischargeFlowDepthOutletFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef dischargeFlowDepthOutletFvPatchVectorField_H
#define dischargeFlowDepthOutletFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
               Class dischargeFlowDepthOutletFvPatch Declaration
\*---------------------------------------------------------------------------*/

class dischargeFlowDepthOutletFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

    word HName_;

public:

   //- Runtime type information
   TypeName("dischargeFlowDepthOutlet");


   // Constructors

        //- Construct from patch and internal field
        dischargeFlowDepthOutletFvPatchVectorField
        (
	    const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        dischargeFlowDepthOutletFvPatchVectorField
        (
	    const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  dischargeFlowDepthOutletFvPatchVectorField
        //  onto a new patch
        dischargeFlowDepthOutletFvPatchVectorField
        (
            const dischargeFlowDepthOutletFvPatchVectorField&,
	    const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        dischargeFlowDepthOutletFvPatchVectorField
        (
            const dischargeFlowDepthOutletFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new dischargeFlowDepthOutletFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        dischargeFlowDepthOutletFvPatchVectorField
        (
            const dischargeFlowDepthOutletFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new dischargeFlowDepthOutletFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

        //- Update the coefficients associated with the patch field
//        virtual void updateCoeffs();
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
