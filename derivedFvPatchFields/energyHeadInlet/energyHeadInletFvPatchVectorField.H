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
    Foam::energyHeadInletFvPatchVectorField

Description

    Sets HU such that following holds:
          energyHead = H + U^2/2g

    Example of the boundary condition specification:
    @verbatim
    inlet
    {
        type            energyHeadInlet;
        energieHead        0.2;        // energy head [m]
        value           uniform (0 0 0); // placeholder
    }
    @endverbatim

Note
    - H at the boundary must be specified (z.B. zeroGradient)
    - Inflow normal to boundary

Author
    KM-Turbulenz GmbH, 2009

SourceFiles
    energyHeadInletFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef energyHeadInletFvPatchVectorField_H
#define energyHeadInletFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
               Class energyHeadInletFvPatch Declaration
\*---------------------------------------------------------------------------*/

class energyHeadInletFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Inlet energy head
        scalar energyHead_;

    //
    word HName_;


public:

   //- Runtime type information
   TypeName("energyHeadInlet");


   // Constructors

        //- Construct from patch and internal field
        energyHeadInletFvPatchVectorField
        (
	    const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        energyHeadInletFvPatchVectorField
        (
	    const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  energyHeadInletFvPatchVectorField
        //  onto a new patch
        energyHeadInletFvPatchVectorField
        (
            const energyHeadInletFvPatchVectorField&,
	    const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        energyHeadInletFvPatchVectorField
        (
            const energyHeadInletFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new energyHeadInletFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        energyHeadInletFvPatchVectorField
        (
            const energyHeadInletFvPatchVectorField&,
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
                new energyHeadInletFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the energy head
            scalar energyHead() const
            {
                return energyHead_;
            }

            //- Return reference to the energy head to allow adjustment
            scalar& energyHead()
            {
                return energyHead_;
            }


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
