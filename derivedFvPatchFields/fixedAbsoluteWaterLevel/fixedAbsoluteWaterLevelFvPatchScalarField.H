/*---------------------------------------------------------------------------*\
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

Class
    Foam::fixedAbsoluteWaterLevelFvPatchScalarField

Description

    Sets a fixed absolute water level H + S on a patch
    
    Example of the boundary condition specification:
    @verbatim
    inlet
    {
        type            fixedAbsoluteWaterLevel;
        fixedH          164.2;     // absolute water level [m]
        value           uniform 0; // placeholder
    }
    @endverbatim

SourceFiles
    fixedAbsoluteWaterLevelFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedAbsoluteWaterLevelFvPatchScalarFields_H
#define fixedAbsoluteWaterLevelFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fixedAbsoluteWaterLevelFvPatch Declaration
\*---------------------------------------------------------------------------*/

class fixedAbsoluteWaterLevelFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- critical flowdepth at inlet
        scalar fixedH_;

public:

    //- Runtime type information
    TypeName("fixedAbsoluteWaterLevel");


    // Constructors

        //- Construct from patch and internal field
        fixedAbsoluteWaterLevelFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedAbsoluteWaterLevelFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedAbsoluteWaterLevelFvPatchScalarField onto
        //  a new patch
        fixedAbsoluteWaterLevelFvPatchScalarField
        (
            const fixedAbsoluteWaterLevelFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedAbsoluteWaterLevelFvPatchScalarField
        (
            const fixedAbsoluteWaterLevelFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedAbsoluteWaterLevelFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedAbsoluteWaterLevelFvPatchScalarField
        (
            const fixedAbsoluteWaterLevelFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedAbsoluteWaterLevelFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the critical flowdepth
            scalar fixedH() const
            {
                return fixedH_;
            }

            //- Return reference to the critical flowdepth to allow adjustment
            scalar& fixedH()
            {
                return fixedH_;
            }
        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

        //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
