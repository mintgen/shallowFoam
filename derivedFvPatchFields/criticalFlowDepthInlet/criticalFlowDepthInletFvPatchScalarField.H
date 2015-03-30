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
    Foam::criticalFlowDepthInletFvPatchScalarField

Description
    
    In order to avoid super critical conditions at the inlet, this BC limits the absolute flow depth H + S to a user defined critical absolute flow depth criticalH.

    Example of the boundary condition specification:
    @verbatim
    inlet
    {
        type            criticalFlowDepthInlet;
        criticalH       164.2;        // critical _absolute_ flow depth [m]
        value           uniform 0; // placeholder
    }
    @endverbatim

SourceFiles
    criticalFlowDepthInletFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef criticalFlowDepthInletFvPatchScalarFields_H
#define criticalFlowDepthInletFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class criticalFlowDepthInletFvPatch Declaration
\*---------------------------------------------------------------------------*/

class criticalFlowDepthInletFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- critical flowdepth at inlet
        scalar criticalH_;

public:

    //- Runtime type information
    TypeName("criticalFlowDepthInlet");


    // Constructors

        //- Construct from patch and internal field
        criticalFlowDepthInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        criticalFlowDepthInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given criticalFlowDepthInletFvPatchScalarField onto
        //  a new patch
        criticalFlowDepthInletFvPatchScalarField
        (
            const criticalFlowDepthInletFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        criticalFlowDepthInletFvPatchScalarField
        (
            const criticalFlowDepthInletFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new criticalFlowDepthInletFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        criticalFlowDepthInletFvPatchScalarField
        (
            const criticalFlowDepthInletFvPatchScalarField&,
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
                new criticalFlowDepthInletFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the critical flowdepth
            scalar criticalH() const
            {
                return criticalH_;
            }

            //- Return reference to the critical flowdepth to allow adjustment
            scalar& criticalH()
            {
                return criticalH_;
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
