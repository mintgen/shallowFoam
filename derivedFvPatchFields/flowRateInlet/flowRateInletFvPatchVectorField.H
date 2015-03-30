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
    Foam::flowRateInletFvPatchVectorField

Description

    Sets HU such that following holds:

          HU = flowRate / widthInflowPatch

    Example of the boundary condition specification:
    @verbatim
    inlet
    {
        type            flowRateInlet;
        flowRate        50;              // Inflow rate [m^3/s]
        value           uniform (0 0 0); // placeholder
    }
    @endverbatim

Author
    KM-Turbulenz GmbH, 2009
    dervied from flowRateInletVelocityFvPatchVectorField.H

SourceFiles
    flowRateInletFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef flowRateInletFvPatchVectorField_H
#define flowRateInletFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
               Class flowRateInletFvPatch Declaration
\*---------------------------------------------------------------------------*/

class flowRateInletFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Inlet EnergieHoehe
        scalar flowRate_;

    //
    word HName_;


public:

   //- Runtime type information
   TypeName("flowRateInlet");


   // Constructors

        //- Construct from patch and internal field
        flowRateInletFvPatchVectorField
        (
	    const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        flowRateInletFvPatchVectorField
        (
	    const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  flowRateInletFvPatchVectorField
        //  onto a new patch
        flowRateInletFvPatchVectorField
        (
            const flowRateInletFvPatchVectorField&,
	    const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        flowRateInletFvPatchVectorField
        (
            const flowRateInletFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new flowRateInletFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        flowRateInletFvPatchVectorField
        (
            const flowRateInletFvPatchVectorField&,
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
                new flowRateInletFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the flux
            scalar flowRate() const
            {
                return flowRate_;
            }

            //- Return reference to the flux to allow adjustment
            scalar& flowRate()
            {
                return flowRate_;
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
