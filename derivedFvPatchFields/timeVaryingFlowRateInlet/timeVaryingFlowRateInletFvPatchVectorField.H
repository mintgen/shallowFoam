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

    Transient version of BC flowRateInlet. Needs a data file in the case folder that contains flow rate [m^3/s] over time [s] (see example below)

    Example of the boundary condition specification:
    @verbatim
    inlet
    {
        type            timeVaryingFlowRateInlet;
        fileName        "timeSeries";              // data file in case folder
        outOfBounds     clamp;                     // see OpenFOAM documentation
        value           uniform (0 0 0);           // placeholder
    }
    @endverbatim

    Example of the data file:
    @verbatim
    (
    (0 0)
    (300 1.75)
    (600 2.5)
    (900 1.5)
    (1200 0)
    )
    @endverbatim

Author
    KM-Turbulenz GmbH, 2009

SourceFiles
    timeVaryingFlowRateInletFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingFlowRateInletFvPatchVectorField_H
#define timeVaryingFlowRateInletFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
               Class timeVaryingFlowRateInletFvPatch Declaration
\*---------------------------------------------------------------------------*/

class timeVaryingFlowRateInletFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        interpolationTable<scalar> timeSeries_;

    //
    word HName_;


public:

   //- Runtime type information
   TypeName("timeVaryingFlowRateInlet");


   // Constructors

        //- Construct from patch and internal field
        timeVaryingFlowRateInletFvPatchVectorField
        (
	    const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeVaryingFlowRateInletFvPatchVectorField
        (
	    const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  flowRateInletFvPatchVectorField
        //  onto a new patch
        timeVaryingFlowRateInletFvPatchVectorField
        (
            const timeVaryingFlowRateInletFvPatchVectorField&,
	    const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        timeVaryingFlowRateInletFvPatchVectorField
        (
            const timeVaryingFlowRateInletFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new timeVaryingFlowRateInletFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        timeVaryingFlowRateInletFvPatchVectorField
        (
            const timeVaryingFlowRateInletFvPatchVectorField&,
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
                new timeVaryingFlowRateInletFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the flux
              const interpolationTable<scalar>& timeSeries() const
            {
                return timeSeries_;
            }

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
