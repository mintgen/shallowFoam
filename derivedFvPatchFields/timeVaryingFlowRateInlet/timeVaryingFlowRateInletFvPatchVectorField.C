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

Author
    KM-Turbulenz GmbH, 2009

\*---------------------------------------------------------------------------*/

#include "timeVaryingFlowRateInletFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
timeVaryingFlowRateInletFvPatchVectorField::
timeVaryingFlowRateInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    timeSeries_(),
    HName_("H")
{}


Foam::
timeVaryingFlowRateInletFvPatchVectorField::
timeVaryingFlowRateInletFvPatchVectorField
(
    const timeVaryingFlowRateInletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    timeSeries_(ptf.timeSeries_),
    HName_(ptf.HName_)
{}


Foam::
timeVaryingFlowRateInletFvPatchVectorField::
timeVaryingFlowRateInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    timeSeries_(dict),
    HName_("H")
{
   if (dict.found("value"))
   {
       fvPatchField<vector>::operator==(Field<vector>("value", dict, p.size()));
   }
   else
   {
       updateCoeffs();
   }

    if (dict.found("H"))
    {
        dict.lookup("H") >> HName_;
    }
}


Foam::
timeVaryingFlowRateInletFvPatchVectorField::
timeVaryingFlowRateInletFvPatchVectorField
(
    const timeVaryingFlowRateInletFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    timeSeries_(ptf.timeSeries_),
    HName_(ptf.HName_)
{}


Foam::
timeVaryingFlowRateInletFvPatchVectorField::
timeVaryingFlowRateInletFvPatchVectorField
(
    const timeVaryingFlowRateInletFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    timeSeries_(ptf.timeSeries_),
    HName_(ptf.HName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeVaryingFlowRateInletFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }


    vectorField n = patch().nf();

    const volScalarField& H = db().lookupObject<volScalarField>(HName_);

    scalar flowRate = timeSeries_(this->db().time().timeOutputValue());

    scalar surfaceSumH = gSum(H.boundaryField()[patch().index()]*patch().magSf());
    scalar Umean = flowRate/surfaceSumH;
    scalarField phip = patch().lookupPatchField<surfaceScalarField, scalar>("phi");
    scalar surfaceSumPhip = gSum(phip*H.boundaryField()[patch().index()]);
    Info << "Inflow rate - Q: " <<  surfaceSumPhip << endl;   

    operator==(-n*Umean*H.boundaryField()[patch().index()]);

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::timeVaryingFlowRateInletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);

    os.writeKeyword("timeSeries") << timeSeries_
        << token::END_STATEMENT << nl;

    if (HName_ != "H")
    {
        os.writeKeyword("H") << HName_ << token::END_STATEMENT << nl;
    }

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       timeVaryingFlowRateInletFvPatchVectorField
   );
}


// ************************************************************************* //
