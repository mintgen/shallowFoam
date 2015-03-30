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
    dervied from flowRateInletVelocityFvPatchVectorField.C
\*---------------------------------------------------------------------------*/

#include "energyHeadInletFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
energyHeadInletFvPatchVectorField::
energyHeadInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    energyHead_(0),
    HName_("H")
{}


Foam::
energyHeadInletFvPatchVectorField::
energyHeadInletFvPatchVectorField
(
    const energyHeadInletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    energyHead_(ptf.energyHead_),
    HName_(ptf.HName_)
{}


Foam::
energyHeadInletFvPatchVectorField::
energyHeadInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    energyHead_(readScalar(dict.lookup("energyHead"))),
    HName_("H")
{

    if (dict.found("H"))
    {
        dict.lookup("H") >> HName_;
    }
}


Foam::
energyHeadInletFvPatchVectorField::
energyHeadInletFvPatchVectorField
(
    const energyHeadInletFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    energyHead_(ptf.energyHead_),
    HName_(ptf.HName_)
{}


Foam::
energyHeadInletFvPatchVectorField::
energyHeadInletFvPatchVectorField
(
    const energyHeadInletFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    energyHead_(ptf.energyHead_),
    HName_(ptf.HName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::energyHeadInletFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    vectorField n = patch().nf();

    const volScalarField& H = db().lookupObject<volScalarField>(HName_);

    scalar surfaceSumH = gSum(H.boundaryField()[patch().index()]*patch().magSf());
//    Info << "surfaceSumH: " << surfaceSumH << endl;

    scalar Hmean = surfaceSumH/gSum(patch().magSf());
//        Info << "Hmean: " << Hmean << endl;
    scalar Umean = sqrt(2*9.81*(mag(energyHead_-Hmean)));

//    Info << "HU " << -sign(energyHead_-Hmean)*n*Hmean*Umean << endl;

    operator==(-sign(energyHead_-Hmean)*n*Hmean*Umean);

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::energyHeadInletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);

    os.writeKeyword("energyHead") << energyHead_
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
       energyHeadInletFvPatchVectorField
   );
}


// ************************************************************************* //
