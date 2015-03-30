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

#include "dischargeFlowDepthOutletFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
dischargeFlowDepthOutletFvPatchVectorField::
dischargeFlowDepthOutletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    HName_("H")
{}


Foam::
dischargeFlowDepthOutletFvPatchVectorField::
dischargeFlowDepthOutletFvPatchVectorField
(
    const dischargeFlowDepthOutletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    HName_(ptf.HName_)
{}


Foam::
dischargeFlowDepthOutletFvPatchVectorField::
dischargeFlowDepthOutletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    HName_("H")
{

    if (dict.found("H"))
    {
        dict.lookup("H") >> HName_;
    }
}


Foam::
dischargeFlowDepthOutletFvPatchVectorField::
dischargeFlowDepthOutletFvPatchVectorField
(
    const dischargeFlowDepthOutletFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    HName_(ptf.HName_)
{}


Foam::
dischargeFlowDepthOutletFvPatchVectorField::
dischargeFlowDepthOutletFvPatchVectorField
(
    const dischargeFlowDepthOutletFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    HName_(ptf.HName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dischargeFlowDepthOutletFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const vectorField n = patch().nf();
    const fvPatchScalarField Sp = patch().lookupPatchField<volScalarField, scalar>("S");
    const scalarField gradSp = Sp.snGrad();
    const scalarField kstp = patch().lookupPatchField<volScalarField, scalar>("kst");
    const scalarField Hp = patch().lookupPatchField<volScalarField, scalar>(HName_);

    const scalarField phip = patch().lookupPatchField<surfaceScalarField, scalar>("phi");
    const scalar surfaceSumPhip = gSum(phip*Hp);

    operator==(n*sqrt(mag(gradSp))*kstp*pow(Hp,(5.0/3.0)));

    Info << "Outflow rate - Q: " <<  surfaceSumPhip << endl;

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::dischargeFlowDepthOutletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);

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
       dischargeFlowDepthOutletFvPatchVectorField
   );
}


// ************************************************************************* //
