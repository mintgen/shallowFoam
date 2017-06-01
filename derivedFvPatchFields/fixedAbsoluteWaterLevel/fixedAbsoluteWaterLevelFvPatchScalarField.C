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

\*---------------------------------------------------------------------------*/

#include "fixedAbsoluteWaterLevelFvPatchScalarField.H"
#include "freestreamFvPatchFields.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedAbsoluteWaterLevelFvPatchScalarField::fixedAbsoluteWaterLevelFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    fixedH_(0)
{}


fixedAbsoluteWaterLevelFvPatchScalarField::fixedAbsoluteWaterLevelFvPatchScalarField
(
    const fixedAbsoluteWaterLevelFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    fixedH_(ptf.fixedH_)    
{}


fixedAbsoluteWaterLevelFvPatchScalarField::fixedAbsoluteWaterLevelFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    fixedH_(readScalar(dict.lookup("fixedH")))
{}


fixedAbsoluteWaterLevelFvPatchScalarField::fixedAbsoluteWaterLevelFvPatchScalarField
(
    const fixedAbsoluteWaterLevelFvPatchScalarField& wbppsf
)
:
    fixedValueFvPatchField<scalar>(wbppsf),
    fixedH_(wbppsf.fixedH_)
{}


fixedAbsoluteWaterLevelFvPatchScalarField::fixedAbsoluteWaterLevelFvPatchScalarField
(
    const fixedAbsoluteWaterLevelFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(wbppsf, iF),
    fixedH_(wbppsf.fixedH_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fixedAbsoluteWaterLevelFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
 
    // Get bottom elevation at patch
    scalarField Sp = patch().lookupPatchField<volScalarField, scalar>("S");
    
    // operator==(max(fixedH_-Sp,0));
    operator==((fixedH_-Sp)*pos(fixedH_-Sp));

    fixedValueFvPatchScalarField::updateCoeffs();
}

void fixedAbsoluteWaterLevelFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);

    os.writeKeyword("fixedH") << fixedH_
        << token::END_STATEMENT << nl;

    writeEntry("value", os);
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, fixedAbsoluteWaterLevelFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
