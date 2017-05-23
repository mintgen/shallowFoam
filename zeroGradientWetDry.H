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

\*---------------------------------------------------------------------------*/

wetDry = pos(H-Hdry);

volScalarField gSpH = g*(S+H);
surfaceScalarField gSpHf = fvc::interpolate(gSpH);
volVectorField gradgSpHf = fvc::grad(gSpHf);

const unallocLabelList& owner = mesh.owner();
const unallocLabelList& neighbour = mesh.neighbour();
const vectorField& Sf = mesh.Sf();
const scalarField& V = mesh.V();

const scalarField& gSpHi = gSpH.internalField();
const scalarField& gSpHfi = gSpHf.internalField();
vectorField& gradgSpHfi = gradgSpHf.ref();

forAll(owner, facei)
{
  label P = owner[facei];
  label N = neighbour[facei];

  if (wetDry[P] != wetDry[N])
	  {
	    gradgSpHfi[P] += Sf[facei]*(gSpHi[P] - gSpHfi[facei]) / V[P];
	    gradgSpHfi[N] -= Sf[facei]*(gSpHi[N] - gSpHfi[facei]) / V[N];
	  }
}


