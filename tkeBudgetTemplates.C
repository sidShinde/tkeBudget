/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "volFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::tkeBudget::processField
(
    const word& fieldName,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvalue
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> FieldType;

    const word scopedName = word("tke_") + fieldName;

    if (obr_.foundObject<FieldType>(scopedName))
    {
        FieldType& fld =
            const_cast<FieldType&>(obr_.lookupObject<FieldType>(scopedName));
        fld == tvalue();
    }
    else if (obr_.found(scopedName))
    {
        WarningInFunction
            << "Cannot store turbulence field " << scopedName
            << " since an object with that name already exists"
            << nl << endl;
    }
    else
    {
      obr_.store
        (
            new FieldType
            (
                IOobject
                (
                    scopedName,
                    obr_.time().timeName(),
                    obr_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                tvalue
            )
        );
    }
}


// returns the convection term
Foam::tmp<Foam::volScalarField>
Foam::functionObjects::tkeBudget::Ck() const
{
  const volVectorField& U = obr_.lookupObject<volVectorField>("U");
  //const volVectorField& UMean = obr_.lookupObject<volVectorField>("UMean");
  //const volVectorField UPrime = U - UMean;

  // trace of Reynolds stress tensor
  const volSymmTensorField& UP2M = obr_.lookupObject<volSymmTensorField>("UPrime2Mean");

  // inner product of uprime
  //const volScalarField ipUPrime( UPrime & UPrime );
  const volScalarField ipUPrime( tr(UP2M) );
  
  // gradient of tIPUPrime
  const volVectorField gradIPUPrime(fvc::grad(ipUPrime));

  // inner product of U and tgradIPUPrime
  const volScalarField ipUUPrime( U & gradIPUPrime );

  return tmp<volScalarField>
  (
      new volScalarField
      (
          IOobject
	  (
	      "Ck",
	      ipUUPrime.mesh().time().timeName(),
	      ipUUPrime.mesh()
	   ),
	  (-0.5)*ipUUPrime,
	  ipUUPrime.boundaryField().types()
       )
   );

}


// returns the production term
Foam::tmp<Foam::volScalarField>
Foam::functionObjects::tkeBudget::Pk() const
{
  const volVectorField& U = obr_.lookupObject<volVectorField>("U");
  const volSymmTensorField& UP2M = obr_.lookupObject<volSymmTensorField>("UPrime2Mean");

  // gradient of velocity
  const volTensorField gradU( fvc::grad(U) );

  // inner product of UP2M and gradU
  const volScalarField ipUP2MGradU( UP2M && gradU );

  return tmp<volScalarField>
  (
      new volScalarField
      (
          IOobject
	  (
	      "Pk",
	      ipUP2MGradU.mesh().time().timeName(),
	      ipUP2MGradU.mesh()
	   ),
	  (-1.0)*ipUP2MGradU,
	  ipUP2MGradU.boundaryField().types()
       )
   );
}


// returns the turbulence transport term
Foam::tmp<Foam::volScalarField>
Foam::functionObjects::tkeBudget::Tk() const
{
  const volVectorField& U = obr_.lookupObject<volVectorField>("U");
  const volVectorField& UMean = obr_.lookupObject<volVectorField>("UMean");
  
  // fluctutation field
  const volVectorField UPrime = U - UMean;

  // inner product of UPrime and UPrime
  const volScalarField ipUPrime( UPrime & UPrime );

  // outer product of fluctuation field
  const volVectorField opUPrime( ipUPrime * UPrime );

  // divergence of opUPrime
  const volScalarField divOpUPrime( fvc::div(opUPrime) );

  return tmp<volScalarField>
  (
      new volScalarField
      (
          IOobject
	  (
	      "Tk",
	      divOpUPrime.mesh().time().timeName(),
	      divOpUPrime.mesh()
	   ),
	  (-0.5)*divOpUPrime,
	  divOpUPrime.boundaryField().types()
       )
   );

}


// returns the viscous diffusion term
Foam::tmp<Foam::volScalarField>
Foam::functionObjects::tkeBudget::Dk() const
{
  const volSymmTensorField& UP2M = obr_.lookupObject<volSymmTensorField>("UPrime2Mean");

  // trace of UP2M
  const volScalarField trUP2M( tr(UP2M) );

  // laplacian of trUP2M
  const volScalarField lapUP2M( fvc::laplacian(trUP2M) );

  return tmp<volScalarField>
  (
      new volScalarField
      (
          IOobject
	  (
	      "Dk",
	      lapUP2M.mesh().time().timeName(),
	      lapUP2M.mesh()
	   ),
	  (1/Re)*(0.5)*lapUP2M,
	  lapUP2M.boundaryField().types()
       )
   );
}


// returns the viscous dissipation term
Foam::tmp<Foam::volScalarField>
Foam::functionObjects::tkeBudget::Epik() const
{
  const volVectorField& U = obr_.lookupObject<volVectorField>("U");
  const volVectorField& UMean = obr_.lookupObject<volVectorField>("UMean");
 
  // fluctuation field
  const volVectorField UPrime = U - UMean;

  // gradient of the fluctuation field
  const volTensorField gradUPrime( fvc::grad(UPrime) );

  // double inner product of gradUPrime
  const volScalarField ipGradUPrime( gradUPrime && gradUPrime );
  
  return tmp<volScalarField>
  (
      new volScalarField
      (
          IOobject
	  (
	      "Epik",
	      ipGradUPrime.mesh().time().timeName(),
	      ipGradUPrime.mesh()
	   ),
	  (-1/Re)*ipGradUPrime,
	  ipGradUPrime.boundaryField().types()
       )
   );
}


// returns the velocity-pressure gradient term
Foam::tmp<Foam::volScalarField>
Foam::functionObjects::tkeBudget::Pik() const
{
  const volVectorField& U = obr_.lookupObject<volVectorField>("U");
  const volVectorField& UMean = obr_.lookupObject<volVectorField>("UMean");
 
  const volScalarField& P = obr_.lookupObject<volScalarField>("p");
  const volScalarField& PMean = obr_.lookupObject<volScalarField>("pMean");

  // fluctuation field
  const volVectorField UPrime = U - UMean;
  const volScalarField PPrime = P - PMean;

  // gradient of the PPrime field
  const volVectorField gradPPrime( fvc::grad(PPrime) );

  // double inner product of UPrime and gradPPrime
  const volScalarField ipUPrimeGradPPrime( UPrime & gradPPrime );
  
  return tmp<volScalarField>
  (
      new volScalarField
      (
          IOobject
	  (
	      "Pik",
	      ipUPrimeGradPPrime.mesh().time().timeName(),
	      ipUPrimeGradPPrime.mesh()
	   ),
	  (-1.0)*ipUPrimeGradPPrime,
	  ipUPrimeGradPPrime.boundaryField().types()
       )
   );
}

// ************************************************************************* //
