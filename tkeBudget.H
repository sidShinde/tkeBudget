/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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

Class
    Foam::functionObjects::turbulenceFields

Group
    grpFieldFunctionObjects

Description
    This function object stores turbulence fields on the mesh database for
    further manipulation.

    Fields are stored as copies of the original, with the prefix
    "tubulenceModel:", e.g.:

    \verbatim
    turbulenceModel:R
    \endverbatim

    Example of function object specification:
    \verbatim
    turbulenceFields1
    {
        type        turbulenceFields;
        libs ("libfieldFunctionObjects.so");
        ...
        fields
        (
            R
            devRhoReff
        );
    }
    \endverbatim

Usage
    \table
        Property     | Description                 | Required | Default value
        type         | type name: processorField   | yes      |
        fields       | fields to store (see below) | yes      |
    \endtable

    Where \c fields can include:
    \plaintable
        k           | turbulence kinetic energy
        epsilon     | turbulence kinetic energy dissipation rate
        omega       | turbulence specific dissipation rate
        nut         | turbulence viscosity (incompressible)
        nuEff       | effective turbulence viscosity (incompressible)
        mut         | turbulence viscosity (compressible)
        muEff       | effective turbulence viscosity (compressible)
        alphat      | turbulence thermal diffusivity (compressible)
        alphaEff    | effective turbulence thermal diffusivity (compressible)
        R           | Reynolds stress tensor
        devReff     | Deviatoric part of the effective Reynolds stress
        devRhoReff  | Divergence of the Reynolds stress
    \endplaintable

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::timeControl

SourceFiles
    tkeBudget.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_tkeBudget_H
#define functionObjects_tkeBudget_H

#include "fvMeshFunctionObject.H"
#include "HashSet.H"
#include "NamedEnum.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                      Class turbulenceFields Declaration
\*---------------------------------------------------------------------------*/

//class turbulenceFields
class tkeBudget
:
    public fvMeshFunctionObject
{
public:
    enum incompressibleField
    {
      ifCk,   // convection term
      ifPk,   // production term
      ifTk,   // turbulence transport
      ifDk,   // viscous diffusion
      ifEpik, // viscous dissipation
      ifPik   // velocity-pressure gradient
    };
    static const NamedEnum<incompressibleField, 6> incompressibleFieldNames_;

protected:

    // Protected data

        //- Fields to load
        wordHashSet fieldSet_;

    // Protected Member Functions

        //- Process the tke field
        template<class Type>
        void processField
        (
            const word& fieldName,
            const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvalue
        );

  // return the convection term Ck
  tmp<volScalarField> Ck() const;

  // return the production term Pk
  tmp<volScalarField> Pk() const;

  // return the production term Tk
  tmp<volScalarField> Tk() const;

  // return the production term Dk
  tmp<volScalarField> Dk() const;

  // return the production term Dk
  tmp<volScalarField> Epik() const;

  // return the production term Dk
  tmp<volScalarField> Pik() const;

private:
  // private variables

    // Private member functions

        //- Disallow default bitwise copy construct
        tkeBudget(const tkeBudget&);

        //- Disallow default bitwise assignment
        void operator=(const tkeBudget&);

public:

    //- Runtime type information
    TypeName("tkeBudget");

    // Constructors

        //- Construct from Time and dictionary
        tkeBudget
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~tkeBudget();

    // Member Functions

        //- Read the controls
        virtual bool read(const dictionary&);

        //- Calculate turbulence fields
        virtual bool execute();

        //- Do nothing.
        //  The turbulence fields are registered and written automatically
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "tkeBudgetTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
