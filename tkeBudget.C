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

\*---------------------------------------------------------------------------*/

//#include "turbulenceFields.H"
#include "tkeBudget.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(tkeBudget, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        tkeBudget,
        dictionary
    );
}
}

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::tkeBudget::incompressibleField,
    6
>::names[] =
{
  "Ck",   // convection term
  "Pk",   // production term
  "Tk",   // turbulence transport
  "Dk",   // viscous diffusion
  "Epik", // viscous dissipation
  "Pik"   // velocity-pressure gradient
};

const Foam::NamedEnum
<
    Foam::functionObjects::tkeBudget::incompressibleField,
    6
> Foam::functionObjects::tkeBudget::incompressibleFieldNames_;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

 Foam::functionObjects::tkeBudget::tkeBudget
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldSet_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::tkeBudget::~tkeBudget()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::tkeBudget::read(const dictionary& dict)
{
    if (dict.found("field"))
    {
        fieldSet_.insert(word(dict.lookup("field")));
    }
    else
    {
        fieldSet_.insert(wordList(dict.lookup("fields")));
    }

    Info<< type() << " " << name() << ": ";
    if (fieldSet_.size())
    {
        Info<< "storing fields:" << nl;
        forAllConstIter(wordHashSet, fieldSet_, iter)
        {
            Info<< "    " <<  iter.key() << nl;
        }
        Info<< endl;
    }
    else
    {
        Info<< "no fields requested to be stored" << nl << endl;
    }

    return true;
}


bool Foam::functionObjects::tkeBudget::execute()
{
  forAllConstIter(wordHashSet, fieldSet_, iter)
    {
      const word& f = iter.key();
      switch (incompressibleFieldNames_[f])
	{

	case ifCk:
	  {
	    processField<scalar>(f, Ck());
	    break;
	  }
	case ifPk:
	  {
	    processField<scalar>(f, Pk());
	    break;
	  }
	case ifTk:
	  {
	    processField<scalar>(f, Dk());
	    break;
	  }
	case ifDk:
	  {
	    processField<scalar>(f, Dk());
	    break;
	  }
	case ifEpik:
	  {
	    processField<scalar>(f, Epik());
	    break;
	  }
	case ifPik:
	  {
	    processField<scalar>(f, Pik());
	    break;
	  }
	default:
	  {
	    FatalErrorInFunction
	      << "Invalid field selection" << abort(FatalError);
	  }

	}
    }
	
  return true;
}


bool Foam::functionObjects::tkeBudget::write()
{
    forAllConstIter(wordHashSet, fieldSet_, iter)
    {
      //const word fieldName = modelName + ':' + iter.key();
      const word fieldName = word("tke_") + iter.key();  
      writeObject(fieldName);
    }

    return true;
}


// ************************************************************************* //
