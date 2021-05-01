#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Ep \n" << endl;
volScalarField Ep
(
    IOobject
    (
        "Ep",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

Info<< "Reading transportProperties\n" << endl;

immiscibleConstitutiveTwoPhaseMixture mixture(U,phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

IOdictionary constitutiveDict
(
    IOobject
    (
        "constitutiveProperties",    // dictionary name
        runTime.constant(),     // dict is found in "constant"
        mesh,                   // registry for the dict
        IOobject::MUST_READ,    // must exist, otherwise failure
        IOobject::NO_WRITE      // dict is only read by the solver
    )
);

dimensionedScalar eps1_("eps", dimensionSet(-1, -3, 4, 0, 0, 2, 0), scalar(0.0));
dimensionedScalar eps2_("eps", dimensionSet(-1, -3, 4, 0, 0, 2, 0), scalar(0.0));

dimensionedScalar sig1_("sig", dimensionSet(-1, -3, 3, 0, 0, 2, 0), scalar(0.0));
dimensionedScalar sig2_("sig", dimensionSet(-1, -3, 3, 0, 0, 2, 0), scalar(0.0));


wordList phaseNames(constitutiveDict.lookup("phases"));

dictionary& subDict1 = constitutiveDict.subDict(phaseNames[0]).subDict("parameters");

eps1_ = subDict1.lookup("eps");
sig1_ = subDict1.lookup("sig");

dictionary& subDict2 = constitutiveDict.subDict(phaseNames[1]).subDict("parameters");

eps2_ = subDict2.lookup("eps");
sig2_ = subDict2.lookup("sig");

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();
const dimensionedScalar& eps1 = eps1_; 
const dimensionedScalar& eps2 = eps2_;
const dimensionedScalar& sig1 = sig1_;
const dimensionedScalar& sig2 = sig2_;

// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();

// ADD
volScalarField eps
(
    IOobject
    (
        "eps",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alpha1*eps1 + alpha2*eps2,
//        scalar(1) / (alpha1/eps1 + (scalar(1) - alpha1)/eps2),
        alpha1.boundaryField().types()
);
    eps.oldTime();

// ADD
volScalarField sig
(
    IOobject
    (
        "sig",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alpha1*sig1 + alpha2*sig2,
//        scalar(1) / (alpha1/sig1 + (scalar(1) - alpha1)/sig2),
        alpha1.boundaryField().types()
);
    sig.oldTime();

// ADD
volScalarField rhoE  //Electric Charge density
(
    IOobject
    (
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
//    fvc::laplacian(eps,Ep)
        mesh
);

volTensorField MST
(
    IOobject
    (
        "MST",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    eps*(fvc::grad(Ep)*fvc::grad(Ep) - scalar(0.5)*magSqr(fvc::grad(Ep))*I)
);

// ADD
volVectorField Ef  //Electric Field
(
    IOobject
    (
        "Ef",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(Ep)
);

// ADD
/*
volVectorField Efp  //Dielectrophoretic Force
(
    IOobject
    (
        "Efp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -scalar(0.5)*magSqr(Ef)*fvc::grad(eps)
);

volVectorField Efc  //Electrophoretic Force
(
    IOobject
    (
        "Efc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::laplacian(eps,Ep)*fvc::grad(Ep)
);

volVectorField Eft //Total Electric Force
(
    IOobject
    (
        "Eft",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -scalar(0.5)*magSqr(Ef)*fvc::grad(eps)-fvc::laplacian(eps,Ep)*fvc::grad(Ep)
);
*/

// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

#include "createMRF.H"
#include "createFvOptions.H"
