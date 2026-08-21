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

#include "createPhi.H"


Info<< "Creating phaseChangeTwoPhaseMixture\n" << endl;
autoPtr<phaseChangeTwoPhaseMixture> mixture =
    phaseChangeTwoPhaseMixture::New(U, phi);

volScalarField& alpha1(mixture->alpha1());
volScalarField& alpha2(mixture->alpha2());

const dimensionedScalar& rho1 = mixture->rho1();
const dimensionedScalar& rho2 = mixture->rho2();


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


// Construct interface from alpha1 distribution
interfaceProperties interface(alpha1, U, mixture());

// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture())
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

Info<< "Creating Wall Distance Field (yn)\n" << endl;

volScalarField yn
(
    IOobject
    (
        "yn",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("yn", dimLength, 0.0)
);

if (yn.headerOk())
{
    Info<< "Reading wall distance field yn\n" << endl;
}
else
{
    Info<< "Calculating wall distance field yn\n" << endl;
    
    wallDist y(mesh);
    yn = y.y();
}

// DEBUG
volScalarField slipTerm
(
    IOobject
    (
        "slipTerm",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE // ★ これで時間ディレクトリに書き出されます
    ),
    mesh,
    dimensionedScalar("slipTerm", dimless/dimTime, 0.0) // 次元は [1/s]
);


volScalarField vDotvmcAlphal
(
    IOobject
    (
        "vDotvmcAlphal",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE // ★ これで時間ディレクトリに書き出されます
    ),
    mesh,
    dimensionedScalar("vDotvmcAlphal", dimless/dimTime, 0.0) // 次元は [1/s]
);
volScalarField TotalSourceTerm
(
    IOobject
    (
        "TotalSourceTerm",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE // ★ これで時間ディレクトリに書き出されます
    ),
    mesh,
    dimensionedScalar("TotalSourceTerm", dimless/dimTime, 0.0) // 次元は [1/s]
);
volScalarField blockFactor
(
    IOobject("blockFactor", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("one", dimless, 1.0)
);
//

const dimensionSet dimStress(dimMass/dimLength/sqr(dimTime));

volSymmTensorField taus
(
    IOobject
    (
        "taus",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedSymmTensor("taus", dimStress, symmTensor::zero)
);

const dimensionSet revdimlength(dimless/dimLength);

volScalarField rRb
(
    IOobject
    (
        "rRb",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", revdimlength, 0.0)
);


volVectorField Vs
(
    IOobject
    (
        "Vs",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField Vv
(
    IOobject
    (
        "Vv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimVelocity, vector::zero)
);

const dimensionSet revdimtime(dimless/dimTime);

volScalarField time_coeff
(
    IOobject
    (
        "time_coeff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", revdimtime, 0.0)
);

volVectorField Vt
(
    IOobject
    (
        "Vt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimVelocity, vector::zero)
);

// volScalarField slipdot
// (
//     IOobject
//     (
//         "slipdot",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::NO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("zero", revdimtime, 0.0)
// );

volVectorField pcoeff
(
    IOobject
    (
        "pcoeff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimVelocity, vector::zero)
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

#include "createFvOptions.H"

IOobject alphaPhi10Header
(
    IOobject::groupName("alphaPhi0", alpha1.group()),
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
);

// MULES flux from previous time-step
surfaceScalarField alphaPhi10
(
    alphaPhi10Header,
    phi*fvc::interpolate(alpha1)
);
