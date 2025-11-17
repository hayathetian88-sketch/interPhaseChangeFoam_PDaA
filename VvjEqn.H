const volScalarField limitedAlpha1
(
    min(max(alpha1, scalar(1e-12)), scalar(1))
);
tmp<volScalarField> tRb = mixture->rRb(limitedAlpha1);
const volScalarField& rRb = tRb();

const dimensionedScalar& nu = mixture->nuw(); 
const dimensionedScalar& kd = mixture->kd();
const dimensionedScalar& kw = mixture->kw();
volVectorField gradp = fvc::grad(p_rgh);
const volScalarField yr = pow(kw,-1.0*yn*rRb);

Fyr = pow4(1 - yr);

volVectorField Vvj(
    IOobject(
        "Vvj",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimVelocity, vector::zero)
);

// Vs calculate
Vvj = Fyr * limitedAlpha1 * ((-16.0) / (3.0 * sqr(rRb) * (rho1 * nu.value() * kd.value()))) * gradp;

// Vs boundary condition
forAll(Vvj.boundaryField(), patchI)
{
    Vvj.boundaryFieldRef()[patchI] = U.boundaryField()[patchI];
}

Vvj.correctBoundaryConditions();

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

taus = (((rho1 * rho2) / rho) * (scalar(1) - limitedAlpha1) / limitedAlpha1) * sqr(Vvj);

volVectorField pPlus
(
    IOobject
    (
        "pPlus",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("pPlus", dimVelocity, vector::zero)
);

pPlus = (((rho1 - rho2) / rho) * (scalar(1) - limitedAlpha1)) * Vvj;

volVectorField Vs(
    IOobject(
        "Vs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimVelocity, vector::zero)
);

Vs = Vvj / limitedAlpha1;

volScalarField slipcoeff
(
    IOobject
    (
        "slipcoeff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (rho1 * rho2) / rho
);