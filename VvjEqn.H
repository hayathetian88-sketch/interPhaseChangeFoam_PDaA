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

// Vs calculate
Vvj = Fyr * limitedAlpha1 * ((-16.0) / (3.0 * sqr(rRb) * (rho1 * nu * kd))) * gradp;

// Vs boundary condition
forAll(Vvj.boundaryField(), patchI)
{
    Vvj.boundaryFieldRef()[patchI] = U.boundaryField()[patchI];
}

Vvj.correctBoundaryConditions();

Info<< "Vvj magnitude min: " 
    << min(mag(Vvj)).value() << " , max: "
    << max(mag(Vvj)).value() << nl << endl;


taus = (((rho1 * rho2) / rho) * (scalar(1) - limitedAlpha1) / limitedAlpha1) * sqr(Vvj);

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
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

slipcoeff = (rho2) / rho;

volScalarField pcoeff
(
    IOobject
    (
        "pcoeff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(-1, 3, 1, 0, 0, 0, 0), 0.0)
);

pcoeff = ((rho1 - rho2) / rho) * Fyr * ((-16.0) / (3.0 * sqr(rRb) * (rho1 * nu * kd)));

Info<< "VvjEqn end" << endl;