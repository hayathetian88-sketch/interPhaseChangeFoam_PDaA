const volScalarField limitedAlpha1
(
    min(max(alpha1, scalar(0)), scalar(1))
);
tmp<volScalarField> tRb = mixture->rRb(limitedAlpha1);
const volScalarField& rRb = tRb();

const dimensionedScalar& nu = mixture->nuw(); 
const dimensionedScalar& kd = mixture->kd(); 
volVectorField gradp = fvc::grad(p_rgh);

// Vs を必ず作る
volVectorField Vs(
    IOobject(
        "Vs",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimVelocity, vector::zero)
);

// Vs の計算
Vs = ((-16.0) / (3.0 * sqr(rRb) * (rho1 * nu.value() * kd.value()))) * gradp;

// Vs の境界条件
forAll(Vs.boundaryField(), patchI)
{
    Vs.boundaryFieldRef()[patchI] = U.boundaryField()[patchI];
}

Vs.correctBoundaryConditions();
Vs.write();

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

taus = (((rho1 * rho2) / rho) * (scalar(1) - limitedAlpha1) * limitedAlpha1) * sqr(Vs);
taus.write();
