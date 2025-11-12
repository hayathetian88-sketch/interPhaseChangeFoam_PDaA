const volScalarField limitedAlpha1
(
    min(max(alpha1, scalar(0)), scalar(1))
);
tmp<volScalarField> tRb = mixture->rRb(limitedAlpha1);
const volScalarField& rRb = tRb();

const dimensionedScalar& nu = mixture->nuw(); 
const dimensionedScalar& kd = mixture->kd(); 

volVectorField gradp = fvc::grad(p_rgh);


if (!mesh.foundObject<volVectorField>("Vs"))
{
    Info<< "Creating volVectorField Vs" << endl;
    volVectorField Vs
    (
        IOobject
        (
            "Vs",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimVelocity, vector::zero)
    );

    // Vs difinition
    Vs = ((-16.0) / (3.0 * sqr(rRb) * (rho1 * nu * kd))) * gradp;

    // Vs boundary condition
    forAll(Vs.boundaryField(), patchI)
    {
        Vs.boundaryFieldRef()[patchI] = U.boundaryField()[patchI];
    }

    Vs.correctBoundaryConditions();
    Vs.write();
}
else
{
    volVectorField& Vs = const_cast<volVectorField&>(
        mesh.lookupObject<volVectorField>("Vs")
    );

    // recalculate
    Vs = ((-16.0) / (3.0 * sqr(rRb) * (rho1 * nu * kd))) * gradp;

    forAll(Vs.boundaryField(), patchI)
    {
        Vs.boundaryFieldRef()[patchI] = U.boundaryField()[patchI];
    }

    Vs.correctBoundaryConditions();
    Vs.write();
}

