volScalarField limitedAlpha1
(
    IOobject
    (
        "limitedAlpha1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,    
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("limitedAlpha1Value", dimless, 0.0)
);
limitedAlpha1 = min(max(alpha1, scalar(1e-12)), scalar(1));

tmp<volScalarField> tRb = mixture->rRb(limitedAlpha1);
rRb = tRb();
Info<< "rRb magnitude min: " 
    << min(rRb).value() << " , max: "
    << max(rRb).value() << nl << endl;

const dimensionedScalar& nu = mixture->nuw(); 
const dimensionedScalar& kd = mixture->kd();
const dimensionedScalar& kw = mixture->kw();

volVectorField gradpRaw = fvc::grad(p_rgh);
volVectorField gradpSmooth = fvc::average(fvc::interpolate(gradpRaw));
//volVectorField gradp = fvc::grad(p_rgh);
volScalarField Re_b = (2.0 / rRb) * mag(Vs) / nu;
volScalarField kd_dynamic = 24.0 * (1.0 + 0.15 * pow(max(Re_b, scalar(1e-5)), 0.687))*pow(limitedAlpha1,-0.5);

time_coeff = (3.0/8.0)*(kd_dynamic*nu)*(rRb*rRb);
Info<< "time_coeff magnitude min: " 
    << min(time_coeff).value() << " , max: "
    << max(time_coeff).value() << nl << endl;
Vt= -1*(16.0/3.0)*(1/(kd_dynamic*rho1*nu))*(1/(rRb*rRb))*gradpSmooth;
dimensionedScalar maxVt("maxVt", dimVelocity, 50.0);
dimensionedScalar vSmall("vSmall", dimVelocity, VSMALL);

Vt = Vt * min(maxVt / (mag(Vt) + vSmall), dimensionedScalar("one", dimless, 1.0));
Info<< "Vt magnitude min: " 
    << min(mag(Vt)).value() << " , max: "
    << max(mag(Vt)).value() << nl << endl;

surfaceScalarField phiG = phi + fvc::flux(Vs);
volTensorField gradUSmooth = fvc::average(fvc::interpolate(fvc::grad(U)));
fvVectorMatrix VsEqn
(
    fvm::ddt(Vs)
  + fvm::div(phiG, Vs) + fvm::SuSp(-fvc::div(phiG), Vs)
  + (Vs & gradUSmooth)                              
 ==
  - fvm::Sp(time_coeff, Vs)
  + time_coeff * Vt
);

VsEqn.relax();
VsEqn.solve();

Info<< "Vs magnitude min: " 
    << min(mag(Vs)).value() << " , max: "
    << max(mag(Vs)).value() << nl << endl;

Vgj = (1 - limitedAlpha1) * Vs ;

slipcoeff = (rho2) / rho;

pcoeff = ((rho1 - rho2) / rho) * Vs;

Info<< "VvjEqn end" << endl;