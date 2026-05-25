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

const dimensionedScalar& nu = mixture->nuw(); 
const dimensionedScalar& kd = mixture->kd();
const dimensionedScalar& kw = mixture->kw();
volVectorField gradp = fvc::grad(p_rgh);

time_coeff = (3/8)*(kd*nu)*(rRb*rRb);
Vt= (16/3)*(1/(kd*rho1*nu))*(1/(rRb*rRb))*gradp;

Vg = U + Vs;

surfaceScalarField phiG = phi + fvc::flux(Vs);

fvVectorMatrix VsEqn
(
    fvm::ddt(Vs)
  + fvm::div(phiG, Vs) - fvm::Sp(fvc::div(phiG), Vs)
  + (Vs & fvc::grad(U))                              
 ==
  - fvm::Sp(time_coeff, Vs)
  + time_coeff * Vt
);

VsEqn.solve();

Info<< "Vs magnitude min: " 
    << min(mag(Vs)).value() << " , max: "
    << max(mag(Vs)).value() << nl << endl;

Vgj = (1 - limitedAlpha1) * Vs ;

slipcoeff = (rho2) / rho;

pcoeff = ((rho1 - rho2) / rho) * Vs;

Info<< "VvjEqn end" << endl;
