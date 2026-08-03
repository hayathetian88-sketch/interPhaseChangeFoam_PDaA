Vv = U + (rho1)/(rho)*alpha1*Vs;
surfaceScalarField phiVv = fvc::interpolate(Vv) & mesh.Sf();
surfaceScalarField phiVvj = fvc::interpolate(alpha1*Vs) & mesh.Sf();
dimensionedScalar aSmall("aSmall", dimless, VSMALL);

fvScalarMatrix lnucEqn
(
    fvm::ddt(lnuc)
  + fv::gaussConvectionScheme<scalar>
    (
        mesh,
        phiVv,
        upwind<scalar>(mesh, phiVv)
    ).fvmDiv(phiVv, lnuc)
  + fvm::SuSp(-fvc::div(phiVv), lnuc)
  ==
  - fvc::div(phiVvj)/ max(alpha1, aSmall)
);

lnucEqn.relax();
lnucEqn.solve();
Info<< "lnuc magnitude min: " 
    << min(mag(lnuc)).value() << " , max: "
    << max(mag(lnuc)).value() << nl << endl;

dimensionedScalar Nnuc_dim("Nnuc_dim", dimless/dimVolume, 1.0);
// Nnucの更新時にも次元を掛ける
Nnuc = exp(lnuc) * Nnuc_dim;
Info<< "Nnuc magnitude min: " 
    << min(mag(Nnuc)).value() << " , max: "
    << max(mag(Nnuc)).value() << nl << endl;