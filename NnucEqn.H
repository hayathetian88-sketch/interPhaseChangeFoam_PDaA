Vv = U + (rho1)/(rho)*alpha1*Vs;
surfaceScalarField phiVv = fvc::interpolate(Vv) & mesh.Sf();

fvScalarMatrix NnucEqn
(
    fvm::ddt(NnucN0)
  + fv::gaussConvectionScheme<scalar>
        (
          mesh,
          phiVv,
          upwind<scalar>(mesh, phiVv)
        ).fvmDiv(phiVv, NnucN0)
  - fvm::Sp(fvc::div(phi), NnucN0)
);

NnucEqn.relax();
NnucEqn.solve();

Info<< "NnucN0 magnitude min: " 
    << min(mag(NnucN0)).value() << " , max: "
    << max(mag(NnucN0)).value() << nl << endl;