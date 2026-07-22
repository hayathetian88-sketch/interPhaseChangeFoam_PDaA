//surfaceScalarField phiVs = fvc::interpolate(Vs) & mesh.Sf();
Vv = U + (rho1)/(rho)*alpha1*Vs;
surfaceScalarField phiVv = fvc::interpolate(Vv) & mesh.Sf();
// surfaceScalarField phiNnuc = phi + phiVs;

// fvScalarMatrix NnucEqn
// (
//     fv::EulerDdtScheme<scalar>(mesh).fvmDdt(NnucN0)
//     + fv::gaussConvectionScheme<scalar>
//       (
//         mesh,
//         phiNnuc,
//         upwind<scalar>(mesh, phiNnuc)
//       ).fvmDiv(phiNnuc, NnucN0)
// );

//surfaceScalarField PhiVv =fvc::interpolate(alpha1)*(phi + phiVs);

fvScalarMatrix NnucEqn
(
    fvm::ddt(NnucN0)
  + fv::gaussConvectionScheme<scalar>
        (
          mesh,
          phiVv,
          upwind<scalar>(mesh, phiVv)
        ).fvmDiv(phiVv, NnucN0)
);

NnucEqn.relax();
NnucEqn.solve();

Info<< "NnucN0 magnitude min: " 
    << min(mag(NnucN0)).value() << " , max: "
    << max(mag(NnucN0)).value() << nl << endl;