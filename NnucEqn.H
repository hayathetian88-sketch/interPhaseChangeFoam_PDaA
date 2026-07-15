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
    fvm::ddt(aNnucN0)
  + fv::gaussConvectionScheme<scalar>
        (
          mesh,
          phiVv,
          upwind<scalar>(mesh, phiVv)
        ).fvmDiv(phiVv, aNnucN0)
);

NnucEqn.relax();
NnucEqn.solve();

Info<< "sum(aNnuc)="
    << fvc::domainIntegrate(aNnucN0)
    << endl;

NnucN0 = max(0.0 , aNnucN0)/alpha1;
