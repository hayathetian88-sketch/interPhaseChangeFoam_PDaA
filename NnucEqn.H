surfaceScalarField phiVs = fvc::interpolate(Vs) & mesh.Sf();
surfaceScalarField phiNnuc = phi + phiVs;

fvScalarMatrix NnucEqn
(
    fv::EulerDdtScheme<scalar>(mesh).fvmDdt(NnucN0)
    + fv::gaussConvectionScheme<scalar>
      (
        mesh,
        phiNnuc,
        upwind<scalar>(mesh, phiNnuc)
      ).fvmDiv(phiNnuc, NnucN0)
);


NnucEqn.relax();
NnucEqn.solve();

NnucN0 = max(0.0 , NnucN0);
