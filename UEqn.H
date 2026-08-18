    limitedAlpha1 = min(max(alpha1, scalar(1e-3)), scalar(1));
    taus = (((rho1 * rho2) / rho) * (scalar(1) - limitedAlpha1) * limitedAlpha1) * sqr(Vs);
    surfaceSymmTensorField tausf = fvc::interpolate(taus);
    surfaceScalarField alphaDownwind = fvc::interpolate(alpha1, -fvc::flux(Vs), "upwind");
    surfaceScalarField tauMask = pos0(scalar(0.999) - alphaDownwind);
    surfaceSymmTensorField maskedTausf = tausf * tauMask;
    volVectorField divTaus = fvc::surfaceIntegrate(mesh.Sf() & maskedTausf);

    Info<< "taus magnitude min: " 
    << min(mag(taus)).value() << " , max: "
    << max(mag(taus)).value() << nl << endl;
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      - fvm::Sp(fvc::ddt(rho) + fvc::div(rhoPhi), U)
      + divTaus
      + turbulence->divDevRhoReff(rho, U)
    );

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    interface.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );
    }
