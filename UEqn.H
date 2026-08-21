    limitedAlpha1 = min(max(alpha1, scalar(1e-3)), scalar(1));
    taus = blockFactor * (((rho1 * rho2) / rho) * (scalar(1) - limitedAlpha1) * limitedAlpha1) * sqr(Vs);
    surfaceSymmTensorField tausf = fvc::interpolate(taus);
    volVectorField divTaus = fvc::surfaceIntegrate(mesh.Sf() & tausf);

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
