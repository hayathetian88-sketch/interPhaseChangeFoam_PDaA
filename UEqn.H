volScalarField cs
    (
        "cs",
        ((rho1*rho2)/rho)*(scalar(1) - limitedAlpha1)*limitedAlpha1
    );

    // (B) V_j * div(c V): 面フラックスで発散を取る（divSchemes の指定が不要）
    volScalarField divCVs("divCVs", fvc::div(fvc::interpolate(cs)*phiVs));

    // (A) c * V_i dV_j/dx_i
    taus = Vs*divCVs + cs*(Vs & fvc::grad(Vs));

    Info<< "taus magnitude min: " 
        << min(mag(taus)).value() << " , max: "
        << max(mag(taus)).value() << nl << endl;

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      - fvm::Sp(fvc::ddt(rho) + fvc::div(rhoPhi), U)
      + taus
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
    