{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phir("phir", phic*interface.nHatf());
    surfaceScalarField phiVs("phiVs", fvc::flux(Vs));
    surfaceScalarField phiVsAlpha2("phiVsAlpha2", -fvc::flux(-phiVs, alpha2, alpharScheme));

    volScalarField slipDensityCorr0
    (
        "slipDensityCorr0",
        - alpha1*alpha2*((rho1-rho2)/rho)*(fvc::grad(alpha1) & Vs)
    );

    Pair<tmp<volScalarField>> vDotAlphal = mixture->vDotAlphal();
    volScalarField vDotcAlphal("vDotcAlphal", vDotAlphal[0]());
    volScalarField vDotvAlphal("vDotvAlphal", vDotAlphal[1]());
    volScalarField vDotcAlphal0("vDotcAlphal0", vDotcAlphal);
    volScalarField vDotvmcAlphal0("vDotvmcAlphal0", vDotvAlphal - vDotcAlphal0);

    tmp<surfaceScalarField> talphaPhi;

    if (MULESCorr)
    {
        // ★ divUは判定式から除外
        TotalSourceTerm = vDotvmcAlphal0*alpha1 + vDotcAlphal0 + slipDensityCorr0;

        volScalarField blockLowP("blockLowP", pos(scalar(1e-9) - alpha1) * neg(TotalSourceTerm));
        blockFactor = scalar(1e-3)*(scalar(1) - blockLowP);

        vDotcAlphal = vDotcAlphal0 * blockFactor;
        vDotvmcAlphal = vDotvmcAlphal0 * blockFactor;
        slipDensityCorr = slipDensityCorr0 * blockFactor;
        // ★ divU *= blockFactor; は削除

        fvScalarMatrix alpha1Eqn
        (
            fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phi,
                upwind<scalar>(mesh, phi)
            ).fvmDiv(phi, alpha1)
          - fvm::Sp(divU, alpha1)
          - fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phiVsAlpha2,
                upwind<scalar>(mesh, phiVsAlpha2)
            ).fvmDiv(phiVsAlpha2, alpha1)
         ==
            fvm::Sp(vDotvmcAlphal, alpha1)
          + vDotcAlphal
          + slipDensityCorr
        );

        alpha1Eqn.solve();

        Info<< "slipDensityCorr Min = " << min(slipDensityCorr).value()
            << "  Max = " << max(slipDensityCorr).value() << nl;

        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;

        talphaPhi = alpha1Eqn.flux();
    }

    volScalarField alpha10("alpha10", alpha1);

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        TotalSourceTerm = vDotvmcAlphal0*alpha1 + vDotcAlphal0 + slipDensityCorr0;

        volScalarField blockLowP("blockLowP", pos(scalar(1e-9) - alpha1) * neg(TotalSourceTerm));
        blockFactor = scalar(1) - blockLowP;

        slipDensityCorr = slipDensityCorr0 * blockFactor;
        vDotvmcAlphal = vDotvmcAlphal0 * blockFactor;
        vDotcAlphal   = vDotcAlphal0   * blockFactor;
        // ★ divU *= blockFactor; は削除

        tmp<surfaceScalarField> talphaPhiCorr
        (
            fvc::flux(phi, alpha1, alphaScheme)
          + fvc::flux(-fvc::flux(-phir, alpha2, alpharScheme), alpha1, alpharScheme)
          - fvc::flux(-fvc::flux(-phiVs, alpha2, alpharScheme), alpha1, alpharScheme)
        );

        if (MULESCorr)
        {
            talphaPhiCorr.ref() -= talphaPhi();

            volScalarField alpha100("alpha100", alpha10);
            alpha10 = alpha1;

            MULES::correct
            (
                geometricOneField(),
                alpha1,
                talphaPhi(),
                talphaPhiCorr.ref(),
                vDotvmcAlphal,
                (
                    divU*(alpha10 - alpha100)
                  - vDotvmcAlphal*alpha10
                  + slipDensityCorr
                )(),
                oneField(),
                zeroField()
            );

            if (aCorr == 0)
            {
                talphaPhi.ref() += talphaPhiCorr();
            }
            else
            {
                alpha1 = 0.5*alpha1 + 0.5*alpha10;
                talphaPhi.ref() += 0.5*talphaPhiCorr();
            }
        }
        else
        {
            MULES::explicitSolve
            (
                geometricOneField(),
                alpha1,
                phi,
                talphaPhiCorr.ref(),
                vDotvmcAlphal,
                (divU*alpha1 + vDotcAlphal + slipDensityCorr)(),
                oneField(),
                UniformField<scalar>(0.01)
            );

            talphaPhi = talphaPhiCorr;
        }

        Info<< "Liquid phase volume fraction = "
            << alpha1.weightedAverage(mesh.V()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;

        alpha1 = (min(max(alpha1, scalar(0.001)), scalar(1)));
        alpha2 = 1.0 - alpha1;
    }

    rhoPhi = talphaPhi()*(rho1 - rho2) + phi*rho2;
    alphaPhi10 = talphaPhi();

    Info<< "slipDensityCorr Min = " << min(slipDensityCorr).value()
        << "  Max = " << max(slipDensityCorr).value() << nl
        << "Liquid phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;

    Info<< "vDotvmcAlphal Min = " << min(vDotvmcAlphal).value()
        << "  Max = " << max(vDotvmcAlphal).value() << nl
        << "vDotcAlphal     Min = " << min(vDotcAlphal).value()
        << "  Max = " << max(vDotcAlphal).value() << nl
        << "blockFactor Min = " << min(blockFactor).value()
        << "  Max = " << max(blockFactor).value() << endl;
}