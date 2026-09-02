{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phir("phir", phic*interface.nHatf());

    // ★ 生の値(ブロック前)を別名で保持する
    Pair<tmp<volScalarField>> vDotAlphal = mixture->vDotAlphal();
    volScalarField vDotcAlphal("vDotcAlphal", vDotAlphal[0]());
    volScalarField vDotvAlphal("vDotvAlphal", vDotAlphal[1]());
    volScalarField vDotcAlphal0("vDotcAlphal0", vDotcAlphal);
    volScalarField vDotvmcAlphal0("vDotvmcAlphal0", vDotvAlphal - vDotcAlphal0);

    volVectorField slipFluxVec("slipFluxVec", (alpha1*alpha2*(rho1/rho))*Vs);
    surfaceScalarField slipFlux("slipFlux", fvc::flux(slipFluxVec));
    volScalarField slipTerm0("slipTerm0", (rho/rho1)*fvc::div(slipFlux));

    // ★ 総ソース量制限のための定数
    dimensionedScalar maxDalpha("maxDalpha", dimless, 0.01);
    dimensionedScalar rDeltaT("rDeltaT", dimless/dimTime, 1.0/runTime.deltaT().value());
    Info<< "maxDalpha = " << maxDalpha.value()
    << ", rDeltaT = " << rDeltaT.value()
    << endl;

    tmp<surfaceScalarField> talphaPhi;

    if (MULESCorr)
    {
        TotalSourceTerm = vDotvmcAlphal0*alpha1 + vDotcAlphal0 + slipTerm0+divU*alpha1;

        // ★ 総量制限のみ(下限ブロックは外す)
        blockFactor = min
        (
            maxDalpha*rDeltaT
          / (mag(TotalSourceTerm) + dimensionedScalar("small", TotalSourceTerm.dimensions(), SMALL)),
            scalar(1)
        )*(1- pos(scalar(0.011) - alpha1) * neg(TotalSourceTerm));

        vDotcAlphal   = vDotcAlphal0   ;
        vDotvmcAlphal = vDotvmcAlphal0 ;
        slipTerm      = slipTerm0      * blockFactor;
        divU *= blockFactor;

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
         ==
            fvm::Sp(vDotvmcAlphal, alpha1)
          + vDotcAlphal
          + slipTerm
        );

        alpha1Eqn.solve();

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
        TotalSourceTerm = vDotvmcAlphal0*alpha1 + vDotcAlphal0 + slipTerm0 + divU*alpha1;

        // ★ 総量制限のみ、毎反復その時点のalpha1で再評価
        blockFactor = min
        (
            maxDalpha*rDeltaT
          / (mag(TotalSourceTerm) + dimensionedScalar("small", TotalSourceTerm.dimensions(), SMALL)),
            scalar(1)
        )*(1- pos(scalar(0.011) - alpha1) * neg(TotalSourceTerm));

        slipTerm      = slipTerm0      * blockFactor;
        vDotvmcAlphal = vDotvmcAlphal0 ;
        vDotcAlphal   = vDotcAlphal0   ;
        divU *= blockFactor;

        tmp<surfaceScalarField> talphaPhiCorr
        (
            fvc::flux
            (
                phi,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, alpha2, alpharScheme),
                alpha1,
                alpharScheme
            )
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
                  + slipTerm
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
                (divU*alpha1 + vDotcAlphal + slipTerm)(),
                oneField(),
                UniformField<scalar>(0.0)
            );

            talphaPhi = talphaPhiCorr;
        }

        Info<< "Liquid phase volume fraction = "
            << alpha1.weightedAverage(mesh.V()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;

        alpha1 = (min(max(alpha1, scalar(0.01)), scalar(1)));
        alpha2 = 1.0 - alpha1;
    }

    rhoPhi = talphaPhi()*(rho1 - rho2) + phi*rho2;

    // Cache alphaPhi
    alphaPhi10 = talphaPhi();

    Info<< "Liquid phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;

    Info
        << "slipTerm Min = " << min(slipTerm).value()
        << "  Max = " << max(slipTerm).value() << nl
        << "vDotvmcAlphal Min = " << min(vDotvmcAlphal).value()
        << "  Max = " << max(vDotvmcAlphal).value() << nl
        << "vDotcAlphal     Min = " << min(vDotcAlphal).value()
        << "  Max = " << max(vDotcAlphal).value() << nl
        << "blockFactor Min = " << min(blockFactor).value()
        << "  Max = " << max(blockFactor).value() << endl;
}