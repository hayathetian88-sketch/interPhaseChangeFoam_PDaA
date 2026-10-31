{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phir("phir", phic*interface.nHatf());

    // 体積平均速度 j のフラックス。
    // pslipf と phiVs は前ステップで整合が取れたものをそのまま使う。
    surfaceScalarField phij("phij", phi + pslipf);
    volScalarField divUj
    (
        "divUj",
        mesh.moving()
      ? fvc::div(phij + mesh.phi())
      : fvc::div(phij)
    );

    surfaceScalarField phiVsAlpha2("phiVsAlpha2", phiVs*fvc::interpolate(alpha2));

    // alpha1 をスリップで運ぶ実効フラックス（液相は j に対して -alpha2*Vs で動く）
    surfaceScalarField phiSlipAdv("phiSlipAdv", -phiVsAlpha2);

    Pair<tmp<volScalarField>> vDotAlphal = mixture->vDotAlphal();
    vDotcAlphal = vDotAlphal[0]();
    vDotvAlphal = vDotAlphal[1]();
    vDotvmcAlphal = vDotvAlphal - vDotcAlphal;

    {
        const scalar dt = runTime.deltaTValue();
        auto co = [&](const surfaceScalarField& f)
        {
            return 0.5*gMax(fvc::surfaceSum(mag(f))().primitiveField()/mesh.V().field())*dt;
        };
        Info<< "Co(phi) = " << co(phi)
            << "  Co(phij) = " << co(phij)
            << "  Co(slip) = " << co(phiVsAlpha2)
            << "  max|Vs| = " << gMax(mag(Vs.primitiveField()))
            << "  max|divUj| = " << gMax(mag(divUj.primitiveField())) << endl;
    }

    tmp<surfaceScalarField> talphaPhi;

    if (MULESCorr)
    {
        // 予測子: 対流・スリップとも upwind（有界な予測解を作る）
        fvScalarMatrix alpha1Eqn
        (
            fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phij,
                upwind<scalar>(mesh, phij)
            ).fvmDiv(phij, alpha1)
          - fvm::Sp(divUj, alpha1)
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phiSlipAdv,
                upwind<scalar>(mesh, phiSlipAdv)
            ).fvmDiv(phiSlipAdv, alpha1)
         ==
            fvm::Sp(vDotvmcAlphal, alpha1)
          + vDotcAlphal
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
        // 高次流束: j による対流 + 圧縮 - スリップ
        tmp<surfaceScalarField> talphaPhiCorr
        (
            fvc::flux(phij, alpha1, alphaScheme)
          + fvc::flux(-fvc::flux(-phir, alpha2, alpharScheme), alpha1, alpharScheme)
          - phiVsAlpha2*fvc::interpolate(alpha1)
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
                    divUj*(alpha10 - alpha100)
                  - vDotvmcAlphal*alpha10
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
                phij,
                talphaPhiCorr.ref(),
                vDotvmcAlphal,
                (divUj*min(alpha1, scalar(1)) + vDotcAlphal)(),
                oneField(),
                zeroField()
            );

            talphaPhi = talphaPhiCorr;
        }

        Info<< "Liquid phase volume fraction = "
            << alpha1.weightedAverage(mesh.V()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;

        // クリップ前後の総量（[0, 1] の範囲外のみクリップ）
        const scalar a1Pre = fvc::domainIntegrate(alpha1).value();

        alpha1 = min(max(alpha1, scalar(0)), scalar(1));
        alpha2 = 1.0 - alpha1;

        Info<< "alpha1 clip change (aCorr " << aCorr << ") = "
            << fvc::domainIntegrate(alpha1).value() - a1Pre << endl;
    }

    // 質量流束: rho1*(alpha1 u1) + rho2*(j - alpha1 u1)
    rhoPhi = talphaPhi()*(rho1 - rho2) + phij*rho2;
    alphaPhi10 = talphaPhi();

    Info<< "Liquid phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;

    Info<< "divUj Min = " << min(divUj).value()
        << "  Max = " << max(divUj).value() << nl
        << "vDotvmcAlphal Min = " << min(vDotvmcAlphal).value()
        << "  Max = " << max(vDotvmcAlphal).value() << nl
        << "vDotcAlphal     Min = " << min(vDotcAlphal).value()
        << "  Max = " << max(vDotcAlphal).value() << nl
        << endl;
}
