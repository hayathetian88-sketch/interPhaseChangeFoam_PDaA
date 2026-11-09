{
    Info<< "Before alphaEqn: max|Vs| = "
        << gMax(mag(Vs.primitiveField())) << endl;

    Info<< "Before alphaEqn: max|phiVs| = "
        << gMax(mag(phiVs.primitiveField())) << endl;

    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phir("phir", phic*interface.nHatf());

    // 体積平均速度 j のフラックス（div(phij) = 相変化膨張率 D を満たす）
    surfaceScalarField phij("phij", phi + pslipf);
    volScalarField divUj
    (
        "divUj",
        mesh.moving()
      ? fvc::div(phij + mesh.phi())
      : fvc::div(phij)
    );

    // 【修正1】alpha2 の補間も中央差分ではなく、スリップの向きに沿った有界スキーム(alpharScheme)で評価する
    // alpha1 のスリップ流束は -alpha1 * alpha2 * phiVs なので、phir と全く同じ構造を持つ
    surfaceScalarField phiVsAlpha2
    (
        "phiVsAlpha2",
        fvc::flux(phiVs, alpha2, alpharScheme)
    );

    // alpha1 をスリップで運ぶ実効フラックス（液相は j に対して -alpha2*Vs で動く）
    surfaceScalarField phiSlipAdv("phiSlipAdv", -phiVsAlpha2);

    // 【修正2】alpha1 を運ぶ全移流フラックス phiAlphaTotal = phij + phiSlipAdv
    surfaceScalarField phiAlphaTotal("phiAlphaTotal", phij + phiSlipAdv);

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
            << "  max|divUj| = " << gMax(mag(divUj.primitiveField())) << endl;
    }
    Info<< "max|alpha2| = " << max(alpha2).value() << endl;
    Info<< "max|phiVsAlpha2| = " << gMax(mag(phiVsAlpha2.primitiveField())) << endl;
    Info<< "max|Vs| = " << gMax(mag(Vs.primitiveField())) << endl;

    tmp<surfaceScalarField> talphaPhi;

    if (MULESCorr)
    {
        // 【修正3】予測子ではスリップ項の非物理的な発散 div(phiSlipAdv) による吸い出しを防ぐため、
        // 予測段階は体積保存が保証されている phij のみで有界な予測解を作り、
        // スリップ流束は直後の MULES::correct で有界性リミッタ付きで組み込む（標準の phir と同じ扱い）
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
        // 【修正4】高次流束: スリップ項も中央差分 interpolate(alpha1) ではなく、
        // fvc::flux(phiSlipAdv, alpha1, alpharScheme) で有界補間する
        tmp<surfaceScalarField> talphaPhiCorr
        (
            fvc::flux(phij, alpha1, alphaScheme)
          + fvc::flux(-fvc::flux(-phir, alpha2, alpharScheme), alpha1, alpharScheme)
          + fvc::flux(phiSlipAdv, alpha1, alpharScheme)
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
            // 【修正5】explicitSolve 時は MULES が talphaPhiCorr 全体をリミットしてくれる
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