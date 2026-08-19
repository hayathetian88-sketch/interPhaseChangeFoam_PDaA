{
    volScalarField blockHigh
    (
        "blockHigh",
        pos(alpha1 - scalar(0.99)) * pos(slipTerm)
    );

    // alpha1<0.01 かつ slipTermが負(さらに0を割らせる方向) → ブロック
    volScalarField blockLow
    (
        "blockLow",
        pos(scalar(0.01) - alpha1) * neg(slipTerm)
    );
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");
    
    surfaceScalarField phir("phir", phic*interface.nHatf());

    Pair<tmp<volScalarField>> vDotAlphal =
        mixture->vDotAlphal();
    const volScalarField& vDotcAlphal = vDotAlphal[0]();
    const volScalarField& vDotvAlphal = vDotAlphal[1]();
    // const volScalarField vDotvmcAlphal(vDotvAlphal - vDotcAlphal);
    vDotvmcAlphal = vDotvAlphal - vDotcAlphal;

    volVectorField slipFluxVec("slipFluxVec", (alpha1*alpha2/rho)*Vs);
    surfaceScalarField slipFlux("slipFlux", fvc::flux(slipFluxVec));
    slipTerm = rho*fvc::div(slipFlux)*(scalar(1) - blockHigh) * (scalar(1) - blockLow);

    tmp<surfaceScalarField> talphaPhi;

    if (MULESCorr)
    {
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

            // Under-relax the correction for all but the 1st corrector
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
                (divU*alpha1 + vDotcAlphal + slipTerm
                )(),
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
        << "  Max = " << max(vDotcAlphal).value() << endl;
}
