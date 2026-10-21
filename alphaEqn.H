{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phir("phir", phic*interface.nHatf());
    surfaceScalarField phiVs("phiVs", fvc::flux(Vs));
    forAll(phiVs.boundaryField(), patchi)
    {
        if (!mesh.boundary()[patchi].coupled())
        {
            phiVs.boundaryFieldRef()[patchi] = 0.0;
        }
    }
    surfaceScalarField phiVsAlpha2("phiVsAlpha2", phiVs*fvc::interpolate(alpha2));

    slipDensityCorr = alpha1*alpha2*((rho1-rho2)/rho)*(fvc::grad(alpha1) & Vs);

    Pair<tmp<volScalarField>> vDotAlphal = mixture->vDotAlphal();
    vDotcAlphal = vDotAlphal[0]();
    vDotvAlphal = vDotAlphal[1]();
    vDotvmcAlphal = vDotvAlphal - vDotcAlphal;

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
          - fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phiVsAlpha2,
                linear<scalar>(mesh)
            ).fvmDiv(phiVsAlpha2, alpha1)
         ==
            fvm::Sp(vDotvmcAlphal, alpha1)
          + vDotcAlphal
          - slipDensityCorr
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
        tmp<surfaceScalarField> talphaPhiCorr
        (
            fvc::flux(phi, alpha1, alphaScheme)
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
                    divU*(alpha10 - alpha100)
                  - vDotvmcAlphal*alpha10
                )(),
                oneField(),
                UniformField<scalar>(1e-3)
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
                (divU*alpha1 + vDotcAlphal-slipDensityCorr)(),
                oneField(),
                UniformField<scalar>(1e-3)
            );

            talphaPhi = talphaPhiCorr;
        }

        Info<< "Liquid phase volume fraction = "
            << alpha1.weightedAverage(mesh.V()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;

        alpha1 = (min(max(alpha1, scalar(1e-3)), scalar(1)));
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
        << endl;
}