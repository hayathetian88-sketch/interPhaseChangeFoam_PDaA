{
    // スリップ流束 pslipf = Cslipf*phiVs を
    //   pslipf = Cslipf*flux(HbyVs) - Cslipf*DVsf*snGrad(p_rgh)*|Sf|
    // と分解し、圧力勾配部分を p_rgh のラプラシアン係数 Dslipf に入れて陰的に解く
    surfaceScalarField Cslipf
    (
        "Cslipf",
        fvc::interpolate((rho1 - rho2)/rho*limitedAlpha1*(1.0 - limitedAlpha1))
    );
    surfaceScalarField phiHbyVs("phiHbyVs", Cslipf*fvc::flux(HbyVs));
    surfaceScalarField Dslipf("Dslipf", Cslipf*DVsf);

    forAll(Cslipf.boundaryField(), patchi)
    {
        if (!mesh.boundary()[patchi].coupled())
        {
            phiHbyVs.boundaryFieldRef()[patchi] = 0.0;
            Dslipf.boundaryFieldRef()[patchi] = 0.0;
        }
    }

    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi)
    );
    adjustPhi(phiHbyA, U, p_rgh);

    surfaceScalarField phig
    (
        (
            interface.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf);

    Pair<tmp<volScalarField>> vDotP = mixture->vDotP();
    vDotcP = vDotP[0]();
    vDotvP = vDotP[1]();

    while (pimple.correctNonOrthogonal())
    {
        volScalarField phaseChangeCoeff("phaseChangeCoeff", vDotvP - vDotcP);
        volScalarField phaseChangeSource("phaseChangeSource", phaseChangeCoeff * p_rgh);
        Info<< "PhaseChange Source (Coeff * p_rgh) min: " 
            << min(phaseChangeSource).value() << " , max: " 
            << max(phaseChangeSource).value() << nl << endl;

        fvScalarMatrix p_rghEqn
        (
            fvc::div(phiHbyA + phiHbyVs) - fvm::laplacian(rAUf + Dslipf, p_rgh)
          - (vDotvP - vDotcP)*(mixture->pSat() - rho*gh)
          + fvm::Sp(vDotvP - vDotcP, p_rgh)
        );

        p_rghEqn.setReference(pRefCell, pRefValue);

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            // -grad(p_rgh).Sf の面値（非直交補正込み）
            surfaceScalarField gradpFlux
            (
                "gradpFlux",
                p_rghEqn.flux()/(rAUf + Dslipf)
            );

            phi = phiHbyA + rAUf*gradpFlux;

            U = HbyA + rAU*fvc::reconstruct((phig + rAUf*gradpFlux)/rAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);

            // 新しい p_rgh でスリップ流束を更新（div(phi + pslipf) が相変化源と厳密に整合）
            phiVs = fvc::flux(HbyVs) + DVsf*gradpFlux;
            forAll(phiVs.boundaryField(), patchi)
            {
                if (!mesh.boundary()[patchi].coupled())
                {
                    phiVs.boundaryFieldRef()[patchi] = 0.0;
                }
            }
            pslipf = phiHbyVs + Dslipf*gradpFlux;

            Vs = HbyVs + DVs*fvc::reconstruct(gradpFlux);
            {
                const dimensionedScalar VsMax("VsMax", dimVelocity, mixture->kw().value());
                Vs *= min
                (
                    VsMax/(mag(Vs) + dimensionedScalar("vSmall", dimVelocity, VSMALL)),
                    dimensionedScalar("one", dimless, 1.0)
                );
            }
            Vs.correctBoundaryConditions();

            pslipdiv = fvc::div(pslipf);
            Info<< "pslipdiv magnitude min: "
                << min(pslipdiv).value() << " , max: "
                << max(pslipdiv).value() << nl << endl;
        }
    }

    p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }
}
