{
    limitedAlpha1 = min(max(alpha1, scalar(1e-5)), scalar(1));

    tmp<volScalarField> tRb = mixture->rRb(limitedAlpha1);
    rRb = tRb();
    Info<< "rRb magnitude min: " 
        << min(rRb).value() << " , max: "
        << max(rRb).value() << nl << endl;

    const dimensionedScalar& nu = mixture->nuw(); 
    const dimensionedScalar& kd = mixture->kd();
    const dimensionedScalar& kw = mixture->kw();
    const dimensionedScalar vSmall("vSmall", dimVelocity, VSMALL);

    // 1. レイノルズ数 Re_b と抗力関数 kd_dynamic
    volScalarField Re_b("Re_b", (2.0 / rRb) * mag(Vs) / nu);
    volScalarField Re_eff("Re_eff", max(Re_b, scalar(1e-5)));

    volScalarField kd_SN
    (
        "kd_SN",
        24.0 * (1.0 + 0.15 * pow(Re_eff, 0.687))
    );
    volScalarField kd_high("kd_high", 0.44 * Re_eff);

    volScalarField kd_dynamic
    (
        "kd_dynamic",
        max(kd_SN, kd_high) * pow(limitedAlpha1, -0.5)
    );

    // 2. 緩和時間逆数 time_coeff [1/s]
    time_coeff =(3.0/8.0)*(kd_dynamic*nu)*(rRb*rRb);

    //    抗力 f = time_coeff(|Vs|)*Vs の Newton 線形化用の指数 nDrag = dln(Cd Re)/dln Re
    //    （Schiller-Naumann 側: 0.15*0.687*Re^0.687/(1+0.15 Re^0.687)、Newton 側(0.44Re): 1）
    //    f ~ (1+nDrag)*time_coeff*Vs - nDrag*time_coeff*Vs*  （Vs* は前反復値。収束すれば元の抗力と一致）
    volScalarField nDrag
    (
        "nDrag",
        pos0(kd_high - kd_SN)
      + neg(kd_high - kd_SN)
       *(0.15*0.687*pow(Re_eff, 0.687))/(1.0 + 0.15*pow(Re_eff, 0.687))
    );

    Info<< "time_coeff magnitude min: " 
        << min(time_coeff).value() << " , max: "
        << max(time_coeff).value() << nl << endl;

    // 3. 圧力駆動力の係数 cGradp = 2 * kd / rho1
    const dimensionedScalar cGradp("cGradp", kd.value() * 2.0 / rho1);

    // セル中心の圧力勾配（ガウス積分）
    volVectorField gradp("gradp", fvc::grad(p_rgh));

    Vt = -cGradp * gradp / (time_coeff + dimensionedScalar("tcSmall", dimless/dimTime, SMALL));
    Info<< "Vt magnitude min: " 
        << min(mag(Vt)).value() << " , max: "
        << max(mag(Vt)).value() << nl << endl;

    // 4. 気相速度フラックス phiv
    surfaceScalarField phiv("phiv", phi + fvc::interpolate((rho1/rho)*limitedAlpha1)*phiVs);
    forAll(phiv.boundaryField(), patchi)
    {
        if (isA<wallFvPatch>(mesh.boundary()[patchi]))
        {
            phiv.boundaryFieldRef()[patchi] = 0.0;
        }
    }

    // 液相速度 ~ 混合速度 (u_l = U - Y2*Vs で Y2 << 1 の分散気泡近似と整合)
    volTensorField gradVl("gradVl", fvc::grad(U));

    // 有限気泡径の補正: 気泡は自身の直径より細かい速度勾配を感じない
    //   |grad U_eff| <= |U|/(2 Rb)  ->  fSize = min(|U| rRb/2 / |grad U|, 1)
    //   小さな気泡・緩やかな勾配では fSize = 1 で元のモデルと同一
    fSize = min
    (
        mag(U)*rRb/2.0
       /(mag(gradVl) + dimensionedScalar(dimless/dimTime, SMALL)),
        scalar(1)
    );

    // 5. スリップ速度輸送方程式
    //    圧力勾配項は行列に含めずに組み立て、Vs = HbyVs - DVs*grad(p_rgh) と分解する。
    //    圧力駆動分は pEqn でラプラシアンとして陰的に解く（陽的に入れると kd が大きいとき
    //    p <-> Vs の反復が発散するため）
    fvVectorMatrix VsEqn
    (
        fvm::ddt(Vs)
      + fvm::div(phiv, Vs)
      + fvm::SuSp(-fvc::div(phiv), Vs)
      + fvm::Sp((1.0 + nDrag)*time_coeff, Vs)
     ==
      - fSize*(Vs & gradVl)
      + nDrag*time_coeff*Vs
    );

    VsEqn.relax();
    solve(VsEqn == -cGradp*gradp);

    // 6. 圧力分解用の係数（キャップ前の Vs で H を評価して分解を厳密に保つ）
    {
        volScalarField rAVs("rAVs", 1.0/VsEqn.A());
        HbyVs = rAVs*VsEqn.H();
        HbyVs.correctBoundaryConditions();
        DVs = cGradp*rAVs;
        DVsf = fvc::interpolate(DVs);
    }

    // --- 診断: Vs を駆動・抑制している各項の大きさ [m/s2] ---
    {
        volScalarField fP("fP_gradp", cGradp*mag(gradp));
        volScalarField fD("fD_drag", time_coeff*mag(Vs));
        volScalarField fS("fS_strain", fSize*mag(Vs & gradVl));
        volScalarField fC("fC_divphiv", max(fvc::div(phiv), dimensionedScalar(dimless/dimTime, Zero))*mag(Vs));
        volScalarField strainRatio("strainRatio", fSize*mag(gradVl)/time_coeff);

        const scalarField magVs(mag(Vs.primitiveField()));
        const scalar maxVs = gMax(magVs);
        forAll(magVs, celli)
        {
            if (magVs[celli] >= maxVs)
            {
                Pout<< "maxVs cell " << celli << " at " << mesh.C()[celli]
                    << " |Vs|=" << magVs[celli]
                    << " alpha1=" << alpha1[celli]
                    << " Rb=" << 1.0/rRb[celli]
                    << " Re=" << Re_b[celli]
                    << " time_coeff=" << time_coeff[celli]
                    << " nDrag=" << nDrag[celli]
                    << " | gradp=" << fP[celli]
                    << " drag=" << fD[celli]
                    << " strain=" << fS[celli]
                    << " fSize=" << fSize[celli]
                    << " divphiv=" << fC[celli]
                    << " strainRatio=" << strainRatio[celli]
                    << endl;
                break;
            }
        }

        if (runTime.writeTime() && pimple.finalIter())
        {
            fP.write();
            fD.write();
            fS.write();
            fC.write();
            strainRatio.write();
            nDrag.write();
        }
    }

    // 7. 面フラックス phiVs（Rhie-Chow 形式。pEqn 後に新しい p_rgh で更新し直される）
    phiVs = fvc::flux(HbyVs) - DVsf*fvc::snGrad(p_rgh)*mesh.magSf();

    forAll(phiVs.boundaryField(), patchi)
    {
        if (!mesh.boundary()[patchi].coupled())
        {
            phiVs.boundaryFieldRef()[patchi] = 0.0;
        }
    }

    // 8. セル中心 Vs の上限キャップ（phiVs には適用しない）
    const dimensionedScalar VsMax("VsMax", dimVelocity, kw.value());
    Vs *= min(VsMax/(mag(Vs) + vSmall), dimensionedScalar("one", dimless, 1.0));
    Vs.correctBoundaryConditions();

    Info<< "Vs magnitude min: " 
        << min(mag(Vs)).value() << " , max: "
        << max(mag(Vs)).value() << nl << endl;

    Info<< "VvjEqn end" << endl;
}