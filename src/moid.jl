"""
    moid(
        saxisA, eccenA, argpeA, omegaA, incliA,
        saxisB, eccenB, argpeB, omegaB, incliB
    ) -> Float64

Calculate MOID (Minimum Orbit Intersection Distance) [AU] between two bodies
whose orbital elements are given. Angles are assumed in [deg].
"""
function moid(
    saxisA::Float64, eccenA::Float64, argpeA::Float64, omegaA::Float64, incliA::Float64,
    saxisB::Float64, eccenB::Float64, argpeB::Float64, omegaB::Float64, incliB::Float64
)::Float64

    # ---- parameters of the program ----
    cstep     = 0.12       # scanning step of true anomaly/longitude in radians
    stepini   = 0.07       # initial step of first tuning in radians
    steptresh = 1e-5       # final step of first tuning in radians (to choose the MOID)
    stepmin   = 1e-14      # threshold step of second tuning in radians

    # ---- constants ----
    pi     = 3.141592653589793
    twopi  = 2.0 * pi
    degrad = pi / 180.0

    # ---- Convert angles from [deg] to [rad] ----
    argpeA  = argpeA  * degrad
    omegaA  = omegaA  * degrad
    incliA  = incliA  * degrad
    argpeB  = argpeB  * degrad
    omegaB  = omegaB  * degrad
    incliB  = incliB  * degrad

    # ---- Transition matrix (c11...c33) for body A ----
    c11 = cos(omegaA)*cos(argpeA) - sin(omegaA)*cos(incliA)*sin(argpeA)
    c12 = sin(omegaA)*cos(argpeA) + cos(omegaA)*cos(incliA)*sin(argpeA)
    c13 = sin(incliA)*sin(argpeA)

    c21 = -cos(omegaA)*sin(argpeA) - sin(omegaA)*cos(incliA)*cos(argpeA)
    c22 = -sin(omegaA)*sin(argpeA) + cos(omegaA)*cos(incliA)*cos(argpeA)
    c23 = sin(incliA)*cos(argpeA)

    c31 = sin(incliA)*sin(omegaA)
    c32 = -sin(incliA)*cos(omegaA)
    c33 = cos(incliA)

    # ---- Calculate new Euler angles for body B using the transition matrix (z1n,z2n,z3n etc.) ----
    sintmpi = sin(incliB)
    costmpi = cos(incliB)
    costmpo = cos(omegaB)
    sintmpo = sin(omegaB)
    costmpa = cos(argpeB)
    sintmpa = sin(argpeB)

    x1 = costmpo*costmpa - sintmpo*costmpi*sintmpa
    x2 = sintmpo*costmpa + costmpo*costmpi*sintmpa
    x3 = sintmpi*sintmpa

    y1 = -costmpo*sintmpa - sintmpo*costmpi*costmpa
    y2 = -sintmpo*sintmpa + costmpo*costmpi*costmpa
    y3 = sintmpi*costmpa

    z1 = sintmpi*sintmpo
    z2 = -sintmpi*costmpo
    z3 = costmpi

    z1n = c11*z1 + c12*z2 + c13*z3
    z2n = c21*z1 + c22*z2 + c23*z3
    z3n = c31*z1 + c32*z2 + c33*z3

    y3n = c31*y1 + c32*y2 + c33*y3
    x3n = c31*x1 + c32*x2 + c33*x3

    # ---- Use atan(y, x) instead of atan2(y, x) ----
    #  (atan2(y, x) 相当の2引数版関数は Julia では atan(y, x) )
    incliB = atan(sqrt(z1n*z1n + z2n*z2n), z3n)      # was atan2( ..., z3n )
    omegaB = -atan(z1n, -z2n)                       # was -atan2(z1n, -z2n)
    argpeB = -atan(x3n, y3n)                        # was -atan2(x3n, y3n)

    # ---- Helpful precalculated values ----
    costmpo = cos(omegaB)
    sintmpo = sin(omegaB)
    sintmpi = sin(incliB)
    costmpi = z3n   # = cos(incliB) と同じ
    sint    = sintmpo * costmpi
    cost    = costmpo * costmpi
    radA    = saxisA * (1.0 - eccenA*eccenA)
    radB    = saxisB * (1.0 - eccenB*eccenB)

    # ---- Prepare arrays ----
    rAt       = zeros(3)
    rBt       = zeros(3)
    Axt       = zeros(3)
    Ayt       = zeros(3)
    Bxt       = zeros(3)
    Byt       = zeros(3)
    Bzt       = zeros(3)
    tmpmoid   = fill(1e6, 10)
    tmptrueB  = fill(0.0, 10)
    tmplongit = fill(0.0, 10)

    # ---- SCANNING PHASE ----
    # ここで使用する変数を事前に宣言・初期化しておく
    trueB_o   = 0.0
    longit_o  = 0.0
    dist_oo   = 0.0
    dist_o    = 1e6

    trueB     = -2.0 * cstep
    moid      = 1e6

    # 初期化
    tmpmoid[1] = 1e6
    tmpmoid[2] = 1e6
    tmpmoid[3] = 1e6
    tmpmoid[4] = 1e6

    # まず初期3点を作るために2点計算する
    for iii in 1:2
        rB = radB / (1.0 + eccenB*cos(trueB))
        sintmp = sin(trueB + argpeB)
        costmp = cos(trueB + argpeB)
        Bz_sq = sintmpi * sintmp
        Bz_sq = Bz_sq * Bz_sq  # square of Z-coordinate for B

        longit = atan(sintmpo*costmp + sintmp*cost,
                      costmpo*costmp - sintmp*sint)

        tmp2  = eccenA*cos(longit)
        rA   = radA / (1.0 + tmp2)
        rA2  = radA / (1.0 - tmp2)
        tmp1 = rB*sqrt(1.0 - Bz_sq)

        if abs(tmp1 - rA) > abs(tmp1 + rA2)
            rA = rA2
            longit -= pi
            tmp1 = tmp1 + rA2
        else
            tmp1 = tmp1 - rA
        end

        dist = rB*rB*Bz_sq + tmp1*tmp1

        if iii == 1
            dist_oo = dist
        else
            dist_o  = dist
            trueB_o = trueB
            longit_o = longit
        end
        trueB += cstep
    end

    # スキャン本体
    nmax = 0
    dist_min = dist_o
    trueBstop = twopi + cstep

    while trueB < trueBstop
        rB = radB / (1.0 + eccenB*cos(trueB))
        sintmp = sin(trueB + argpeB)
        costmp = cos(trueB + argpeB)
        Bz_sq = sintmpi*sintmp
        Bz_sq = Bz_sq * Bz_sq

        longit = atan(sintmpo*costmp + sintmp*cost,
                      costmpo*costmp - sintmp*sint)

        tmp2  = eccenA*cos(longit)
        rA   = radA / (1.0 + tmp2)
        rA2  = radA / (1.0 - tmp2)
        tmp1 = rB*sqrt(1.0 - Bz_sq)

        if abs(tmp1 - rA) > abs(tmp1 + rA2)
            rA = rA2
            longit -= pi
            tmp1 = tmp1 + rA2
        else
            tmp1 = tmp1 - rA
        end

        dist = rB*rB*Bz_sq + tmp1*tmp1

        # 最小点判定
        if (dist_o <= dist) && (dist_o <= dist_oo)
            nmax += 1
            tmptrueB[nmax]  = trueB_o
            tmplongit[nmax] = longit_o
            tmpmoid[nmax]   = dist_o
        end

        if dist_min > dist
            dist_min = dist
        end

        dist_oo = dist_o
        trueB_o = trueB
        longit_o = longit
        dist_o = dist

        trueB += cstep
    end

    # ---- WATER PROCEDURE ----
    if nmax < 2
        nmax = 4
        for iii in 1:4
            tmptrueB[iii] = (0.25 + 0.5*iii)*pi
            sintmp = sin(tmptrueB[iii] + argpeB)
            costmp = cos(tmptrueB[iii] + argpeB)
            tmplongit[iii] = atan(sintmpo*costmp + sintmp*cost,
                                  costmpo*costmp - sintmp*sint)
            tmpmoid[iii] = 1e6
        end
    end

    # ---- PARALLEL TUNING ----
    for jjj in 1:(nmax+1)
        if jjj <= nmax
            moid_val   = tmpmoid[jjj]
            trueB_m    = tmptrueB[jjj]
            longit_m   = tmplongit[jjj]
            step       = stepini
            threshold  = steptresh
        else
            # nmax+1 回目 → これまでの検出点のうち最小のものを微調整
            if nmax == 2
                if abs(tmpmoid[1] - tmpmoid[2]) < 1e-4
                    # Water procedure に戻す
                    nmax = 1
                    break
                else
                    moid_val = min(tmpmoid[1], tmpmoid[2])
                    idx = tmpmoid[1] < tmpmoid[2] ? 1 : 2
                    trueB_m  = tmptrueB[idx]
                    longit_m = tmplongit[idx]
                end
            else
                moid_val = 1e6
                trueB_m  = 0.0
                longit_m = 0.0
                for iii in 1:(nmax-1)
                    if tmpmoid[iii] < moid_val
                        moid_val   = tmpmoid[iii]
                        trueB_m    = tmptrueB[iii]
                        longit_m   = tmplongit[iii]
                    end
                end
            end
            step      = 2.0*stepini
            threshold = stepmin
        end

        # セットアップ
        rBt[2] = radB / (1.0 + eccenB*cos(trueB_m))
        sintmp = sin(trueB_m + argpeB)
        costmp = cos(trueB_m + argpeB)
        Bxt[2] = costmpo*costmp - sintmp*sint
        Byt[2] = sintmpo*costmp + sintmp*cost
        Bzt[2] = sintmpi*sintmp

        rAt[2] = radA / (1.0 + eccenA*cos(longit_m))
        Axt[2] = cos(longit_m)
        Ayt[2] = sin(longit_m)

        aleft  = true
        aright = true
        bleft  = true
        bright = true

        moidtmp = moid_val

        while step >= threshold
            lpoints = 0
            j1min = 1; j1max = 3
            i1min = 1; i1max = 3

            calc1 = false
            calc2 = false
            calc3 = false
            calc4 = false

            # bleft
            if bleft
                rBt[1] = radB / (1.0 + eccenB*cos(trueB_m - step))
                sintmp = sin(trueB_m - step + argpeB)
                costmp = cos(trueB_m - step + argpeB)
                Bxt[1] = costmpo*costmp - sintmp*sint
                Byt[1] = sintmpo*costmp + sintmp*cost
                Bzt[1] = sintmpi*sintmp
                lpoints += 1
            end

            # bright
            if bright
                rBt[3] = radB / (1.0 + eccenB*cos(trueB_m + step))
                sintmp = sin(trueB_m + step + argpeB)
                costmp = cos(trueB_m + step + argpeB)
                Bxt[3] = costmpo*costmp - sintmp*sint
                Byt[3] = sintmpo*costmp + sintmp*cost
                Bzt[3] = sintmpi*sintmp
                lpoints += 1
            end

            # aleft
            if aleft
                rAt[1] = radA / (1.0 + eccenA*cos(longit_m - step))
                Axt[1] = cos(longit_m - step)
                Ayt[1] = sin(longit_m - step)
                lpoints += 1
            end

            # aright
            if aright
                rAt[3] = radA / (1.0 + eccenA*cos(longit_m + step))
                Axt[3] = cos(longit_m + step)
                Ayt[3] = sin(longit_m + step)
                lpoints += 1
            end

            j1_t = 2
            i1_t = 2

            # lpoints に応じて i1min,i1max, j1min, j1max, calc1..4 の ON/OFF を決める
            if lpoints == 1
                if aleft;  i1max = 1; end
                if aright; i1min = 3; end
                if bleft;  j1max = 1; end
                if bright; j1min = 3; end
            elseif lpoints == 2
                if aleft && bright;  calc1 = true; end
                if aleft && bleft;   calc2 = true; end
                if aright && bright; calc3 = true; end
                if aright && bleft;  calc4 = true; end
            end

            moid_loc = moidtmp
            for j1 in j1min:j1max
                for i1 in i1min:i1max
                    # lpoints==2 の場合、不要な組み合わせはスキップ
                    if lpoints == 2
                        if i1 != 1
                            if ((j1 != 3) && calc1) || ((j1 != 1) && calc2)
                                continue
                            end
                        end
                        if i1 != 3
                            if ((j1 != 3) && calc3) || ((j1 != 1) && calc4)
                                continue
                            end
                        end
                    end
                    # 中央(2,2) は既に計算済みとしてスキップ
                    if i1 == 2 && j1 == 2
                        continue
                    end

                    Dx = rBt[j1]*Bxt[j1] - rAt[i1]*Axt[i1]
                    Dy = rBt[j1]*Byt[j1] - rAt[i1]*Ayt[i1]
                    Dz = rBt[j1]*Bzt[j1]
                    dist = Dx*Dx + Dy*Dy + Dz*Dz

                    if dist < moid_loc
                        moid_loc = dist
                        j1_t = j1
                        i1_t = i1
                    end
                end
            end

            if (j1_t != 2) || (i1_t != 2)
                # best が中央(2,2)以外
                aleft  = false
                aright = false
                bleft  = false
                bright = false

                # 選ばれた i1_t, j1_t によって A/B を一歩ずらす
                if i1_t != 2
                    if i1_t == 1
                        aleft = true
                        longit_m -= step
                        rAt[3] = rAt[2]; Axt[3] = Axt[2]; Ayt[3] = Ayt[2]
                        rAt[2] = rAt[1]; Axt[2] = Axt[1]; Ayt[2] = Ayt[1]
                    else
                        aright = true
                        longit_m += step
                        rAt[1] = rAt[2]; Axt[1] = Axt[2]; Ayt[1] = Ayt[2]
                        rAt[2] = rAt[3]; Axt[2] = Axt[3]; Ayt[2] = Ayt[3]
                    end
                end

                if j1_t != 2
                    if j1_t == 1
                        bleft = true
                        trueB_m -= step
                        rBt[3] = rBt[2]; Bxt[3] = Bxt[2]; Byt[3] = Byt[2]; Bzt[3] = Bzt[2]
                        rBt[2] = rBt[1]; Bxt[2] = Bxt[1]; Byt[2] = Byt[1]; Bzt[2] = Bzt[1]
                    else
                        bright = true
                        trueB_m += step
                        rBt[1] = rBt[2]; Bxt[1] = Bxt[2]; Byt[1] = Byt[2]; Bzt[1] = Bzt[2]
                        rBt[2] = rBt[3]; Bxt[2] = Bxt[3]; Byt[2] = Byt[3]; Bzt[2] = Bzt[3]
                    end
                end
                moidtmp = moid_loc
            else
                # best が中央(2,2) → stepを小さく
                aleft  = true
                aright = true
                bleft  = true
                bright = true
                step   *= 0.15
            end
        end

        if jjj <= nmax
            tmpmoid[jjj]   = moidtmp
            tmptrueB[jjj]  = trueB_m
            tmplongit[jjj] = longit_m
        end
    end

    # 最終的に tmpmoid の最小値をとって sqrt
    bestsq = 1e6
    for k in 1:length(tmpmoid)
        if tmpmoid[k] < bestsq
            bestsq = tmpmoid[k]
        end
    end

    return sqrt(bestsq)
end
