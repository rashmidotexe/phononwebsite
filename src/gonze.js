function determinant3x3(matrix) {
    return matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
}

function invert3x3(matrix) {
    let det = determinant3x3(matrix);
    if (Math.abs(det) < 1e-16) {
        return null;
    }
    let invDet = 1 / det;
    return [
        [
            (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) * invDet,
            (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]) * invDet,
            (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) * invDet,
        ],
        [
            (matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2]) * invDet,
            (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) * invDet,
            (matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]) * invDet,
        ],
        [
            (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) * invDet,
            (matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1]) * invDet,
            (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) * invDet,
        ],
    ];
}

function multiplyRowByMatrix(row, matrix) {
    return [
        row[0] * matrix[0][0] + row[1] * matrix[1][0] + row[2] * matrix[2][0],
        row[0] * matrix[0][1] + row[1] * matrix[1][1] + row[2] * matrix[2][1],
        row[0] * matrix[0][2] + row[1] * matrix[1][2] + row[2] * matrix[2][2],
    ];
}

function multiplyMatrixByVector(matrix, vector) {
    return [
        matrix[0][0] * vector[0] + matrix[0][1] * vector[1] + matrix[0][2] * vector[2],
        matrix[1][0] * vector[0] + matrix[1][1] * vector[1] + matrix[1][2] * vector[2],
        matrix[2][0] * vector[0] + matrix[2][1] * vector[1] + matrix[2][2] * vector[2],
    ];
}

function dot3(left, right) {
    return left[0] * right[0] + left[1] * right[1] + left[2] * right[2];
}

function complexScaleAssign(target, scale) {
    target[0] *= scale;
    target[1] *= scale;
}

function createComplexBlockTensor(natoms) {
    let tensor = [];
    for (let i = 0; i < natoms; i++) {
        let row = [];
        for (let j = 0; j < natoms; j++) {
            let block = [];
            for (let alpha = 0; alpha < 3; alpha++) {
                let components = [];
                for (let beta = 0; beta < 3; beta++) {
                    components.push([0, 0]);
                }
                block.push(components);
            }
            row.push(block);
        }
        tensor.push(row);
    }
    return tensor;
}

function getDielectricPart(q, dielectric) {
    let total = 0;
    for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
            total += q[i] * dielectric[i][j] * q[j];
        }
    }
    return total;
}

function getKKTensor(gVector, qCart, qDirectionCart, dielectric, lambda, tolerance) {
    let qK = [
        gVector[0] + qCart[0],
        gVector[1] + qCart[1],
        gVector[2] + qCart[2],
    ];
    let norm = Math.sqrt(dot3(qK, qK));
    let kk = [
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
    ];

    if (norm < tolerance) {
        if (!qDirectionCart) {
            return kk;
        }
        let dielectricPart = getDielectricPart(qDirectionCart, dielectric);
        if (!(Math.abs(dielectricPart) > 1e-16)) {
            return kk;
        }
        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 3; j++) {
                kk[i][j] = qDirectionCart[i] * qDirectionCart[j] / dielectricPart;
            }
        }
        return kk;
    }

    let dielectricPart = getDielectricPart(qK, dielectric);
    if (!(Math.abs(dielectricPart) > 1e-16)) {
        return kk;
    }
    let l2 = 4 * lambda * lambda;
    let prefactor = Math.exp(-dielectricPart / l2) / dielectricPart;
    for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
            kk[i][j] = qK[i] * qK[j] * prefactor;
        }
    }
    return kk;
}

function getDdTensorPart(gList, qCart, qDirectionCart, dielectric, positionsCar, lambda, tolerance) {
    let natoms = positionsCar.length;
    // sum over G of KK(G + q) exp(2 pi i G.(r_i - r_j)); the phase factorizes into
    // exp(2 pi i G.r_i) exp(-2 pi i G.r_j), so only one cos/sin per atom and G
    let sumReal = new Float64Array(natoms * natoms * 9);
    let sumImag = new Float64Array(natoms * natoms * 9);
    let cosAtom = new Float64Array(natoms);
    let sinAtom = new Float64Array(natoms);

    for (let gIndex = 0; gIndex < gList.length; gIndex++) {
        let gVector = gList[gIndex];
        let kk = getKKTensor(gVector, qCart, qDirectionCart, dielectric, lambda, tolerance);
        if (kk.every((row) => row.every((value) => value === 0))) {
            continue;
        }
        for (let i = 0; i < natoms; i++) {
            let phase = 2 * Math.PI * dot3(positionsCar[i], gVector);
            cosAtom[i] = Math.cos(phase);
            sinAtom[i] = Math.sin(phase);
        }
        for (let i = 0; i < natoms; i++) {
            for (let j = 0; j < natoms; j++) {
                let cosPhase = cosAtom[i] * cosAtom[j] + sinAtom[i] * sinAtom[j];
                let sinPhase = sinAtom[i] * cosAtom[j] - cosAtom[i] * sinAtom[j];
                let offset = (i * natoms + j) * 9;
                for (let alpha = 0; alpha < 3; alpha++) {
                    for (let beta = 0; beta < 3; beta++) {
                        let value = kk[alpha][beta];
                        sumReal[offset + alpha * 3 + beta] += value * cosPhase;
                        sumImag[offset + alpha * 3 + beta] += value * sinPhase;
                    }
                }
            }
        }
    }

    let ddPart = createComplexBlockTensor(natoms);
    for (let i = 0; i < natoms; i++) {
        for (let j = 0; j < natoms; j++) {
            let offset = (i * natoms + j) * 9;
            for (let alpha = 0; alpha < 3; alpha++) {
                for (let beta = 0; beta < 3; beta++) {
                    ddPart[i][j][alpha][beta][0] = sumReal[offset + alpha * 3 + beta];
                    ddPart[i][j][alpha][beta][1] = sumImag[offset + alpha * 3 + beta];
                }
            }
        }
    }
    return ddPart;
}

function multiplyBorns(ddIn, born) {
    let natoms = born.length;
    let dd = createComplexBlockTensor(natoms);

    for (let i = 0; i < natoms; i++) {
        for (let j = 0; j < natoms; j++) {
            for (let alpha = 0; alpha < 3; alpha++) {
                for (let beta = 0; beta < 3; beta++) {
                    let target = dd[i][j][alpha][beta];
                    for (let alphaPrime = 0; alphaPrime < 3; alphaPrime++) {
                        for (let betaPrime = 0; betaPrime < 3; betaPrime++) {
                            let zz = born[i][alphaPrime][alpha] * born[j][betaPrime][beta];
                            target[0] += ddIn[i][j][alphaPrime][betaPrime][0] * zz;
                            target[1] += ddIn[i][j][alphaPrime][betaPrime][1] * zz;
                        }
                    }
                }
            }
        }
    }

    return dd;
}

function getQCart(direction, primitiveLattice) {
    let inverseLattice = invert3x3(primitiveLattice);
    if (!inverseLattice) {
        return null;
    }
    return multiplyMatrixByVector(inverseLattice, direction);
}

export function getGonzeReciprocalCorrection(payload, qpoint, qDirection = null) {
    if (!payload.nac || payload.nac.method !== 'gonze') {
        return null;
    }

    let qCart = getQCart(qpoint, payload.primitive_lattice);
    if (!qCart) {
        return null;
    }
    let qDirectionCart = qDirection ? getQCart(qDirection, payload.primitive_lattice) : null;
    let ddPart = getDdTensorPart(
        payload.nac.g_list,
        qCart,
        qDirectionCart,
        payload.nac.dielectric,
        payload.nac.positions_car,
        payload.nac.lambda,
        payload.nac.q_direction_tolerance || 1e-5
    );
    let dd = multiplyBorns(ddPart, payload.nac.born);
    let natoms = payload.masses.length;

    for (let i = 0; i < natoms; i++) {
        for (let alpha = 0; alpha < 3; alpha++) {
            for (let beta = 0; beta < 3; beta++) {
                dd[i][i][alpha][beta][0] -= payload.nac.dd_q0.real[i][alpha][beta];
                dd[i][i][alpha][beta][1] -= payload.nac.dd_q0.imag[i][alpha][beta];
            }
        }
    }

    for (let i = 0; i < natoms; i++) {
        for (let j = 0; j < natoms; j++) {
            let scale = payload.nac.nac_factor / Math.sqrt(payload.masses[i] * payload.masses[j]);
            for (let alpha = 0; alpha < 3; alpha++) {
                for (let beta = 0; beta < 3; beta++) {
                    complexScaleAssign(dd[i][j][alpha][beta], scale);
                }
            }
        }
    }

    return dd;
}

// ---------------------------------------------------------------------------
// real-space part of the Ewald dipole-dipole sum (Gonze & Lee, PRB 55, 10355
// (1997), eqs. 71-75). phonopy keeps only the damped reciprocal-space sum, with an
// Ewald parameter lambda chosen per cell: that sum alone depends on lambda, the
// rest of the dipole interaction staying in the short-range force constants.
// Adding the real-space remainder gives the full dipole-dipole term, independent
// of lambda, so a mixture can subtract and rebuild all of it.
// ---------------------------------------------------------------------------

function erfc(x) {
    /*
    complementary error function to ~1e-15: Taylor series of erf below 3, Lentz
    continued fraction above
    */
    if (x < 0) {
        return 2 - erfc(-x);
    }
    if (x < 3) {
        let term = x;
        let sum = x;
        let x2 = x * x;
        for (let n = 1; n < 200; n++) {
            term *= -x2 / n;
            let add = term / (2 * n + 1);
            sum += add;
            if (Math.abs(add) < 1e-17 * Math.abs(sum)) {
                break;
            }
        }
        return 1 - (2 / Math.sqrt(Math.PI)) * sum;
    }
    // erfc(x) = exp(-x^2)/sqrt(pi) * 1/(x + (1/2)/(x + 1/(x + (3/2)/(x + ...))))
    let tiny = 1e-300;
    let f = x;
    let c = x;
    let d = 0;
    for (let n = 1; n < 300; n++) {
        let a = n / 2;
        d = x + a * d;
        d = Math.abs(d) < tiny ? tiny : d;
        c = x + a / c;
        c = Math.abs(c) < tiny ? tiny : c;
        d = 1 / d;
        let delta = c * d;
        f *= delta;
        if (Math.abs(delta - 1) < 1e-16) {
            break;
        }
    }
    return Math.exp(-x * x) / Math.sqrt(Math.PI) / f;
}

function getRealSpaceTerms(payload) {
    /*
    the real-space Ewald terms of the payload, cached on it: for every atom pair
    and lattice vector R with d = R + r_j - r_i inside the cutoff, the tensor
    Z_i^T H(d) Z_j / sqrt(m_i m_j) (units of the reciprocal part), and the q = 0
    row sums that keep the acoustic sum rule
    */
    if (payload.nac._realSpace) {
        return payload.nac._realSpace;
    }
    let nac = payload.nac;
    let lattice = payload.primitive_lattice;
    let positions = nac.positions_car;
    let born = nac.born;
    let natoms = positions.length;
    let dielectric = nac.dielectric;
    let inverse = invert3x3(dielectric);
    let sqrtDet = Math.sqrt(determinant3x3(dielectric));
    // lambda multiplies G without 2 pi: the Gaussian width in real space is 2 pi lambda
    let ewald = 2 * Math.PI * nac.lambda;
    // erfc(4.5) = 2e-10, exp(-4.5^2) = 2e-9: terms beyond are below 1e-9 of the sum
    let yMax = 4.5;
    // cutoff in the metric d.eps^-1.d; in cartesian length at most yMax/ewald*sqrt(eps_max)
    let epsMax = Math.max(dielectric[0][0], dielectric[1][1], dielectric[2][2]);
    let reach = (yMax / ewald) * Math.sqrt(epsMax);
    let volume = Math.abs(determinant3x3(lattice));
    let cross = (p, q) => [p[1] * q[2] - p[2] * q[1], p[2] * q[0] - p[0] * q[2], p[0] * q[1] - p[1] * q[0]];
    let maxOffset = 0;
    for (let i = 0; i < natoms; i++) {
        for (let j = 0; j < natoms; j++) {
            let d0 = positions[j][0] - positions[i][0];
            let d1 = positions[j][1] - positions[i][1];
            let d2 = positions[j][2] - positions[i][2];
            maxOffset = Math.max(maxOffset, Math.sqrt(d0 * d0 + d1 * d1 + d2 * d2));
        }
    }
    let span = [0, 1, 2].map((k) => {
        let normal = cross(lattice[(k + 1) % 3], lattice[(k + 2) % 3]);
        let height = volume / Math.sqrt(dot3(normal, normal));
        return Math.ceil((reach + maxOffset) / height);
    });
    let vectors = [];
    let terms = [];
    let scale = nac.unit_conversion;
    let rowSums = [];
    for (let i = 0; i < natoms; i++) {
        rowSums.push(new Float64Array(9));
    }
    let inv = inverse;
    let factor = (ewald * ewald * ewald) / sqrtDet;
    let H = new Float64Array(9);
    let ZH = new Float64Array(9);
    for (let a = -span[0]; a <= span[0]; a++) {
        for (let b = -span[1]; b <= span[1]; b++) {
            for (let c = -span[2]; c <= span[2]; c++) {
                let R0 = a * lattice[0][0] + b * lattice[1][0] + c * lattice[2][0];
                let R1 = a * lattice[0][1] + b * lattice[1][1] + c * lattice[2][1];
                let R2 = a * lattice[0][2] + b * lattice[1][2] + c * lattice[2][2];
                let index = -1;
                for (let i = 0; i < natoms; i++) {
                    for (let j = 0; j < natoms; j++) {
                        let d0 = R0 + positions[j][0] - positions[i][0];
                        let d1 = R1 + positions[j][1] - positions[i][1];
                        let d2 = R2 + positions[j][2] - positions[i][2];
                        let r0 = inv[0][0] * d0 + inv[0][1] * d1 + inv[0][2] * d2;
                        let r1 = inv[1][0] * d0 + inv[1][1] * d1 + inv[1][2] * d2;
                        let r2 = inv[2][0] * d0 + inv[2][1] * d1 + inv[2][2] * d2;
                        let length = Math.sqrt(d0 * r0 + d1 * r1 + d2 * r2);
                        let y = ewald * length;
                        if (length < 1e-8 || y > yMax) {
                            continue;
                        }
                        let gauss = (2 / Math.sqrt(Math.PI)) * Math.exp(-y * y);
                        let tail = erfc(y);
                        let along = (3 * tail / (y * y * y) + gauss * (3 / (y * y) + 2)) / (length * length);
                        let isotropic = tail / (y * y * y) + gauss / (y * y);
                        // H in the cartesian axes, then Z_i^T H Z_j
                        let reduced = [r0, r1, r2];
                        for (let p = 0; p < 3; p++) {
                            for (let r = 0; r < 3; r++) {
                                H[p * 3 + r] = factor * (along * reduced[p] * reduced[r] - isotropic * inv[p][r]);
                            }
                        }
                        let Zi = born[i];
                        let Zj = born[j];
                        // (Z_i^T H)[alpha][r]
                        for (let alpha = 0; alpha < 3; alpha++) {
                            for (let r = 0; r < 3; r++) {
                                ZH[alpha * 3 + r] = Zi[0][alpha] * H[r] + Zi[1][alpha] * H[3 + r] + Zi[2][alpha] * H[6 + r];
                            }
                        }
                        let block = new Float64Array(9);
                        let weight = scale / Math.sqrt(payload.masses[i] * payload.masses[j]);
                        for (let alpha = 0; alpha < 3; alpha++) {
                            for (let beta = 0; beta < 3; beta++) {
                                block[alpha * 3 + beta] = (ZH[alpha * 3] * Zj[0][beta] + ZH[alpha * 3 + 1] * Zj[1][beta] +
                                    ZH[alpha * 3 + 2] * Zj[2][beta]) * weight;
                            }
                        }
                        if (index < 0) {
                            index = vectors.length;
                            vectors.push([R0, R1, R2]);
                        }
                        terms.push({ i: i, j: j, r: index, block: block });
                        // acoustic sum rule: q = 0 sums with the unweighted charges
                        let sumWeight = scale / Math.sqrt(payload.masses[i] * payload.masses[i]);
                        for (let k = 0; k < 9; k++) {
                            rowSums[i][k] += block[k] / weight * sumWeight;
                        }
                    }
                }
            }
        }
    }
    payload.nac._realSpace = { vectors: vectors, terms: terms, rowSums: rowSums };
    return payload.nac._realSpace;
}

export function getGonzeRealSpaceRemainder(payload, qpoint) {
    /*
    the rest of the Ewald dipole-dipole sum beside the reciprocal part of
    getGonzeReciprocalCorrection, for the same payload and lambda:
    -sum_R Z_i^T H(R + r_j - r_i) Z_j exp(2 pi i q.(R + r_j - r_i)), with the q = 0
    row sums taken off the diagonal blocks (acoustic sum rule; the self term of the
    Ewald sum cancels there). Same units and block layout as the reciprocal part.
    */
    let qCart = getQCart(qpoint, payload.primitive_lattice);
    let { vectors, terms, rowSums } = getRealSpaceTerms(payload);
    let positions = payload.nac.positions_car;
    let natoms = positions.length;
    let cosR = new Float64Array(vectors.length);
    let sinR = new Float64Array(vectors.length);
    for (let n = 0; n < vectors.length; n++) {
        let phase = 2 * Math.PI * dot3(qCart, vectors[n]);
        cosR[n] = Math.cos(phase);
        sinR[n] = Math.sin(phase);
    }
    let cosAtom = positions.map((r) => Math.cos(2 * Math.PI * dot3(qCart, r)));
    let sinAtom = positions.map((r) => Math.sin(2 * Math.PI * dot3(qCart, r)));
    let out = createComplexBlockTensor(natoms);
    for (let t = 0; t < terms.length; t++) {
        let { i, j, r, block } = terms[t];
        // exp(2 pi i q.(R + r_j - r_i)) = e(R) e(r_j) / e(r_i)
        let cosA = cosAtom[j] * cosAtom[i] + sinAtom[j] * sinAtom[i];
        let sinA = sinAtom[j] * cosAtom[i] - cosAtom[j] * sinAtom[i];
        let cos = cosR[r] * cosA - sinR[r] * sinA;
        let sin = sinR[r] * cosA + cosR[r] * sinA;
        let target = out[i][j];
        for (let alpha = 0; alpha < 3; alpha++) {
            for (let beta = 0; beta < 3; beta++) {
                let value = block[alpha * 3 + beta];
                target[alpha][beta][0] -= value * cos;
                target[alpha][beta][1] -= value * sin;
            }
        }
    }
    for (let i = 0; i < natoms; i++) {
        let sums = rowSums[i];
        for (let alpha = 0; alpha < 3; alpha++) {
            for (let beta = 0; beta < 3; beta++) {
                // hermitian part, as for phonopy's dd_q0
                out[i][i][alpha][beta][0] += (sums[alpha * 3 + beta] + sums[beta * 3 + alpha]) / 2;
            }
        }
    }
    return out;
}

export function getGonzeDipoleDipole(payload, qpoint, qDirection = null) {
    // the full (lambda-independent) Ewald dipole-dipole term: reciprocal + real space
    let reciprocal = getGonzeReciprocalCorrection(payload, qpoint, qDirection);
    let remainder = getGonzeRealSpaceRemainder(payload, qpoint);
    for (let i = 0; i < reciprocal.length; i++) {
        for (let j = 0; j < reciprocal.length; j++) {
            for (let alpha = 0; alpha < 3; alpha++) {
                for (let beta = 0; beta < 3; beta++) {
                    reciprocal[i][j][alpha][beta][0] += remainder[i][j][alpha][beta][0];
                    reciprocal[i][j][alpha][beta][1] += remainder[i][j][alpha][beta][1];
                }
            }
        }
    }
    return reciprocal;
}

// ---------------------------------------------------------------------------
// the same dipole-dipole term, added straight into a dynamical matrix
// ({real, imag} rows of 3 natoms): for the large cells of the mixing page. Every
// term of the reciprocal sum is K(G + q) = k k^T p with k = G + q (or the q
// direction at G = -q), so with u_i = Z_i^T k exp(2 pi i G.r_i) it is
// p u_i u_j^*: hermitian, only i <= j is summed, and Z_i^T G and the phases of
// the atoms are kept per payload (they do not depend on q).
// ---------------------------------------------------------------------------

function getReciprocalCache(payload) {
    if (payload.nac._reciprocal) {
        return payload.nac._reciprocal;
    }
    let nac = payload.nac;
    let positions = nac.positions_car;
    let born = nac.born;
    let natoms = positions.length;
    let gList = nac.g_list;
    let nG = gList.length;
    // Z_i^T G (3 per atom) and exp(2 pi i G.r_i) (cos, sin) per G and atom
    let zg = new Float64Array(nG * natoms * 3);
    let phase = new Float64Array(nG * natoms * 2);
    for (let g = 0; g < nG; g++) {
        let G = gList[g];
        for (let i = 0; i < natoms; i++) {
            let Z = born[i];
            let o = (g * natoms + i) * 3;
            for (let a = 0; a < 3; a++) {
                zg[o + a] = Z[0][a] * G[0] + Z[1][a] * G[1] + Z[2][a] * G[2];
            }
            let r = positions[i];
            let angle = 2 * Math.PI * (G[0] * r[0] + G[1] * r[1] + G[2] * r[2]);
            phase[(g * natoms + i) * 2] = Math.cos(angle);
            phase[(g * natoms + i) * 2 + 1] = Math.sin(angle);
        }
    }
    nac._reciprocal = { zg: zg, phase: phase };
    return nac._reciprocal;
}

function addReciprocalDipoleDipole(payload, qpoint, qDirection, matrix, factor) {
    // factor x getGonzeReciprocalCorrection(payload, qpoint, qDirection), added to matrix
    let nac = payload.nac;
    let natoms = payload.masses.length;
    let qCart = getQCart(qpoint, payload.primitive_lattice);
    let qDirectionCart = qDirection ? getQCart(qDirection, payload.primitive_lattice) : null;
    let dielectric = nac.dielectric;
    let born = nac.born;
    let tolerance = nac.q_direction_tolerance || 1e-5;
    let l2 = 4 * nac.lambda * nac.lambda;
    let gList = nac.g_list;
    let { zg, phase } = getReciprocalCache(payload);
    // Z_i^T q, and Z_i^T of the q direction
    let zq = new Float64Array(natoms * 3);
    let zd = qDirectionCart ? new Float64Array(natoms * 3) : null;
    for (let i = 0; i < natoms; i++) {
        let Z = born[i];
        for (let a = 0; a < 3; a++) {
            zq[i * 3 + a] = Z[0][a] * qCart[0] + Z[1][a] * qCart[1] + Z[2][a] * qCart[2];
            if (zd) {
                zd[i * 3 + a] = Z[0][a] * qDirectionCart[0] + Z[1][a] * qDirectionCart[1] + Z[2][a] * qDirectionCart[2];
            }
        }
    }
    let size = 3 * natoms;
    let accRe = new Float64Array(size * size);
    let accIm = new Float64Array(size * size);
    let uRe = new Float64Array(size);
    let uIm = new Float64Array(size);
    for (let g = 0; g < gList.length; g++) {
        let G = gList[g];
        let k0 = G[0] + qCart[0];
        let k1 = G[1] + qCart[1];
        let k2 = G[2] + qCart[2];
        let p;
        let useDirection = false;
        if (Math.sqrt(k0 * k0 + k1 * k1 + k2 * k2) < tolerance) {
            if (!qDirectionCart) {
                continue;
            }
            let part = getDielectricPart(qDirectionCart, dielectric);
            if (!(Math.abs(part) > 1e-16)) {
                continue;
            }
            p = 1 / part;
            useDirection = true;
        } else {
            let part = getDielectricPart([k0, k1, k2], dielectric);
            if (!(Math.abs(part) > 1e-16)) {
                continue;
            }
            p = Math.exp(-part / l2) / part;
        }
        if (p === 0) {
            continue;
        }
        let sqrtP = Math.sqrt(Math.abs(p));
        let sign = p < 0 ? -1 : 1;
        for (let i = 0; i < natoms; i++) {
            let c = phase[(g * natoms + i) * 2];
            let s = phase[(g * natoms + i) * 2 + 1];
            for (let a = 0; a < 3; a++) {
                let v = useDirection ? zd[i * 3 + a] : zg[(g * natoms + i) * 3 + a] + zq[i * 3 + a];
                uRe[i * 3 + a] = v * c * sqrtP;
                uIm[i * 3 + a] = v * s * sqrtP;
            }
        }
        // acc_xy += sign u_x conj(u_y), for the blocks i <= j
        for (let i = 0; i < natoms; i++) {
            for (let x = 3 * i; x < 3 * i + 3; x++) {
                let ar = uRe[x] * sign;
                let ai = uIm[x] * sign;
                let row = x * size;
                for (let y = 3 * i; y < size; y++) {
                    let br = uRe[y];
                    let bi = uIm[y];
                    accRe[row + y] += ar * br + ai * bi;
                    accIm[row + y] += ai * br - ar * bi;
                }
            }
        }
    }
    let ddQ0 = nac.dd_q0;
    for (let i = 0; i < natoms; i++) {
        for (let j = i; j < natoms; j++) {
            let scale = factor * nac.nac_factor / Math.sqrt(payload.masses[i] * payload.masses[j]);
            for (let a = 0; a < 3; a++) {
                for (let b = 0; b < 3; b++) {
                    let re = accRe[(3 * i + a) * size + 3 * j + b];
                    let im = accIm[(3 * i + a) * size + 3 * j + b];
                    if (i === j) {
                        re -= ddQ0.real[i][a][b];
                        im -= ddQ0.imag[i][a][b];
                        matrix.real[3 * i + a][3 * j + b] += scale * re;
                        matrix.imag[3 * i + a][3 * j + b] += scale * im;
                    } else {
                        matrix.real[3 * i + a][3 * j + b] += scale * re;
                        matrix.imag[3 * i + a][3 * j + b] += scale * im;
                        // the block j i: the hermitian conjugate
                        matrix.real[3 * j + b][3 * i + a] += scale * re;
                        matrix.imag[3 * j + b][3 * i + a] -= scale * im;
                    }
                }
            }
        }
    }
}

export function addGonzeRealSpaceRemainder(payload, qpoint, matrix, factor) {
    // factor x getGonzeRealSpaceRemainder(payload, qpoint), added to matrix
    let qCart = getQCart(qpoint, payload.primitive_lattice);
    let { vectors, terms, rowSums } = getRealSpaceTerms(payload);
    let positions = payload.nac.positions_car;
    let natoms = positions.length;
    let cosR = new Float64Array(vectors.length);
    let sinR = new Float64Array(vectors.length);
    for (let n = 0; n < vectors.length; n++) {
        let angle = 2 * Math.PI * dot3(qCart, vectors[n]);
        cosR[n] = Math.cos(angle);
        sinR[n] = Math.sin(angle);
    }
    let cosAtom = positions.map((r) => Math.cos(2 * Math.PI * dot3(qCart, r)));
    let sinAtom = positions.map((r) => Math.sin(2 * Math.PI * dot3(qCart, r)));
    for (let t = 0; t < terms.length; t++) {
        let { i, j, r, block } = terms[t];
        let cosA = cosAtom[j] * cosAtom[i] + sinAtom[j] * sinAtom[i];
        let sinA = sinAtom[j] * cosAtom[i] - cosAtom[j] * sinAtom[i];
        let cos = (cosR[r] * cosA - sinR[r] * sinA) * factor;
        let sin = (sinR[r] * cosA + cosR[r] * sinA) * factor;
        for (let alpha = 0; alpha < 3; alpha++) {
            let rowRe = matrix.real[3 * i + alpha];
            let rowIm = matrix.imag[3 * i + alpha];
            for (let beta = 0; beta < 3; beta++) {
                let value = block[alpha * 3 + beta];
                rowRe[3 * j + beta] -= value * cos;
                rowIm[3 * j + beta] -= value * sin;
            }
        }
    }
    for (let i = 0; i < natoms; i++) {
        let sums = rowSums[i];
        for (let alpha = 0; alpha < 3; alpha++) {
            for (let beta = 0; beta < 3; beta++) {
                matrix.real[3 * i + alpha][3 * i + beta] += factor * (sums[alpha * 3 + beta] + sums[beta * 3 + alpha]) / 2;
            }
        }
    }
}

export function addGonzeDipoleDipole(payload, qpoint, qDirection, matrix, factor) {
    // factor x getGonzeDipoleDipole(payload, qpoint, qDirection), added to matrix
    addReciprocalDipoleDipole(payload, qpoint, qDirection, matrix, factor);
    addGonzeRealSpaceRemainder(payload, qpoint, matrix, factor);
}
